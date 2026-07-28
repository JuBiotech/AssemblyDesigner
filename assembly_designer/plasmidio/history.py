"""History / provenance tracking for DNA assembly workflows (build & plot assembly DAGs)."""

from __future__ import annotations

import re
import statistics
from collections.abc import Iterable, Mapping, MutableMapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Protocol

import networkx as nx
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.figure import Figure

from ._shared import _ADAPTER_2L, _UNS_NUMBER


class PCRResultLike(Protocol):
    """Minimal interface we rely on from a PCR result object."""

    @property
    def product(self) -> SeqRecord:  # pragma: no cover - typing-only
        """PCR product sequence."""
        ...


@dataclass(slots=True)
class HistoryNode:
    """A single node in the assembly DAG."""

    nid: str
    label: str
    kind: str  # "Source" | "PCR" | "GoldenGate" | "Gibson" | "Other"
    length: int
    meta: dict[str, Any] = field(default_factory=dict)
    record: Optional[SeqRecord] = None


@dataclass(slots=True)
class HistoryEdge:
    """A directed labeled edge between two nodes."""

    src: str
    dst: str
    label: str = ""


class AssemblyHistory:
    """Provenance tracker for DNA assembly workflows.

    Tracks sources (parts, TU templates, adapters), intermediate steps
    (PCR, Golden Gate), and final constructs (Gibson), and renders a
    readable DAG with minimal crossings. Can also export Mermaid.

    Notes
    -----
    * Node ``kind`` controls color in the plot.
    * Use ``add_tu_with_optional_gg`` to render TU provenance:
      parts → TU (Golden Gate result) → PCR → Gibson.
    """

    def __init__(self, title: str = "Assembly History") -> None:
        self.title = title
        self._counter = 0
        self.nodes: dict[str, HistoryNode] = {}
        self.edges: list[HistoryEdge] = []

    # --------------------------- Node/ID helpers ------------------------

    def _new_id(self) -> str:
        self._counter += 1
        return f"N{self._counter}"

    def add_source(self, label: str, record: SeqRecord, **meta: Any) -> str:
        """Add a source node (e.g., part catalog item, vector, adapter)."""
        nid = self._new_id()
        self.nodes[nid] = HistoryNode(
            nid, label, "Source", len(record), dict(meta), record
        )
        return nid

    def add_pcr(
        self,
        label: str,
        *,
        template_node: Optional[str] = None,
        result: Optional[PCRResultLike] = None,
        product: Optional[SeqRecord] = None,
        fwd_primer: Optional[str] = None,
        rev_primer: Optional[str] = None,
        **meta: Any,
    ) -> str:
        """Register a PCR product.

        Parameters
        ----------
        label
            Node label to display.
        template_node
            Node id of the template source.
        result
            PCR result object exposing a ``product: SeqRecord`` attribute.
        product
            Explicit PCR product (alternative to ``result``).
        fwd_primer, rev_primer
            Optional primer sequences; added into node meta.

        Returns
        -------
        str
            Node id of the PCR node.
        """
        if result is not None:
            rec = result.product
            # Attach common fields if present (gracefully).
            meta.setdefault("fwd_start", getattr(result, "fwd_start", None))
            meta.setdefault("rev_end", getattr(result, "rev_end", None))
            meta.setdefault("fwd_anneal_len", getattr(result, "fwd_anneal_len", None))
            meta.setdefault("rev_anneal_len", getattr(result, "rev_anneal_len", None))
        elif product is not None:
            rec = product
        else:
            raise ValueError(
                "Provide either `result` (PCRResultLike) or `product` (SeqRecord)."
            )

        if fwd_primer:
            meta.setdefault("fwd_primer", fwd_primer)
        if rev_primer:
            meta.setdefault("rev_primer", rev_primer)

        nid = self._new_id()
        self.nodes[nid] = HistoryNode(nid, label, "PCR", len(rec), dict(meta), rec)
        if template_node:
            self.edges.append(HistoryEdge(template_node, nid, "PCR"))
        return nid

    def add_golden_gate(
        self,
        label: str,
        product: SeqRecord,
        inputs: Iterable[str],
        enzyme: str = "BsaI",
        **meta: Any,
    ) -> str:
        """Add a Golden Gate node and connect input sources.

        Parameters
        ----------
        label
            Node label to display.
        product
            TU result record after Golden Gate.
        inputs
            Iterable of node ids to connect into the Golden Gate node.
        enzyme
            Enzyme label shown in edge captions.

        Returns
        -------
        str
            Node id of the Golden Gate node.
        """
        nid = self._new_id()
        imeta = dict(meta)
        imeta.setdefault("enzyme", enzyme)
        self.nodes[nid] = HistoryNode(
            nid, label, "GoldenGate", len(product), imeta, product
        )
        for src in inputs:
            self.edges.append(HistoryEdge(src, nid, f"Golden Gate ({enzyme})"))
        return nid

    def add_gibson(
        self,
        label: str,
        product: SeqRecord,
        inputs: Iterable[str],
        min_overlap: int = 40,
        circularize: bool = True,
        **meta: Any,
    ) -> str:
        """Add a Gibson node and connect input fragments.

        Parameters
        ----------
        label
            Node label to display.
        product
            Final Gibson result record.
        inputs
            Iterable of node ids to connect into the Gibson node.
        min_overlap
            Minimum overlap length (bp) displayed in edge caption.
        circularize
            Whether the final is circularized; stored in meta.

        Returns
        -------
        str
            Node id of the Gibson node.
        """
        imeta = dict(meta)
        imeta.setdefault("min_overlap", min_overlap)
        imeta.setdefault("circularize", circularize)

        nid = self._new_id()
        self.nodes[nid] = HistoryNode(
            nid, label, "Gibson", len(product), imeta, product
        )
        for src in inputs:
            self.edges.append(HistoryEdge(src, nid, f"Gibson (≥{min_overlap} bp)"))
        return nid

    def add_other(
        self, label: str, product: SeqRecord, inputs: Iterable[str], **meta: Any
    ) -> str:
        """Add a generic processing node."""
        nid = self._new_id()
        edge_label = str(meta.get("edge_label", "step"))
        self.nodes[nid] = HistoryNode(
            nid, label, "Other", len(product), dict(meta), product
        )
        for src in inputs:
            self.edges.append(HistoryEdge(src, nid, edge_label))
        return nid

    # --------------------------- TU helpers -----------------------------

    @staticmethod
    def _canonical(text: str) -> str:
        """Lowercase alphanumeric canonical form."""
        return re.sub(r"[^a-z0-9]", "", text.lower())

    @staticmethod
    def _adapter_code_from_name(name: str) -> Optional[str]:
        """Extract a two-letter adapter code (AB/BC/CD/DE) from a part name."""
        match = _ADAPTER_2L.search(str(name))
        if not match:
            return None
        code = match.group(1)
        return code if len(code) == 2 and code.isalpha() else None

    @staticmethod
    def _uns_tags_from_text(text: str) -> list[str]:
        """Return UNS tags found in text (e.g., 'UNS1', 'UNS10')."""
        return [f"UNS{n}" for n in _UNS_NUMBER.findall(str(text))]

    @staticmethod
    def _uns_tags_from_record_name(rec: SeqRecord) -> list[str]:
        """Infer UNS tags from record id/name/description."""
        pool = [
            getattr(rec, "id", ""),
            getattr(rec, "name", ""),
            getattr(rec, "description", ""),
        ]
        seen: set[str] = set()
        out: list[str] = []
        for txt in pool:
            for tag in AssemblyHistory._uns_tags_from_text(txt):
                if tag not in seen:
                    seen.add(tag)
                    out.append(tag)
        return out

    @staticmethod
    def _uns_adapters_from_features(
        rec: SeqRecord, *, default_len: int = 40
    ) -> list[tuple[str, int]]:
        """Collect UNS tags from features with lengths where unambiguous.

        If a feature references exactly one UNS tag → use its span length.
        If a feature lists multiple UNS tags → fall back to ``default_len``.
        Order is preserved; duplicates removed.
        """
        order: list[str] = []
        lengths: dict[str, int] = {}
        for feat in getattr(rec, "features", []):
            texts: list[str] = []
            for key in (
                "label",
                "note",
                "gene",
                "product",
                "locus_tag",
                "name",
                "ApEinfo_label",
            ):
                val = feat.qualifiers.get(key, [])
                if isinstance(val, list):
                    texts.extend(map(str, val))
                elif val:
                    texts.append(str(val))
            texts.append(getattr(feat, "type", ""))

            found: list[str] = []
            for t in texts:
                found += AssemblyHistory._uns_tags_from_text(t)
            if not found:
                continue

            for tag in found:
                if tag not in order:
                    order.append(tag)

            if len(found) == 1:
                tag0 = found[0]
                span = int(feat.location.end) - int(feat.location.start)
                lengths.setdefault(tag0, max(0, span))

        out: list[tuple[str, int]] = []
        seen: set[str] = set()
        for tag in order:
            if tag in seen:
                continue
            seen.add(tag)
            out.append((tag, lengths.get(tag, default_len)))
        return out

    @staticmethod
    def _uns_from_record(
        rec: SeqRecord,
        *,
        tus_row: Optional[Mapping[str, Any]] = None,
        default_len: int = 40,
    ) -> list[tuple[str, int]]:
        """Collect UNS tags from features + id/name/description + optional TUs row."""
        order: list[str] = []
        lengths: dict[str, int] = {}

        # 1) Features
        for feat in getattr(rec, "features", []):
            texts: list[str] = []
            for key in (
                "label",
                "note",
                "gene",
                "product",
                "locus_tag",
                "name",
                "ApEinfo_label",
            ):
                val = feat.qualifiers.get(key, [])
                if isinstance(val, list):
                    texts.extend(map(str, val))
                elif val:
                    texts.append(str(val))
            texts.append(getattr(feat, "type", ""))

            found: list[str] = []
            for t in texts:
                found += AssemblyHistory._uns_tags_from_text(t)
            if not found:
                continue

            for tag in found:
                if tag not in order:
                    order.append(tag)
            if len(found) == 1:
                tag0 = found[0]
                span = int(feat.location.end) - int(feat.location.start)
                lengths.setdefault(tag0, max(0, span))

        # 2) id/name/description
        for txt in (
            getattr(rec, "id", ""),
            getattr(rec, "name", ""),
            getattr(rec, "description", ""),
        ):
            for tag in AssemblyHistory._uns_tags_from_text(txt):
                if tag not in order:
                    order.append(tag)

        # 3) TUs row (UNS_Context)
        if tus_row and str(tus_row.get("UNS_Context", "")).strip():
            for tag in AssemblyHistory._uns_tags_from_text(str(tus_row["UNS_Context"])):
                if tag not in order:
                    order.append(tag)

        out: list[tuple[str, int]] = []
        seen: set[str] = set()
        for tag in order:
            if tag in seen:
                continue
            seen.add(tag)
            out.append((tag, lengths.get(tag, default_len)))
        return out

    @staticmethod
    def _resolve_part_from_catalogs(
        catalogs: Optional[Mapping[str, Mapping[str, SeqRecord]]],
        category: str,
        requested: str,
    ) -> Optional[SeqRecord]:
        """Resolve a part from catalogs with exact/CI/canonical matching."""
        if not catalogs:
            return None

        pool = (
            catalogs.get(category)
            or catalogs.get(category.capitalize())
            or catalogs.get(category.upper())
            or catalogs.get(category.lower())
        )
        if not pool:
            return None

        if requested in pool:
            return pool[requested]

        lower = {k.lower(): k for k in pool}
        if requested.lower() in lower:
            return pool[lower[requested.lower()]]

        canon = {AssemblyHistory._canonical(k): k for k in pool}
        key = canon.get(AssemblyHistory._canonical(requested))
        return pool.get(key) if key else None

    @staticmethod
    def _feature_matches_token(feat: Any, token: str) -> bool:
        """Broad token match across many qualifiers and the feature type."""
        tok_cands = [
            token,
            token.replace("_(A)", ""),
            token.replace("_CD", ""),
            token.split("_")[0],  # e.g. "P45_AB_(A)" -> "P45"
        ]
        toks = [AssemblyHistory._canonical(t) for t in tok_cands if t]
        texts: list[str] = []
        for key in (
            "label",
            "note",
            "gene",
            "product",
            "locus_tag",
            "name",
            "ApEinfo_label",
        ):
            val = feat.qualifiers.get(key, [])
            if isinstance(val, list):
                texts.extend(map(str, val))
            elif val is not None:
                texts.append(str(val))
        texts.append(getattr(feat, "type", ""))
        texts_c = [AssemblyHistory._canonical(t) for t in texts]
        return any(any(t in txt or txt in t for txt in texts_c) for t in toks)

    @staticmethod
    def slice_part_from_tu(tu_rec: SeqRecord, token: str) -> Optional[SeqRecord]:
        """Slice a part from a TU by matching a token to a TU feature."""
        for feat in getattr(tu_rec, "features", []):
            if AssemblyHistory._feature_matches_token(feat, token):
                start, end = int(feat.location.start), int(feat.location.end)
                if end > start:
                    part = tu_rec[start:end]  # keeps inner features
                    part.id = part.name = str(token)
                    part.description = f"{token} (slice from TU)"
                    return part
        return None

    @staticmethod
    def _find_feature_by_token(tu_rec: SeqRecord, token: str) -> Any | None:
        """Return first feature that matches the token or ``None``."""
        for feat in getattr(tu_rec, "features", []):
            if AssemblyHistory._feature_matches_token(feat, token):
                return feat
        return None

    @staticmethod
    def _slice_between(
        tu_rec: SeqRecord, left_feat: Any, right_feat: Any, name: str
    ) -> Optional[SeqRecord]:
        """Slice region between two TU features (left.end → right.start)."""
        if left_feat is None or right_feat is None:
            return None
        start = int(left_feat.location.end)
        end = int(right_feat.location.start)
        if end > start:
            part = tu_rec[start:end]
            part.id = part.name = name
            part.description = f"{name} (inferred by landmarks)"
            return part
        return None

    def add_tu_with_optional_gg(
        self,
        *,
        cid: str,
        tu_role: str,
        tu_rec: SeqRecord,
        tus_df: Optional[pd.DataFrame] = None,
        catalogs: Optional[Mapping[str, Mapping[str, SeqRecord]]] = None,
        enzyme: str = "BsaI",
        include_adapters: bool = True,
        adapter_style: str = "edge",  # "edge" | "node"
        include_uns: bool = True,
        uns_style: str = "node",  # "node" | "edge"
        default_uns_len: int = 40,
    ) -> str:
        """Add a TU node with upstream Golden Gate provenance.

        Parameters
        ----------
        cid
            Construct ID used to select the TUs row (if provided).
        tu_role
            TU role for this record, e.g., ``"TU1"``.
        tu_rec
            The TU sequence record (post Golden Gate) that serves as the PCR template.
        tus_df
            Optional TUs sheet with columns:
            ``ConstructID``, ``TU``, ``Promoter``, ``RBS``, ``Gene``, ``Terminator``,
            and optionally ``UNS_Context``. If present, part names are used to create
            source nodes and connect them via a Golden Gate step to the TU node.
        catalogs
            Optional nested mapping ``{category -> {name -> SeqRecord}}`` to resolve
            part sequences before falling back to slicing from the TU.
        enzyme
            Type IIS enzyme label for Golden Gate edges (purely cosmetic).
        include_adapters
            If ``True``, show Type IIS 2-letter adapter codes (AB/BC/…)
            extracted from part names.
        adapter_style
            Either ``"edge"`` (put code on the edge) or ``"node"`` (insert tiny
            adapter nodes).
        include_uns
            If ``True``, show UNS adapters carried by the TU (useful for Gibson).
        uns_style
            Either ``"node"`` (mini nodes for UNS arms) or ``"edge"`` (edge labels).
        default_uns_len
            Fallback UNS length (bp) when not annotated in features.

        Returns
        -------
        str
            Node id of the TU node. If no TUs row is available, the TU is added as a
            plain source node (no Golden Gate inputs shown).
        """
        # --- Collect Golden Gate inputs from the TUs row (if available) -----------
        part_nodes: list[tuple[str, Optional[str]]] = []  # (node_id, adapter_code)
        row_dict: dict[str, Any] | None = None

        if tus_df is not None and not tus_df.empty:
            row = tus_df[
                (tus_df["ConstructID"].astype(str) == str(cid))
                & (tus_df["TU"].astype(str).str.upper() == str(tu_role).upper())
            ]
            if not row.empty:
                row_dict = row.iloc[0].to_dict()
                for category in ("Promoter", "RBS", "Gene", "Terminator"):
                    name_obj = (row_dict or {}).get(category)
                    if name_obj is None or pd.isna(name_obj):
                        continue
                    name = str(name_obj)

                    # Resolve part: catalog → slice from TU → placeholder
                    seq = (
                        self._resolve_part_from_catalogs(catalogs, category, name)
                        or self.slice_part_from_tu(tu_rec, name)
                        or SeqRecord(
                            Seq(""),
                            id=f"{tu_role}_{category}",
                            description="placeholder",
                        )
                    )
                    node_id = self.add_source(f"{tu_role} {category}: {name}", seq)

                    code = (
                        self._adapter_code_from_name(name) if include_adapters else None
                    )
                    part_nodes.append((node_id, code))

        # --- If we have inputs, build a Golden Gate TU node -----------------------
        if part_nodes:
            gg_id = self._new_id()
            self.nodes[gg_id] = HistoryNode(
                gg_id,
                f"{tu_role} (GG result)",
                "GoldenGate",
                len(tu_rec),
                {"enzyme": enzyme},
                tu_rec,
            )

            # Connect parts → TU (with optional adapter display)
            for src_id, code in part_nodes:
                if include_adapters and code:
                    if adapter_style == "node":
                        a_id = self._new_id()
                        self.nodes[a_id] = HistoryNode(
                            a_id,
                            f"Adapter {code} (4 bp)",
                            "Source",
                            4,
                            {"adapter": code},
                            None,
                        )
                        self.edges.append(
                            HistoryEdge(src_id, a_id, f"Golden Gate ({enzyme})")
                        )
                        self.edges.append(HistoryEdge(a_id, gg_id, ""))
                    else:
                        self.edges.append(
                            HistoryEdge(
                                src_id, gg_id, f"Golden Gate ({enzyme}, {code})"
                            )
                        )
                else:
                    self.edges.append(
                        HistoryEdge(src_id, gg_id, f"Golden Gate ({enzyme})")
                    )

            # Attach UNS adapters present on the TU (for downstream Gibson)
            if include_uns:
                uns_pairs = self._uns_from_record(
                    tu_rec, tus_row=row_dict, default_len=default_uns_len
                )
                if uns_pairs:
                    if uns_style == "node":
                        for tag, length in uns_pairs:
                            u_id = self._new_id()
                            self.nodes[u_id] = HistoryNode(
                                u_id,
                                f"{tag} adapter ({int(length)} bp)",
                                "Source",
                                int(length),
                                {"uns": tag},
                                None,
                            )
                            # UNS arms conceptually feed into the TU product
                            self.edges.append(
                                HistoryEdge(u_id, gg_id, "Gibson adapter")
                            )
                    else:
                        # Edge-label style (less explicit)
                        for tag, length in uns_pairs:
                            self.edges.append(
                                HistoryEdge(gg_id, gg_id, f"{tag} ({int(length)} bp)")
                            )

            return gg_id

        # --- No TUs row/inputs → treat TU as a plain template source --------------
        return self.add_source(f"{tu_role} (template)", tu_rec)

    # --------------------------- Plot / Export --------------------------

    def plot(
        self,
        figsize: tuple[float, float] = (12.0, 7.0),
        seed: int = 42,
        layout: str = "grouped",  # "grouped" | "grid" | "dot" | "spring"
        x_gap: float = 3.0,
        y_gap: float = 1.6,
        curve: float = 0.12,
        save_path: Optional[Path] = None,
        show: bool = True,
    ) -> Figure:
        """Render the assembly DAG and optionally save the figure.

        Parameters
        ----------
        figsize : tuple of float, optional
            Matplotlib figure size (width, height) in inches. Default is ``(12.0, 7.0)``.
        seed : int, optional
            Random seed used by the spring layout. Ignored by other layouts.
        layout : {"grouped", "grid", "dot", "spring"}, optional
            Layout strategy for node positions:

            - ``"grouped"``: layered columns (Source → GoldenGate → PCR → Gibson)
            with biologically sensible ordering (default).
            - ``"grid"``: simple layered grid.
            - ``"dot"``: Graphviz DOT layout (requires ``pygraphviz``); falls back
            to ``"grid"`` if DOT is unavailable.
            - ``"spring"``: force-directed layout.

        x_gap : float, optional
            Horizontal spacing between columns (layered layouts). Default ``3.0``.
        y_gap : float, optional
            Vertical spacing between nodes within a column. Default ``1.6``.
        curve : float, optional
            Edge curvature (rad) to reduce label collisions. Default ``0.12``.
        save_path : pathlib.Path or None, optional
            If given, the figure is saved **before** display to avoid blank/white PNGs
            on some backends. Saved with ``dpi=160`` and ``bbox_inches='tight'``.
        show : bool, optional
            If ``True`` (default), call ``plt.show()`` to display the figure. Set to
            ``False`` when batch-exporting images without popping a window.

        Returns
        -------
        matplotlib.figure.Figure
            The created figure. Useful for further customization or testing.

        Notes
        -----
        - Colors are keyed by node kind: Source (grey), GoldenGate (orange),
        PCR (blue), Gibson (green), Other (light grey).
        - When ``layout="grouped"``, PCR nodes that are not fed by a TU (“orphan PCRs”
        like a backbone PCR) are aligned to the median Y of TU-PCR rows for a tidy
        look.
        - Saving is performed **before** ``plt.show()`` to prevent empty PNGs.
        """

        # Build the graph
        G = nx.DiGraph()  # noqa: N806
        for nid, node in self.nodes.items():
            G.add_node(nid, kind=node.kind, label=node.label, length=node.length)
        for e in self.edges:
            G.add_edge(e.src, e.dst, label=e.label)

        # ---- helpers for layout -------------------------------------------------
        def _role_rank(label: str) -> tuple[int, int]:
            lab = label.lower()
            m_uns = re.search(r"uns\s*([0-9]+)", lab, flags=re.IGNORECASE)
            uns_num = int(m_uns.group(1)) if m_uns else 9999
            if "promoter" in lab:
                return (0, 0)
            if "rbs" in lab:
                return (1, 0)
            if "gene" in lab and "adapter" not in lab:
                return (2, 0)
            if "terminator" in lab:
                return (3, 0)
            if "adapter" in lab and "uns" in lab:
                return (4, uns_num)
            if "adapter" in lab:
                return (5, 0)
            return (6, 0)

        def _tu_index(label: str) -> int:
            m = re.search(r"\btu\s*([0-9]+)", label, flags=re.IGNORECASE)
            return int(m.group(1)) if m else 999

        def _grid_positions() -> dict[str, tuple[float, float]]:
            order = {"Source": 0, "GoldenGate": 1, "PCR": 2, "Gibson": 3, "Other": 1}
            cols: dict[int, list[str]] = {}
            for n, d in G.nodes(data=True):
                cols.setdefault(order.get(d.get("kind", "Other"), 1), []).append(n)
            pos: dict[str, tuple[float, float]] = {}
            for ix, col in enumerate(sorted(cols)):
                nodes = cols[col]
                nodes.sort(key=lambda k: G.nodes[k].get("label", ""))
                n = len(nodes)
                for iy, node in enumerate(nodes):
                    y = -(iy - (n - 1) / 2) * y_gap
                    x = ix * x_gap
                    pos[node] = (x, y)
            return pos

        def _grouped_positions() -> dict[str, tuple[float, float]]:
            X_SRC, X_GG, X_PCR, X_GIB = (  # noqa: N806
                0.0 * x_gap,
                1.0 * x_gap,
                2.0 * x_gap,
                3.0 * x_gap,
            )

            gg_nodes = [n for n in G.nodes if G.nodes[n]["kind"] == "GoldenGate"]
            gg_nodes.sort(
                key=lambda n: (_tu_index(G.nodes[n]["label"]), G.nodes[n]["label"])
            )

            in_to: dict[str, list[str]] = {gg: [] for gg in gg_nodes}
            for u, v in G.edges():
                if v in in_to and G.nodes[u]["kind"] == "Source":
                    in_to[v].append(u)

            gg_to_pcr: dict[str, str] = {}
            for u, v in G.edges():
                if G.nodes[u]["kind"] == "GoldenGate" and G.nodes[v]["kind"] == "PCR":
                    gg_to_pcr[u] = v

            pos: dict[str, tuple[float, float]] = {}
            y_cursor = 0.0

            for gg in gg_nodes:
                sources = in_to.get(gg, [])
                sources.sort(key=lambda n: _role_rank(G.nodes[n]["label"]))
                n_src = max(1, len(sources))
                start_y = y_cursor - (n_src - 1) * (y_gap / 2.0)

                for i, s in enumerate(sources):
                    pos[s] = (X_SRC, start_y + i * y_gap)

                gg_y = start_y + (n_src - 1) * (y_gap / 2.0)
                pos[gg] = (X_GG, gg_y)

                if gg in gg_to_pcr:
                    pos[gg_to_pcr[gg]] = (X_PCR, gg_y)

                y_cursor -= (n_src + 2) * y_gap

            # Orphan PCRs (e.g., Backbone PCR)
            pcr_nodes = [n for n in G.nodes if G.nodes[n]["kind"] == "PCR"]
            placed_pcrs = set(gg_to_pcr.values())
            orphan_pcrs = [p for p in pcr_nodes if p not in placed_pcrs]
            if orphan_pcrs:
                gg_pcr_ys = [pos[p][1] for p in placed_pcrs if p in pos]
                target_y = (
                    statistics.median(gg_pcr_ys)
                    if gg_pcr_ys
                    else (y_cursor - 1.0 * y_gap)
                )
                for p in orphan_pcrs:
                    pos[p] = (X_PCR, target_y)
                    for src, _ in G.in_edges(p):
                        if G.nodes[src]["kind"] == "Source":
                            pos[src] = (X_SRC, target_y)

            # Final Gibson near centroid of its inputs
            gibsons = [n for n in G.nodes if G.nodes[n]["kind"] == "Gibson"]
            if gibsons:
                final = gibsons[0]
                preds = [u for u, v in G.edges() if v == final]
                ys = [pos[p][1] for p in preds if p in pos]
                pos[final] = (X_GIB, sum(ys) / len(ys) if ys else 0.0)

            # Fallback: place any unpositioned node near its first successor
            for n in G.nodes:
                if n in pos:
                    continue
                kind = G.nodes[n]["kind"]
                x = {
                    "Source": X_SRC,
                    "GoldenGate": X_GG,
                    "PCR": X_PCR,
                    "Gibson": X_GIB,
                }.get(kind, X_GG)
                succ = [v for _, v in G.out_edges(n)]
                y = pos[succ[0]][1] if succ and succ[0] in pos else 0.0
                pos[n] = (x, y)

            return pos

        # Choose layout
        if layout == "dot":
            try:
                pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
            except Exception:
                pos = _grid_positions()
        elif layout == "spring":
            pos = nx.spring_layout(G, seed=seed, k=1.0)
        elif layout == "grouped":
            pos = _grouped_positions()
        else:
            pos = _grid_positions()

        # Draw
        kind_colors = {
            "Source": "#bdbdbd",
            "GoldenGate": "#fdae6b",
            "PCR": "#9ecae1",
            "Gibson": "#a1d99b",
            "Other": "#d9d9d9",
        }
        node_colors = [kind_colors.get(G.nodes[n]["kind"], "#cccccc") for n in G.nodes]
        labels = {
            n: f"{G.nodes[n]['label']}\n{G.nodes[n]['length']} bp" for n in G.nodes
        }
        edge_labels = {(u, v): d.get("label", "") for u, v, d in G.edges(data=True)}

        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=figsize)
        nx.draw(
            G,
            pos,
            with_labels=False,
            node_color=node_colors,
            node_size=1600,
            arrows=True,
            connectionstyle=f"arc3,rad={curve}",
        )
        nx.draw_networkx_labels(
            G,
            pos,
            labels=labels,
            font_size=8,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.8),
        )
        nx.draw_networkx_edge_labels(
            G,
            pos,
            edge_labels=edge_labels,
            font_size=7,
            rotate=False,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7),
        )
        plt.title(self.title)
        plt.axis("off")
        plt.tight_layout()

        # Save BEFORE show to avoid blank images on some backends
        if save_path is not None:
            fig.savefig(save_path, dpi=160, bbox_inches="tight", facecolor="white")

        if show:
            plt.show()

        return fig

    def to_mermaid(self) -> str:
        """Export the DAG as a Mermaid flowchart."""

        def esc(s: str) -> str:
            return s.replace('"', r"\"")

        lines = ["flowchart LR"]
        for nid, n in self.nodes.items():
            lines.append(f'{nid}["{esc(n.label)}\\n{n.length} bp\\n{n.kind}"]')
        for e in self.edges:
            edge_label = f"|{esc(e.label)}|" if e.label else ""
            lines.append(f"{e.src} -->{edge_label} {e.dst}")
        return "\n".join(lines)


def build_histories_for_all_constructs(
    *,
    constructs_df: pd.DataFrame,
    assembly_df: pd.DataFrame,
    tus_df: Optional[pd.DataFrame],
    vector_rec: SeqRecord,
    products: dict[tuple[str, str], SeqRecord],
    pcr_results: dict[tuple[str, str], Any],
    finals: dict[str, SeqRecord],
    plot: bool = True,
    save_png_dir: Optional[Path] = None,
) -> dict[str, AssemblyHistory]:
    """Build one :class:`AssemblyHistory` per ConstructID in the sheets.

    Parameters
    ----------
    constructs_df : pandas.DataFrame
        ``Constructs`` sheet. If present, columns ``MinOverlap`` and
        ``Circularize`` are used for the Gibson node meta.
    assembly_df : pandas.DataFrame
        ``Assembly`` sheet with columns ``ConstructID``, ``FragmentRole``,
        and ``Order`` (assembly order).
    tus_df : pandas.DataFrame or None
        Optional ``TUs`` sheet (Promoter/RBS/Gene/Terminator/UNS_Context)
        used to enrich TU provenance (Golden Gate inputs / UNS adapters).
    vector_rec : Bio.SeqRecord.SeqRecord
        Backbone template record used for the backbone PCR.
    products : dict[tuple[str, str], SeqRecord]
        PCR products keyed by ``(cid, role)`` where ``role`` is
        ``"Backbone"``, ``"TU1"``, ``"TU2"``, …
        Each product may carry ``record.annotations['pcr_meta']`` with
        ``fwd_primer`` / ``rev_primer``.
    pcr_results : dict[tuple[str, str], Any]
        Optional raw PCR result objects keyed by ``(cid, role)``; if an object
        exposes a ``.product`` attribute it is passed into ``add_pcr`` so
        anneal metadata (if any) can be attached.
    finals : dict[str, SeqRecord]
        Final Gibson assemblies keyed by ``cid``.
    plot : bool, default True
        If ``True``, display a figure per construct via ``AssemblyHistory.plot()``.
    save_png_dir : pathlib.Path or None, default None
        If set, save a PNG per construct into this folder. Images are saved
        **before** showing to avoid blank/white PNGs on some backends.

    Returns
    -------
    dict[str, AssemblyHistory]
        Mapping ``ConstructID -> AssemblyHistory``.
    """
    if save_png_dir is not None:
        save_png_dir.mkdir(parents=True, exist_ok=True)

    def _norm_role(s: str) -> str:
        return str(s).strip().upper()

    def _roles_for_construct(df: pd.DataFrame, cid: str) -> list[str]:
        rows = df[df["ConstructID"].astype(str) == str(cid)]
        return rows.sort_values("Order")["FragmentRole"].astype(str).tolist()

    def _is_tu(role: str) -> bool:
        r = _norm_role(role)
        return r.startswith("TU") and r[2:].isdigit()

    def _pm(rec: SeqRecord) -> MutableMapping[str, Any]:
        return rec.annotations.get("pcr_meta", {}) or {}

    histories: dict[str, AssemblyHistory] = {}

    for cid in sorted(assembly_df["ConstructID"].astype(str).unique()):
        roles = _roles_for_construct(assembly_df, cid)
        if not roles:
            continue

        # Per-construct Gibson settings (optional)
        rowc = constructs_df[constructs_df["ConstructID"].astype(str) == str(cid)]
        min_ov = (
            int(rowc["MinOverlap"].iloc[0])
            if (not rowc.empty and "MinOverlap" in rowc)
            else 40
        )
        circ = (
            bool(rowc["Circularize"].iloc[0])
            if (not rowc.empty and "Circularize" in rowc)
            else True
        )

        hist = AssemblyHistory(title=f"{cid} assembly")
        node_by_role: dict[str, str] = {}
        pcr_node_by_role: dict[str, str] = {}

        # Backbone source
        node_by_role["BACKBONE"] = hist.add_source("Backbone (vector)", vector_rec)

        # TU provenance (Golden Gate)
        for role in roles:
            if not _is_tu(role):
                continue
            rkey = _norm_role(role)
            # Use the TU PCR product as a proxy for the TU record (fine for provenance)
            prod = products.get((cid, role)) or products.get((cid, rkey))
            if prod is None:
                raise KeyError(f"Missing TU PCR product for ({cid}, {role}).")

            node_by_role[rkey] = hist.add_tu_with_optional_gg(
                cid=cid,
                tu_role=role,
                tu_rec=prod,
                tus_df=tus_df,
                catalogs=None,
                enzyme="BsaI",
                include_adapters=True,
                adapter_style="edge",
                include_uns=True,
                uns_style="node",
                default_uns_len=40,
            )

        # PCR nodes (one per role, incl. Backbone)
        for role in roles:
            rkey = _norm_role(role)
            prod = products.get((cid, role)) or products.get((cid, rkey))
            if prod is None:
                raise KeyError(f"Missing PCR product for ({cid}, {role}).")

            meta = _pm(prod)
            res = pcr_results.get((cid, role)) or pcr_results.get((cid, rkey))

            template_node: Optional[str]
            if rkey == "BACKBONE":
                template_node = node_by_role["BACKBONE"]
                label = "Backbone PCR (with overhangs)"
            else:
                template_node = node_by_role.get(rkey)
                label = f"{role} PCR"

            use_result = hasattr(res, "product")
            n_pcr = hist.add_pcr(
                label,
                template_node=template_node,
                result=res if use_result else None,
                product=prod,
                fwd_primer=meta.get("fwd_primer"),
                rev_primer=meta.get("rev_primer"),
            )
            pcr_node_by_role[rkey] = n_pcr

        # Gibson: inputs in assembly order
        gibson_inputs: list[str] = []
        for role in roles:
            rkey = _norm_role(role)
            nid = pcr_node_by_role.get(rkey)
            if nid:
                gibson_inputs.append(nid)

        final_rec = finals.get(cid)
        if final_rec is None:
            raise KeyError(f"No final Gibson record for {cid}.")

        hist.add_gibson(
            f"{cid} final construct",
            final_rec,
            inputs=gibson_inputs,
            min_overlap=min_ov,
            circularize=circ,
        )

        # Render &/or save (save *before* show)
        save_path = (
            (save_png_dir / f"{cid}_history.png") if save_png_dir is not None else None
        )
        if plot or save_path is not None:
            # show equals `plot`; if just saving in batch, set show=False
            hist.plot(
                layout="grouped",
                figsize=(12.0, 7.0),
                x_gap=3.2,
                y_gap=1.7,
                curve=0.12,
                save_path=save_path,
                show=plot,
            )

        histories[cid] = hist

    return histories
