from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import FancyArrowPatch, Rectangle

LOGGER = logging.getLogger(__name__)


class LiquidHandlingWorkflow:
    """Container for liquid-handling workflow steps and simple visualization."""

    def __init__(self, default_plan: Optional[dict[str, Any]] = None) -> None:
        """Initialize an empty workflow.

        Parameters
        ----------
        default_plan : dict, optional
            Optional default plan structure (not used internally, kept for
            compatibility).
        """
        self.steps: list[dict[str, Any]] = []
        self.default_plan: dict[str, Any] = default_plan or {}
        self.fake_labware: dict[str, Any] = {}

    def add_step(
        self,
        process_label: str,
        source_plate: str = "hier bitte Info",
        destination_plate: str = "hier bitte Info",
        liquid_handling: str = "hier bitte Info",
        sample_handling: str = "transfer",
        fake_source_plate: Optional[str] = None,
        start: str = "No",
        additional_liquid_factor: Optional[float] = None,
    ) -> None:
        """Append a step to the workflow.

        Parameters
        ----------
        process_label : str
            Human-readable step label (e.g., "Step 1 dilute template").
        source_plate : str, optional
            Source plate name. Default ``"hier bitte Info"``.
        destination_plate : str, optional
            Destination plate name. Default ``"hier bitte Info"``.
        liquid_handling : str, optional
            Operation (e.g., ``"sample_distribution"``, ``"transfer"``).
            Default ``"hier bitte Info"``.
        sample_handling : str, optional
            Type of handling (e.g., ``"sample_distribution"``, ``"transfer"``).
            Default ``"transfer"``.
        fake_source_plate : str or None, optional
            Optional fake source plate for simulation. Default ``None``.
        start : str, optional
            Whether this is a starting step (``"Yes"``/``"No"``). Default ``"No"``.
        additional_liquid_factor : float or None, optional
            Optional factor for adjusted volumes. Default ``None``.
        """
        self.steps.append(
            {
                "process_label": process_label,
                "source_plate": source_plate,
                "destination_plate": destination_plate,
                "liquid_handling": liquid_handling,
                "sample_handling": sample_handling,
                "fake_source_plate": fake_source_plate,
                "start": start,
                "additional_liquid_factor": additional_liquid_factor,
            }
        )
        LOGGER.info("Added step: %s", process_label)

    def remove_step(self, process_label: str) -> None:
        """Remove a step by its ``process_label``.

        Parameters
        ----------
        process_label : str
            Label of the step to remove.

        Raises
        ------
        ValueError
            If no step with the given label exists.
        """
        for step in list(self.steps):
            if step.get("process_label") == process_label:
                self.steps.remove(step)
                LOGGER.info("Removed step: %s", process_label)
                return
        raise ValueError(f"Step with label {process_label!r} not found.")

    def list_steps(self) -> list[dict[str, Any]]:
        """Return all steps as a list of dictionaries."""
        return self.steps

    def generate_flowchart(
        self,
        save_path: Optional[str] = None,
        export_format: str = "pdf",
        fig_size: tuple[float, float] = (22.0, 13.0),
        dpi: int = 300,
        font_family: str = "DejaVu Sans",
        node_font_size: int = 20,
        edge_font_size: int = 15,
        title: str = "Liquid Handling Workflow",
        title_font_size: int = 24,
        box_width: float = 4.5,
        box_height: float = 1.8,
        node_facecolor: str = "#CFE2F3",
        node_edgecolor: str = "black",
        node_linewidth: float = 1.5,
        edge_color: str = "black",
        edge_width: float = 1.6,
        arrow_size: int = 22,
        color_by_plate_type: bool = True,
        plate_colors: Optional[dict[str, str]] = None,
        label_box_alpha: float = 0.95,
        show_axes: bool = False,
        transparent: bool = False,
        pad_inches: float = 0.4,
        show_plot: bool = True,
        max_labels_per_edge: Optional[int] = None,
        spring_k: float = 2.2,
        spring_iterations: int = 400,
        spring_seed: int = 42,
        scale: float = 12.0,
    ) -> None:
        """
        Generate a flowchart visualization of a liquid handling workflow.

        The workflow is represented as a directed graph where plates are nodes and
        liquid handling steps are edges. Multiple steps between the same plates are
        aggregated into a single edge label. Self-transfers are visualized as
        semicircular arrows above the corresponding plate.

        Parameters
        ----------
        save_path : str or None, optional
            Output file path. If no suffix is provided, the extension defined by
            `export_format` is used. If None, the plot is not saved.
        export_format : str, optional
            File format for export when `save_path` has no suffix (e.g. "pdf", "png").
        fig_size : tuple of float, optional
            Figure size in inches (width, height).
        dpi : int, optional
            Resolution of the output image.
        font_family : str, optional
            Font family used for all text.
        node_font_size : int, optional
            Font size for plate labels.
        edge_font_size : int, optional
            Font size for edge labels.
        title : str, optional
            Title of the flowchart.
        title_font_size : int, optional
            Font size of the title.
        box_width : float, optional
            Width of the plate rectangles.
        box_height : float, optional
            Height of the plate rectangles.
        node_facecolor : str, optional
            Default fill color for nodes.
        node_edgecolor : str, optional
            Border color of nodes.
        node_linewidth : float, optional
            Line width of node borders.
        edge_color : str, optional
            Color of edges.
        edge_width : float, optional
            Line width of edges.
        arrow_size : int, optional
            Size of arrow heads.
        color_by_plate_type : bool, optional
            If True, node colors are inferred from plate names.
        plate_colors : dict of str to str, optional
            Custom mapping of plate types to colors.
        label_box_alpha : float, optional
            Transparency of label background boxes.
        show_axes : bool, optional
            If False, axes are hidden.
        transparent : bool, optional
            Save figure with transparent background.
        pad_inches : float, optional
            Padding used when saving the figure.
        show_plot : bool, optional
            If True, display the plot (unless running in non-interactive backend).
        max_labels_per_edge : int or None, optional
            Maximum number of labels shown per edge. Additional labels are summarized.
        spring_k : float, optional
            Spring layout spacing parameter.
        spring_iterations : int, optional
            Number of iterations for layout optimization.
        spring_seed : int, optional
            Seed for deterministic layout.
        scale : float, optional
            Scaling factor for layout coordinates.

        Returns
        -------
        None
            The function renders and optionally saves the flowchart.

        Raises
        ------
        RuntimeError
            If flowchart generation fails due to unexpected errors.
        """

        try:

            def _resolve_path(path: str, fmt: str) -> Path:
                p = Path(path)
                if not p.suffix:
                    p = p.with_suffix(f".{fmt.lower()}")

                # Always save into the results folder unless already included
                if p.parts and p.parts[0] == "results":
                    final_path = p
                else:
                    final_path = Path("results") / p.name

                final_path.parent.mkdir(parents=True, exist_ok=True)
                return final_path

            def _labels(labels: list[str]) -> str:
                if max_labels_per_edge is None or len(labels) <= max_labels_per_edge:
                    return "\n".join(labels)
                visible = labels[:max_labels_per_edge]
                remaining = len(labels) - max_labels_per_edge
                return "\n".join(visible + [f"... (+{remaining} more)"])

            def _infer_color(name: str) -> str:
                if plate_colors:
                    for key, color in plate_colors.items():
                        if key.lower() in name.lower():
                            return color

                n = name.lower()
                if "water" in n:
                    return "#F9E79F"
                if "mastermix" in n:
                    return "#D6EAF8"
                if "destination" in n:
                    return "#D5F5E3"
                if "source" in n:
                    return "#D4E6F1"
                if "stock" in n:
                    return "#FAD7A0"
                return node_facecolor

            # ---------------------------
            # Collect data
            # ---------------------------
            plates: set[str] = set()
            transfers: dict[tuple[str, str], list[str]] = defaultdict(list)
            loops: dict[str, list[str]] = defaultdict(list)

            for step in self.steps:
                src = str(step.get("source_plate", "")).strip()
                dst = str(step.get("destination_plate", "")).strip()
                label = str(step.get("process_label", "")).strip()

                if src:
                    plates.add(src)
                if dst:
                    plates.add(dst)

                if src and dst:
                    if src == dst:
                        if label:
                            loops[src].append(label)
                    else:
                        if label:
                            transfers[(src, dst)].append(label)

            # ---------------------------
            # Figure setup
            # ---------------------------
            plt.rcParams["font.family"] = font_family
            fig, ax = plt.subplots(figsize=fig_size)

            ax.set_title(title, fontsize=title_font_size)

            # ---------------------------
            # Handle empty workflow
            # ---------------------------
            if not plates:
                if not show_axes:
                    ax.axis("off")

                fig.tight_layout()

                if save_path:
                    p = _resolve_path(save_path, export_format)
                    fig.savefig(
                        p,
                        dpi=dpi,
                        bbox_inches="tight",
                        pad_inches=pad_inches,
                        transparent=transparent,
                    )
                    LOGGER.info("Flowchart saved to %s", p)

                if show_plot and "agg" not in plt.get_backend().lower():
                    plt.show()
                else:
                    plt.close(fig)
                return

            # ---------------------------
            # Build graph layout
            # ---------------------------
            graph = nx.DiGraph()
            graph.add_nodes_from(plates)
            for src, dst in transfers:
                graph.add_edge(src, dst)

            pos = nx.spring_layout(
                graph,
                seed=spring_seed,
                k=spring_k,
                iterations=spring_iterations,
                scale=scale,
            )

            # ---------------------------
            # Draw nodes
            # ---------------------------
            for plate in plates:
                x, y = pos[plate]

                rect = Rectangle(
                    (x - box_width / 2, y - box_height / 2),
                    box_width,
                    box_height,
                    facecolor=(
                        _infer_color(plate) if color_by_plate_type else node_facecolor
                    ),
                    edgecolor=node_edgecolor,
                    linewidth=node_linewidth,
                    zorder=3,
                )
                ax.add_patch(rect)

                ax.text(
                    x,
                    y,
                    plate,
                    ha="center",
                    va="center",
                    fontsize=node_font_size,
                    fontweight="bold",
                    zorder=4,
                )

            # ---------------------------
            # Draw edges
            # ---------------------------
            for (src, dst), labels in transfers.items():
                x1, y1 = pos[src]
                x2, y2 = pos[dst]

                arrow = FancyArrowPatch(
                    (x1, y1),
                    (x2, y2),
                    arrowstyle="-|>",
                    mutation_scale=arrow_size,
                    linewidth=edge_width,
                    color=edge_color,
                    zorder=2,
                )
                ax.add_patch(arrow)

                ax.text(
                    (x1 + x2) / 2,
                    (y1 + y2) / 2,
                    _labels(labels),
                    ha="center",
                    va="center",
                    fontsize=edge_font_size,
                    bbox={
                        "facecolor": "white",
                        "edgecolor": "gray",
                        "alpha": label_box_alpha,
                        "boxstyle": "round,pad=0.3",
                    },
                    zorder=5,
                )

            # ---------------------------
            # Draw self-loops
            # ---------------------------
            for plate, labels in loops.items():
                x, y = pos[plate]

                loop = FancyArrowPatch(
                    (x - box_width * 0.25, y + box_height / 2),
                    (x + box_width * 0.25, y + box_height / 2),
                    connectionstyle="arc3,rad=-2.0",
                    arrowstyle="-|>",
                    mutation_scale=arrow_size + 5,
                    linewidth=edge_width,
                    color=edge_color,
                    zorder=4,
                )
                ax.add_patch(loop)

                ax.text(
                    x,
                    y + box_height * 1.3,
                    _labels(labels),
                    ha="center",
                    va="center",
                    fontsize=edge_font_size,
                    bbox={
                        "facecolor": "white",
                        "edgecolor": "gray",
                        "alpha": label_box_alpha,
                        "boxstyle": "round,pad=0.3",
                    },
                    zorder=5,
                )

            # ---------------------------
            # Final layout
            # ---------------------------
            xs = [x for x, _ in pos.values()]
            ys = [y for _, y in pos.values()]
            ax.set_xlim(min(xs) - 4, max(xs) + 4)
            ax.set_ylim(min(ys) - 3, max(ys) + 3)

            if not show_axes:
                ax.axis("off")

            fig.tight_layout()

            # ---------------------------
            # Save
            # ---------------------------
            if save_path:
                p = _resolve_path(save_path, export_format)
                fig.savefig(
                    p,
                    dpi=dpi,
                    bbox_inches="tight",
                    pad_inches=pad_inches,
                    transparent=transparent,
                )
                LOGGER.info("Flowchart saved to %s", p)

            # ---------------------------
            # Show / close
            # ---------------------------
            if show_plot and "agg" not in plt.get_backend().lower():
                plt.show()
            else:
                plt.close(fig)

        except Exception as err:
            raise RuntimeError(f"Flowchart generation failed: {err}") from err

    # def generate_flowchart(
    #     self,
    #     save_path: Optional[str] = None,
    #     fig_size: tuple[int, int] = (14, 10),
    #     node_size: int = 5000,
    #     font_size: int = 10,
    #     dpi: int = 100,
    # ) -> None:
    #     """Create a plate-centric flowchart from the workflow.

    #     Plates are nodes; edges connect source→destination. Multiple worklists
    #     between the same pair are combined as a multi-line edge label.

    #     Parameters
    #     ----------
    #     save_path : str or None, optional
    #         If provided, write the figure to this path. Otherwise show it.
    #     fig_size : tuple of int, optional
    #         Matplotlib figure size (width, height). Default ``(14, 10)``.
    #     node_size : int, optional
    #         Node size passed to `networkx.draw`. Default ``5000``.
    #     font_size : int, optional
    #         Font size for labels. Default ``10``.
    #     dpi : int, optional
    #         Figure DPI for saving. Default ``100``.
    #     """
    #     try:
    #         graph = nx.DiGraph()

    #         # Add nodes for plates
    #         for step in self.steps:
    #             src = step.get("source_plate", "")
    #             dst = step.get("destination_plate", "")
    #             if src and src not in graph:
    #                 graph.add_node(src, label=src)
    #             if dst and dst not in graph:
    #                 graph.add_node(dst, label=dst)

    #         # Aggregate edge labels per (src, dst)
    #         edge_labels: dict[tuple[str, str], list[str]] = {}
    #         for step in self.steps:
    #             src = step.get("source_plate", "")
    #             dst = step.get("destination_plate", "")
    #             label = str(step.get("process_label", ""))
    #             if not src or not dst:
    #                 continue
    #             edge_labels.setdefault((src, dst), []).append(label)

    #         # Create edges with combined labels
    #         for (src, dst), labels in edge_labels.items():
    #             combined = "\n".join(labels)
    #             graph.add_edge(src, dst, label=combined)

    #         # Layout and drawing
    #         pos = nx.spring_layout(graph, seed=42, k=0.5)

    #         plt.figure(figsize=fig_size)
    #         nx.draw(
    #             graph,
    #             pos,
    #             with_labels=True,
    #             labels={n: d.get("label", n) for n, d in graph.nodes(data=True)},
    #             node_size=node_size,
    #             node_shape="s",
    #             node_color="lightblue",
    #             font_size=font_size,
    #             font_weight="bold",
    #             edge_color="gray",
    #         )

    #         # Edge labels (handle self-loops separately)
    #         for (u, v), text in nx.get_edge_attributes(graph, "label").items():
    #             if u == v:
    #                 x, y = pos[u][0], pos[u][1] + 0.22
    #                 plt.text(
    #                     x,
    #                     y,
    #                     text,
    #                     fontsize=font_size,
    #                     ha="center",
    #                     va="center",
    #                     bbox={"facecolor": "white", "edgecolor": "black", "alpha": 0.8},
    #                 )
    #             else:
    #                 x = (pos[u][0] + pos[v][0]) / 2.0
    #                 y = (pos[u][1] + pos[v][1]) / 2.0 + 0.05
    #                 plt.text(
    #                     x,
    #                     y,
    #                     text,
    #                     fontsize=font_size,
    #                     ha="center",
    #                     va="center",
    #                     bbox={"facecolor": "white", "edgecolor": "black", "alpha": 0.8},
    #                 )

    #         plt.title("Workflow Flowchart", fontsize=font_size + 4)

    #         if save_path:
    #             plt.savefig(save_path, bbox_inches="tight", dpi=dpi)
    #             LOGGER.info("Flowchart saved as %s", save_path)
    #         else:
    #             plt.show()
    #     except Exception as err:  # pylint: disable=broad-except
    #         LOGGER.exception("Failed to generate flowchart: %s", err)

    # def to_mermaid(self, direction: str = "LR", code_block: bool = True) -> str:
    #     """Return a Mermaid flowchart definition for the workflow.

    #     Plates become nodes; each process step adds an edge from source to
    #     destination. Multiple steps between the same plates are combined as a
    #     multi-line edge label.

    #     Parameters
    #     ----------
    #     direction : {"LR", "RL", "TB", "BT"}, optional
    #         Layout direction (Left→Right, Right→Left, Top→Bottom, Bottom→Top).
    #         Default is ``"LR"``.
    #     code_block : bool, optional
    #         If ``True``, wrap the flowchart in a Markdown mermaid code fence.
    #         Default is ``True``.

    #     Returns
    #     -------
    #     str
    #         Mermaid definition, suitable for Markdown rendering if
    #         ``code_block=True``.

    #     Notes
    #     -----
    #     - Node IDs are sanitized to alphanumeric/underscore and forced to
    #     start with a letter to satisfy Mermaid identifier rules.
    #     - Edge labels use ``<br/>`` to display multiple step labels.
    #     """

    #     def _sanitize_id(name: str) -> str:
    #         # keep alnum and underscore; ensure starts with a letter
    #         import re as _re  # local to avoid polluting module namespace

    #         cleaned = _re.sub(r"\W+", "_", str(name))
    #         if not cleaned or not cleaned[0].isalpha():
    #             cleaned = f"N_{cleaned or 'node'}"
    #         return cleaned

    #     # Collect nodes (plates)
    #     nodes: set[str] = set()
    #     for step in self.steps:
    #         src = step.get("source_plate", "")
    #         dst = step.get("destination_plate", "")
    #         if src:
    #             nodes.add(src)
    #         if dst:
    #             nodes.add(dst)

    #     # Aggregate edge labels per (src, dst)
    #     edge_labels: dict[tuple[str, str], list[str]] = {}
    #     for step in self.steps:
    #         src = step.get("source_plate", "")
    #         dst = step.get("destination_plate", "")
    #         if not src or not dst:
    #             continue
    #         edge_labels.setdefault((src, dst), []).append(
    #             str(step.get("process_label", ""))
    #         )

    #     # Build Mermaid lines
    #     direction = direction.upper()
    #     if direction not in {"LR", "RL", "TB", "BT"}:
    #         direction = "LR"

    #     lines: list[str] = [f"flowchart {direction}"]

    #     # Declare nodes with labels
    #     for label in sorted(nodes):
    #         nid = _sanitize_id(label)
    #         # square nodes for plates
    #         lines.append(f'  {nid}["{label}"]')

    #     # Declare edges
    #     for (src, dst), labels in edge_labels.items():
    #         sid = _sanitize_id(src)
    #         did = _sanitize_id(dst)
    #         # join multiple process labels
    #         joined = "<br/>".join(label for label in labels if label)
    #         if joined:
    #             lines.append(f"  {sid} -->|{joined}| {did}")
    #         else:
    #             lines.append(f"  {sid} --> {did}")

    #     mermaid = "\n".join(lines)
    #     if code_block:
    #         return f"```mermaid\n{mermaid}\n```"
    #     return mermaid
