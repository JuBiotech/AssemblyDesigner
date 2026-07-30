from __future__ import annotations

import csv
import importlib
import io
import os
import zipfile
from pathlib import Path
from typing import Any

# Headless plotting for CI
import matplotlib
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

# =========================
# Dynamic import helpers
# =========================

REQUIRED_EXPORTS = {
    "load_dna_file",
    "filter_features",
    "snapgene_to_seqrecord",
    "load_parts_from_folders",
    "generate_assembly_reports",
    "organize_assembly_reports",
    "delete_all_zip_files",
    "run_all_functions",
    "remove_near_duplicate_features",
    "get_ladder",
    "perform_restriction_digest",
    "read_plasmid_file",
    "plot_gel",
    "digest_and_plot_plasmids",
    "_safe_filename",
}


def _try_import(name: str):
    try:
        return importlib.import_module(name)
    except Exception:
        return None


def _module_has_api(m) -> bool:
    return m is not None and all(hasattr(m, attr) for attr in REQUIRED_EXPORTS)


def _import_plasmidio():
    """
    Import the plasmidio module. Prefer env var AD_PLASMIDIO_MODULE.
    If the package lacks re-exports, fall back to its implementation submodule.
    """
    candidates = []
    env = os.environ.get("AD_PLASMIDIO_MODULE")
    if env:
        candidates.extend([env, env + ".plasmidio"])

    # Common locations
    bases = [
        "assembly_designer.plasmidio",
        "assembly_designer.manager.plasmidio",
    ]
    for base in bases:
        candidates.extend([base, base + ".plasmidio"])

    # Try each, preferring ones that expose the required API
    for name in candidates:
        m = _try_import(name)
        if _module_has_api(m):
            return m

    # As a last resort, return the first importable module and let AttributeError surface
    for name in candidates:
        m = _try_import(name)
        if m is not None:
            return m

    raise ImportError(
        "Could not import plasmidio module. Set AD_PLASMIDIO_MODULE, e.g. "
        "'assembly_designer.plasmidio', or re-export functions in __init__.py."
    )


mod = _import_plasmidio()

# Submodules for monkeypatching attributes like `dc` / `snapgene_file_to_dict` that
# now live in their own per-concern files rather than a single implementation module.
import assembly_designer.plasmidio.assembly as assembly_mod  # noqa: E402
import assembly_designer.plasmidio.io_ as io_mod  # noqa: E402
import assembly_designer.plasmidio.pcr as pcr_mod  # noqa: E402
import assembly_designer.plasmidio.s1_doc as s1_doc_mod  # noqa: E402

# Bind API under test
load_dna_file = mod.load_dna_file
filter_features = mod.filter_features
snapgene_to_seqrecord = mod.snapgene_to_seqrecord
load_parts_from_folders = mod.load_parts_from_folders
generate_assembly_reports = mod.generate_assembly_reports
organize_assembly_reports = mod.organize_assembly_reports
delete_all_zip_files = mod.delete_all_zip_files
run_all_functions = mod.run_all_functions
remove_near_duplicate_features = mod.remove_near_duplicate_features
get_ladder = mod.get_ladder
perform_restriction_digest = mod.perform_restriction_digest
read_plasmid_file = mod.read_plasmid_file
plot_gel = mod.plot_gel
digest_and_plot_plasmids = mod.digest_and_plot_plasmids
_safe_filename = mod._safe_filename  # testing private util intentionally

# Optional exports used by additional tests (guarded if absent)
longest_overlap = getattr(mod, "longest_overlap", None)
merge_with_gibson_features = getattr(mod, "merge_with_gibson_features", None)
simulate_pcr = getattr(mod, "simulate_pcr", None)
simulate_pcr_overhangs = getattr(mod, "simulate_pcr_overhangs", None)
AssemblyHistory = getattr(mod, "AssemblyHistory", None)
build_histories_for_all_constructs = getattr(
    mod, "build_histories_for_all_constructs", None
)
make_template_lookup = getattr(mod, "make_template_lookup", None)
run_3g_batch_safe = getattr(mod, "run_3g_batch_safe", None)
BuildStatus = getattr(mod, "BuildStatus", None)


# =========================
# Local helpers
# =========================


def _mk_record(seq: str, *, rid: str = "rec", rname: str = "rec") -> SeqRecord:
    rec = SeqRecord(Seq(seq), id=rid, name=rname)
    rec.annotations["molecule_type"] = "DNA"
    return rec


def _write_gb(path: Path, record: SeqRecord) -> None:
    with path.open("w") as fh:
        SeqIO.write(record, fh, "genbank")


# =========================
# Existing (your) tests
# =========================


def test_read_plasmid_file_sets_unknown_id_name(tmp_path: Path) -> None:
    gb = tmp_path / "unk.gb"
    rec = _mk_record("ATGCATGCAT", rid="Unknown", rname="Unknown")
    _write_gb(gb, rec)
    out = read_plasmid_file(str(gb))
    assert out.id == "unk" and out.name == "unk"


def test_remove_near_duplicate_features(tmp_path: Path) -> None:
    gb = tmp_path / "dup.gb"
    rec = _mk_record("A" * 50, rid="dup", rname="dup")
    rec.features = [
        SeqFeature(
            FeatureLocation(5, 20, strand=1), type="gene", qualifiers={"label": ["z"]}
        ),
        SeqFeature(
            FeatureLocation(6, 21, strand=1), type="gene", qualifiers={"label": ["z"]}
        ),
    ]
    _write_gb(gb, rec)
    cleaned = remove_near_duplicate_features(gb, tolerance=3)
    assert cleaned.exists()
    got = SeqIO.read(str(cleaned), "genbank")
    assert len(got.features) == 1


def test_load_parts_from_folders(tmp_path: Path) -> None:
    f1 = tmp_path / "catA"
    f2 = tmp_path / "catB"
    f1.mkdir()
    f2.mkdir()
    _write_gb(f1 / "A1.gb", _mk_record("ATGCATGC", rid="a1", rname="a1"))
    _write_gb(f2 / "B1.gb", _mk_record("ATGCAT", rid="b1", rname="b1"))

    loaded = load_parts_from_folders([str(f1), str(f2)])
    assert set(loaded.keys()) == {"catA", "catB"}
    assert "A1" in loaded["catA"] and "B1" in loaded["catB"]
    assert str(loaded["catA"]["A1"].seq) == "ATGCATGC"


def test_safe_filename_sanitization_and_truncation() -> None:
    unsafe = 'weird name<>:"/\\|?* and spaces'
    out = _safe_filename(unsafe)
    assert all(ch not in out for ch in '<>:"/\\|?* ')
    long = "x" * 200
    out2 = _safe_filename(long, maxlen=32)
    assert len(out2.split("-")[-1]) == 8 and len(out2) <= 32


def test_get_ladder_and_digest() -> None:
    sizes = get_ladder("2kb")
    assert isinstance(sizes, list) and all(isinstance(s, int) for s in sizes)
    with pytest.raises(ValueError):
        get_ladder("nope")

    # If no cuts, a single fragment == len(seq). RestrictionBatch.search expects a Seq.
    seq = Seq("A" * 1234)
    frags = perform_restriction_digest(seq, enzymes=[])  # empty batch
    assert frags == [len(seq)]


def test_plot_gel_saves_image(tmp_path: Path) -> None:
    out = tmp_path / "gel.png"
    fragments_by_lane: dict[str, list[int]] = {"Sample1": [500, 250, 125]}
    plot_gel(
        fragments_by_lane,
        lane_labels=["Ladder", "Sample1"],
        ladder_type="1kb",
        save_path=str(out),
    )
    assert out.exists() and out.stat().st_size > 0


def test_digest_and_plot_plasmids(tmp_path: Path) -> None:
    gb = tmp_path / "p.gb"
    _write_gb(gb, _mk_record("A" * 1000, rid="p", rname="p"))
    out = tmp_path / "gel2.png"
    digest_and_plot_plasmids([(str(gb), [])], ladder_type="1kb", save_path=str(out))
    assert out.exists() and out.stat().st_size > 0


def test_snapgene_to_seqrecord_mocked(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def _fake_reader(_: str) -> dict[str, Any]:
        return {
            "seq": "ATGCATGC",
            "features": [
                {
                    "start": 0,
                    "end": 4,
                    "strand": "+",
                    "type": "CDS",
                    "qualifiers": {"label": ["abc"]},
                },
                {
                    "start": 0,
                    "end": 4,
                    "strand": "+",
                    "type": "CDS",
                    "qualifiers": {"label": ["abc"]},
                },
            ],
            "dna": {"topology": "circular"},
            "notes": {"author": "unit-test"},
            "id": "snap1",
            "name": "snap1",
        }

    # Patch attribute on the actual implementation module
    monkeypatch.setattr(io_mod, "snapgene_file_to_dict", _fake_reader, raising=True)
    fake_path = tmp_path / "x.dna"
    fake_path.write_bytes(b"\x00\x01")
    rec = snapgene_to_seqrecord(str(fake_path))
    assert isinstance(rec, SeqRecord)
    assert rec.annotations.get("topology") == "circular"
    assert len(rec.features) == 1  # duplicates filtered
    assert str(rec.seq) == "ATGCATGC"


class _DummySim:
    def __init__(self, zip_target: Path, gb_payload: SeqRecord):
        self._zip_target = zip_target
        self._gb_payload = gb_payload

    def write_report(self, out_path: str) -> None:
        out = Path(out_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(out, "w") as zf:
            gb_bytes = io.StringIO()
            SeqIO.write(self._gb_payload, gb_bytes, "genbank")
            zf.writestr("constructed/construct.gb", gb_bytes.getvalue())


class _DummyAssembly:
    def __init__(self, parts: list[str], name: str, gb_payload: SeqRecord):
        self._parts = parts
        self._name = name
        self._gb_payload = gb_payload

    def simulate(self, sequence_repository: Any) -> _DummySim:  # noqa: ANN401
        zip_target = Path(f"{self._name}_report.zip")
        return _DummySim(zip_target, self._gb_payload)


class _DummyDC:
    class SequenceRepository:
        def __init__(self, collections: dict[str, dict[str, SeqRecord]]):
            self.collections = collections

    def __init__(self):
        self.Type2sRestrictionAssembly = self._factory

    def _factory(self, parts: list[str], name: str) -> _DummyAssembly:
        payload = _mk_record("ATGCATGCATGC", rid=name, rname=name)
        return _DummyAssembly(parts, name, payload)


def _write_minimal_gb(dirpath: Path, name: str) -> Path:
    p = dirpath / f"{name}.gb"
    rec = _mk_record("ATGCATGC", rid=name, rname=name)
    rec.features = [
        SeqFeature(FeatureLocation(0, 4), type="gene", qualifiers={"label": ["g"]})
    ]
    _write_gb(p, rec)
    return p


def test_run_all_functions_wrapper(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    cat1 = tmp_path / "C1"
    cat2 = tmp_path / "C2"
    cat1.mkdir()
    cat2.mkdir()
    _write_minimal_gb(cat1, "A")
    _write_minimal_gb(cat2, "B")

    monkeypatch.setattr(assembly_mod, "dc", _DummyDC(), raising=True)
    outputs = run_all_functions(
        [str(cat1), str(cat2)], report_dir=str(tmp_path / "R"), delete_zip=False
    )
    assert outputs and all(p.exists() for p in outputs)


def test_delete_all_zip_files(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    z1 = tmp_path / "x.zip"
    z2 = tmp_path / "y.zip"
    z1.write_bytes(b"X")
    z2.write_bytes(b"Y")
    delete_all_zip_files(str(tmp_path))
    out = capsys.readouterr().out
    assert "Deleted: x.zip" in out and "Deleted: y.zip" in out
    assert not z1.exists() and not z2.exists()


def _make_report_zip(zip_path: Path, *, files: dict[str, bytes]) -> Path:
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w") as zf:
        for name, content in files.items():
            zf.writestr(name, content)
    return zip_path


def test_organize_assembly_reports_skips_failed_assembly(tmp_path: Path) -> None:
    """A failed DNAcauldron assembly writes error.csv and no assembled construct.

    organize_assembly_reports() must not fall back to an input part (e.g. the
    bare backbone) and save it under the assembly's name, since that would make
    a failed assembly look like it succeeded. Regression test for the bug where
    Golden Gate combinations with an incompatible backbone still produced a
    "constructed" GenBank file in reports/Assembly.
    """
    report_dir = tmp_path / "reports"
    _make_report_zip(
        report_dir / "asmA_report.zip",
        files={
            "provided_parts_records/backbone.gb": b"LOCUS backbone\n//\n",
            "error.csv": (
                "assembly_name;message;suggestion;data\n"
                "asmA;Wrong number of constructs;expected_: 1,found: 0\n"
            ).encode(),
        },
    )

    extracted = organize_assembly_reports(report_dir, delete_zip=False)

    assert extracted == []
    assert not (report_dir / "Assembly" / "asmA.gb").exists()


def test_organize_assembly_reports_extracts_successful_assembly(
    tmp_path: Path,
) -> None:
    """Sanity check: a zip with a real construct is unaffected by the error.csv
    guard and still extracts the exact-match GenBank file."""
    report_dir = tmp_path / "reports"
    _make_report_zip(
        report_dir / "asmB_report.zip",
        files={
            "asmB.gb": b"LOCUS asmB construct\n//\n",
            "provided_parts_records/backbone.gb": b"LOCUS backbone\n//\n",
        },
    )

    extracted = organize_assembly_reports(report_dir, delete_zip=False)

    assert [p.name for p in extracted] == ["asmB.gb"]
    assert (
        report_dir / "Assembly" / "asmB.gb"
    ).read_bytes() == b"LOCUS asmB construct\n//\n"


# =========================
# Additional tests
# =========================


@pytest.mark.skipif(longest_overlap is None, reason="longest_overlap not exported")
def test_longest_overlap_basic() -> None:
    assert longest_overlap("ABCD", "CDEF", 1) >= 2  # 'CD'
    assert longest_overlap("XYZ", "ABC", 1) == 0


@pytest.mark.skipif(
    merge_with_gibson_features is None, reason="merge_with_gibson_features not exported"
)
def test_merge_with_gibson_features_linear_with_feature_preservation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Provide a simple feature-preserving concat implementation via monkeypatch
    def _concat(parts: list[SeqRecord], *, new_id: str) -> SeqRecord:
        assert len(parts) == 2
        a, b = parts
        out = SeqRecord(Seq(str(a.seq) + str(b.seq)), id=new_id, name=new_id)
        out.annotations["molecule_type"] = "DNA"

        def _strand_of(feat: SeqFeature):
            # robust: in case strand is missing → None
            loc = getattr(feat, "location", None)
            return getattr(loc, "strand", None)

        shift = len(a)
        feats: list[SeqFeature] = []

        # Adopt features from A 1:1
        for f in getattr(a, "features", []):
            feats.append(
                SeqFeature(
                    FeatureLocation(
                        int(f.location.start),
                        int(f.location.end),
                        strand=_strand_of(f),
                    ),
                    type=f.type,
                    qualifiers=dict(f.qualifiers),
                )
            )

        # Take features from B, but shift them by len(a)
        for f in getattr(b, "features", []):
            s = int(f.location.start) + shift
            e = int(f.location.end) + shift
            feats.append(
                SeqFeature(
                    FeatureLocation(s, e, strand=_strand_of(f)),
                    type=f.type,
                    qualifiers=dict(f.qualifiers),
                )
            )

        out.features = feats
        return out

    monkeypatch.setattr(pcr_mod, "_concat_records_with_features", _concat, raising=True)

    a = _mk_record("AAAACCC", rid="A", rname="A")
    b = _mk_record("CCCGGG", rid="B", rname="B")
    # add simple features
    a.features = [
        SeqFeature(
            FeatureLocation(0, 4), type="misc_feature", qualifiers={"label": ["partA"]}
        )
    ]
    b.features = [
        SeqFeature(
            FeatureLocation(3, 6), type="misc_feature", qualifiers={"label": ["partB"]}
        )
    ]

    out = merge_with_gibson_features([a, b], min_overlap=3, circularize=False)
    assert str(out.seq) == "AAAACCCGGG"
    # features preserved and shifted
    assert any("partA" in f.qualifiers.get("label", [""])[0] for f in out.features)
    assert any("partB" in f.qualifiers.get("label", [""])[0] for f in out.features)


@pytest.mark.skipif(AssemblyHistory is None, reason="AssemblyHistory not available")
def test_assembly_history_plot_and_export(tmp_path: Path) -> None:
    # Minimal history: source -> PCR -> Gibson
    hist = AssemblyHistory(title="demo")
    src = _mk_record("AAAA", rid="s", rname="s")
    src_node = hist.add_source("Backbone (vector)", src)
    pcr_prod = _mk_record("AAAATTTT", rid="p", rname="p")
    n_pcr = hist.add_pcr("Backbone PCR", template_node=src_node, product=pcr_prod)
    final = _mk_record("AAAATTTTCCCC", rid="f", rname="f")
    hist.add_gibson("final", final, inputs=[n_pcr], min_overlap=4, circularize=True)
    # Plot (headless) and save figure
    hist.plot(layout="grouped", figsize=(6, 4), x_gap=2.0, y_gap=1.2, curve=0.1)
    out = tmp_path / "hist.png"
    plt.gcf().savefig(out, dpi=120, bbox_inches="tight")
    assert out.exists() and out.stat().st_size > 0
    # Mermaid export should be a non-empty string
    mermaid = hist.to_mermaid()
    assert isinstance(mermaid, str) and "flowchart" in mermaid


@pytest.mark.skipif(
    simulate_pcr_overhangs is None, reason="simulate_pcr_overhangs not exported"
)
def test_simulate_pcr_overhangs_included() -> None:
    template = _mk_record("ATGCGTAAAGTTAGC", rid="tpl2", rname="tpl2")
    overhang = "AAAA"
    fwd_anneal = "ATGCGT"
    fwd = overhang + fwd_anneal
    rev = "GCTAAC"  # RC for "...GTTAGC"

    res = simulate_pcr_overhangs(
        template,
        fwd_primer=fwd,
        rev_primer=rev,
        include_overhangs=True,
        product_id="amp",
        trim_terminal_N=False,
        min_anneal=6,  # Annealing-Lengths in this test
    )

    prod = res.product if hasattr(res, "product") else res
    seq = str(prod.seq)

    # The Overhang + Anneal section are indeed intended to be located at the front of the product
    assert seq.startswith(overhang + fwd_anneal)


@pytest.mark.skipif(
    build_histories_for_all_constructs is None, reason="history builder not exported"
)
def test_build_histories_for_all_constructs_minimal(tmp_path: Path) -> None:
    # Sheets
    constructs_df = pd.DataFrame(
        [{"ConstructID": "C1", "MinOverlap": 20, "Circularize": True}]
    )
    assembly_df = pd.DataFrame(
        [
            {"ConstructID": "C1", "Order": 1, "FragmentRole": "Backbone"},
            {"ConstructID": "C1", "Order": 2, "FragmentRole": "TU1"},
        ]
    )
    tus_df = pd.DataFrame(
        [
            {
                "ConstructID": "C1",
                "TU": "TU1",
                "Promoter": "P1_AB",
                "RBS": "R1_BC",
                "Gene": "G1_CD",
                "Terminator": "T1_DE",
                "UNS_Context": "UNS1_UNS10",
            }
        ]
    )

    vector_rec = _mk_record("A" * 60, rid="vec", rname="vec")
    tu_rec = _mk_record("A" * 40 + "C" * 20, rid="tu1", rname="tu1")

    products = {
        ("C1", "Backbone"): _mk_record("A" * 60, rid="bb", rname="bb"),
        ("C1", "TU1"): tu_rec,
    }
    pcr_results = {}  # optional
    finals = {"C1": _mk_record("A" * 60 + "C" * 20, rid="final", rname="final")}

    histories = build_histories_for_all_constructs(
        constructs_df=constructs_df,
        assembly_df=assembly_df,
        tus_df=tus_df,
        vector_rec=vector_rec,
        products=products,
        pcr_results=pcr_results,
        finals=finals,
        plot=False,
        save_png_dir=None,
    )
    assert "C1" in histories
    assert hasattr(histories["C1"], "nodes") and hasattr(histories["C1"], "edges")


@pytest.mark.skipif(
    make_template_lookup is None, reason="make_template_lookup not exported"
)
def test_make_template_lookup_selects_tu_by_uns(tmp_path: Path) -> None:
    reports = tmp_path / "reports"
    asm_dir = reports / "Assembly"
    asm_dir.mkdir(parents=True, exist_ok=True)

    # TU record with UNS1 mentioned in description & features
    tu = _mk_record("ATGCATGCATGC", rid="TU1", rname="TU1")
    tu.description = "TU with UNS1 and UNS10"
    tu.features = [
        SeqFeature(
            FeatureLocation(0, 4),
            type="misc_feature",
            qualifiers={"label": ["UNS1 adapter"]},
        )
    ]
    tu_path = asm_dir / "TU1_final.gb"
    _write_gb(tu_path, tu)

    # Vector
    vector_rec = _mk_record("A" * 50, rid="vec", rname="vec")

    # Minimal TUs sheet row with UNS_Context
    tus_df = pd.DataFrame(
        [{"ConstructID": "C1", "TU": "TU1", "UNS_Context": "UNS1_UNS10"}]
    )

    lookup = make_template_lookup(
        reports_dir=str(reports),
        vector_rec=vector_rec,
        tus_df=tus_df,
    )

    # TU1 should resolve to our TU record
    tu_sel = lookup("C1", "TU1")
    assert isinstance(tu_sel, SeqRecord)
    assert (
        str(tu_sel.id) == "TU1_final"
        or str(tu_sel.name) == "TU1_final"
        or "UNS1" in tu_sel.description
    )

    # Backbone should return the given vector
    bb_sel = lookup("C1", "Backbone")
    assert str(bb_sel.seq) == str(vector_rec.seq)


# -------- S1 documentation tests -------------------------------------------------

# Bind S1 API (skip if not exported/reachable)
BUILD_S1 = getattr(mod, "build_s1_documentation", None)
RESOLVE_DONOR = getattr(s1_doc_mod, "resolve_donor", None)


@pytest.mark.skipif(BUILD_S1 is None, reason="build_s1_documentation not exported")
def test_build_s1_writes_excel_in_cwd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Excel should be written to CWD (not the GenBank folder)."""
    gb_dir = tmp_path / "gb"
    out_dir = tmp_path / "out"
    gb_dir.mkdir()
    out_dir.mkdir()

    # Make a minimal GB with clear Amp marker
    rec = _mk_record("ATGCATGCATGC", rid="X1", rname="X1")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["organism"] = "Escherichia coli"  # donor from record
    rec.features = [
        SeqFeature(
            FeatureLocation(0, len(rec)),
            type="source",
            qualifiers={"organism": ["Escherichia coli"]},
        ),
        SeqFeature(
            FeatureLocation(0, 6), type="promoter", qualifiers={"label": ["J23106"]}
        ),
        SeqFeature(
            FeatureLocation(6, 9),
            type="misc_binding",
            qualifiers={"label": ["B0032 RBS"]},
        ),
        SeqFeature(FeatureLocation(9, 12), type="CDS", qualifiers={"gene": ["panD"]}),
        SeqFeature(
            FeatureLocation(0, 3),
            type="misc_feature",
            qualifiers={"note": ["ampicillin resistance"]},
        ),
    ]
    gb_path = gb_dir / "J23106_B0032_panD.gb"
    _write_gb(gb_path, rec)

    # CWD should be out_dir
    monkeypatch.chdir(out_dir)
    out_file = BUILD_S1(str(gb_dir))
    assert out_file == out_dir / "S1_Dokumentation.xlsx"
    assert out_file.exists()

    df = pd.read_excel(out_file)
    assert not df.empty
    assert set(
        ["Nr.", "Spender", "Ausgangsvektor", "Gen/Promotor/RBS", "Bezeichnung"]
    ).issubset(df.columns)
    # donor from record
    assert df.loc[0, "Spender"] == "Escherichia coli"
    # vector detection saw Amp
    assert "DVA" in str(df.loc[0, "Ausgangsvektor"])
    # parts from features (promoter/RBS/gene)
    text = str(df.loc[0, "Gen/Promotor/RBS"])
    assert "J23106" in text and "B0032" in text and "panD" in text
    # Bezeichnung is the stem
    assert df.loc[0, "Bezeichnung"] == "J23106_B0032_panD"


@pytest.mark.skipif(RESOLVE_DONOR is None, reason="resolve_donor not available")
def test_resolve_donor_from_annotations_and_source(tmp_path: Path) -> None:
    """Donor should be found in annotations['organism'] or /source feature."""
    # Case 1: annotations.organism
    r1 = _mk_record("AAAAAA", rid="A", rname="A")
    r1.annotations["organism"] = "Bacillus subtilis"
    assert RESOLVE_DONOR(r1, file_stem="A") == "Bacillus subtilis"

    # Case 2: via source feature
    r2 = _mk_record("AAAAAA", rid="B", rname="B")
    r2.annotations.pop("organism", None)
    r2.features = [
        SeqFeature(
            FeatureLocation(0, len(r2)),
            type="source",
            qualifiers={"organism": ["Corynebacterium glutamicum"]},
        )
    ]
    assert RESOLVE_DONOR(r2, file_stem="B") == "Corynebacterium glutamicum"


@pytest.mark.skipif(BUILD_S1 is None, reason="build_s1_documentation not exported")
def test_build_s1_include_glob_filters(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    gb_dir = tmp_path / "gb"
    gb_dir.mkdir()
    # two files, only one matches glob
    _write_gb(gb_dir / "TU1_part1.gb", _mk_record("AAAA", rid="r1", rname="r1"))
    _write_gb(gb_dir / "TU2_part2.gb", _mk_record("AAAA", rid="r2", rname="r2"))
    monkeypatch.chdir(tmp_path)
    out_file = BUILD_S1(str(gb_dir), include_glob="*TU1*.gb*")
    df = pd.read_excel(out_file)
    assert len(df) == 1
    assert df.loc[0, "Bezeichnung"] == "TU1_part1"


@pytest.mark.skipif(BUILD_S1 is None, reason="build_s1_documentation not exported")
def test_build_s1_uses_curated_mapping_over_anything(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CSV mapping (by stem / gene / accession) should override."""
    gb_dir = tmp_path / "gb"
    out_dir = tmp_path / "out"
    gb_dir.mkdir()
    out_dir.mkdir()

    # Record with no donor inside
    rec = _mk_record("ATGCATGC", rid="R", rname="R")
    rec.annotations.pop("organism", None)
    rec.annotations["accessions"] = ["AB123456"]  # for accession-based mapping
    rec.features = [
        SeqFeature(FeatureLocation(0, 4), type="CDS", qualifiers={"gene": ["fooX"]}),
    ]
    stem = "my_construct_01"
    _write_gb(gb_dir / f"{stem}.gb", rec)

    # Curated mapping CSV covering three key types; stem takes precedence in our lookup order
    csv_path = tmp_path / "donor_mapping.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["stem", "donor", "gene", "accession"])
        w.writeheader()
        w.writerow({"stem": stem, "donor": "StemDonor", "gene": "", "accession": ""})
        w.writerow({"stem": "", "donor": "GeneDonor", "gene": "foox", "accession": ""})
        w.writerow(
            {"stem": "", "donor": "AccDonor", "gene": "", "accession": "AB123456"}
        )

    monkeypatch.chdir(out_dir)
    out_file = BUILD_S1(str(gb_dir), donor_mapping_csv=csv_path, prefer_features=True)
    df = pd.read_excel(out_file)
    assert df.loc[0, "Spender"] == "StemDonor"  # stem match used


@pytest.mark.skipif(RESOLVE_DONOR is None, reason="resolve_donor not available")
def test_resolve_donor_ncbi_fallback_mocked(monkeypatch: pytest.MonkeyPatch) -> None:
    """NCBI fallback should be used only if enabled, and we mock Entrez to avoid network."""
    # Record without local donor, but with an accession to drive lookup
    rec = _mk_record("ATGC", rid="Y", rname="Y")
    rec.annotations.pop("organism", None)
    rec.annotations["accessions"] = ["XYZ999"]

    # Build a dummy Entrez shim that returns a summary containing Organism
    class _DummyHandle:
        def close(self):
            pass

    class _DummyEntrez:
        email = None
        api_key = None

        @staticmethod
        def esummary(db: str, id: str, retmode: str):  # noqa: A003
            return _DummyHandle()

        @staticmethod
        def read(handle):  # noqa: A003
            return [{"Organism": "Lactococcus lactis"}]

    # Patch the Entrez symbol used *inside the implementation module*
    monkeypatch.setattr(s1_doc_mod, "Entrez", _DummyEntrez, raising=True)
    # Avoid real sleeping if your code rate-limits
    monkeypatch.setattr(
        s1_doc_mod,
        "time",
        type("T", (), {"sleep": staticmethod(lambda *_: None)}),
        raising=False,
    )

    donor_disabled = RESOLVE_DONOR(
        rec, file_stem="Y", enable_ncbi_fallback=False, ncbi_email="x@y.z"
    )
    assert donor_disabled is None  # fallback off → still unknown

    donor_enabled = RESOLVE_DONOR(
        rec,
        file_stem="Y",
        enable_ncbi_fallback=True,
        ncbi_email="x@y.z",
        ncbi_api_key=None,
        ncbi_cache=None,
    )
    assert donor_enabled == "Lactococcus lactis"


@pytest.mark.skipif(BUILD_S1 is None, reason="build_s1_documentation not exported")
def test_build_s1_prefers_features_over_filename(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """When prefer_features=True, feature labels beat filename heuristics."""
    gb_dir = tmp_path / "gb"
    gb_dir.mkdir()
    # filename suggests J23107 / B0015, but features say J23106 / B0032 / gene panD
    rec = _mk_record("ATGCATGCATGC", rid="rec1", rname="rec1")
    rec.features = [
        SeqFeature(
            FeatureLocation(0, 6), type="promoter", qualifiers={"label": ["J23106"]}
        ),
        SeqFeature(
            FeatureLocation(6, 9), type="misc_binding", qualifiers={"label": ["B0032"]}
        ),
        SeqFeature(FeatureLocation(9, 12), type="CDS", qualifiers={"gene": ["panD"]}),
    ]
    gb = gb_dir / "J23107_B0015_foo.gb"
    _write_gb(gb, rec)

    monkeypatch.chdir(tmp_path)
    out = BUILD_S1(str(gb_dir), prefer_features=True)
    df = pd.read_excel(out)
    text = str(df.loc[0, "Gen/Promotor/RBS"])
    assert "J23106" in text and "B0032" in text and "panD" in text


def _mk_feature(start: int, end: int, ftype: str, label: str) -> SeqFeature:
    """Small helper to build a GenBank feature with a visible display name."""
    return SeqFeature(
        FeatureLocation(start, end, strand=1),
        type=ftype,
        qualifiers={"label": [label], "note": [label]},
    )


def _write_gb(path: Path, rec: SeqRecord) -> None:
    """Write a SeqRecord to GenBank with minimal required annotations."""
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations.setdefault("topology", "circular")
    SeqIO.write(rec, str(path), "genbank")


@pytest.fixture()
def tiny_library(tmp_path: Path) -> dict[str, Any]:
    """
    Build a tiny library with **three** plasmids:
      P1: promoter J23100 + RBS B0032m + CDS ecPanD + terminator B0015
      P2: promoter J23100 + RBS B0033m + CDS ecPanD + terminator B0015
      P3: promoter J23100 + RBS B0034m + CDS ecPanD + terminator B0015 (no matching read)
    And a FASTA with two reads that uniquely match P1 and P2.
    """
    plasmids_dir = tmp_path / "gb"
    plasmids_dir.mkdir()

    # segment lengths
    lp, lr, lc, lt = 50, 30, 150, 40

    # basic building blocks
    promoter = ("ATG" * 30)[:lp]
    rbs_32 = "A" * lr
    rbs_33 = "C" * lr
    rbs_34 = "G" * lr  # P3 only
    cds = ("ATGC" * 100)[:lc]
    term = ("GC" * 100)[:lt]

    # sequences
    seq_p1 = promoter + rbs_32 + cds + term
    seq_p2 = promoter + rbs_33 + cds + term
    seq_p3 = promoter + rbs_34 + cds + term

    def build_record(pid: str, seq: str, rbs_label: str) -> SeqRecord:
        rec = SeqRecord(Seq(seq), id=pid, name=pid, description=pid)
        # features (indexes in 0-based half-open coordinates)
        rec.features.append(_mk_feature(0, lp, "promoter", "J23100"))
        rec.features.append(_mk_feature(lp, lp + lr, "RBS", rbs_label))
        rec.features.append(_mk_feature(lp + lr, lp + lr + lc, "CDS", "ecPanD"))
        rec.features.append(
            _mk_feature(lp + lr + lc, lp + lr + lc + lt, "terminator", "B0015")
        )
        return rec

    # write plasmids
    _write_gb(plasmids_dir / "P1.gb", build_record("P1", seq_p1, "B0032m"))
    _write_gb(plasmids_dir / "P2.gb", build_record("P2", seq_p2, "B0033m"))
    _write_gb(plasmids_dir / "P3.gb", build_record("P3", seq_p3, "B0034m"))

    # reads: slice around promoter+RBS+start-CDS so that RBS disambiguates P1 vs P2
    reads_fa = tmp_path / "reads.fasta"
    read1 = seq_p1[5 : lp + lr + 60]  # matches P1
    read2 = seq_p2[3 : lp + lr + 55]  # matches P2
    read3 = "TTT" * 60  # no match anywhere

    with reads_fa.open("w") as fh:
        fh.write(">EF_read_P1\n" + read1 + "\n")
        fh.write(">EF_read_P2\n" + read2 + "\n")
        fh.write(">EF_nomatch\n" + read3 + "\n")

    return {
        "plasmids_dir": plasmids_dir,
        "reads_path": reads_fa,
        "expected": {"EF_read_P1": "P1", "EF_read_P2": "P2"},
    }


# =========================
# run_3g_batch_safe (Golden Gate -> PCR -> Gibson orchestrator)
# =========================
#
# These exercise the real PCR (simulate_pcr / simulate_pcr_overhangs) and Gibson
# (merge_with_gibson_features) stages end-to-end. Only the Golden Gate stage
# (DNAcauldron-backed generate_recipe_assemblies/organize_assembly_reports) is
# monkeypatched, since it has heavy external dependencies and is orthogonal to
# the orchestration logic under test here. `_concat_records_with_features` is
# monkeypatched with a plain string-concat stand-in, same as
# test_merge_with_gibson_features_linear_with_feature_preservation above.

# Junction/middle blocks for a 2-fragment (Backbone + TU1) circular assembly.
# backbone_frag = J2 + BB_MID + J1 ; tu1_frag = J1 + TU_MID + J2
# so merging them circularly reproduces J2 + BB_MID + J1 + TU_MID.
_J1 = "TGCAGGTACCTGACTGGA"
_J2 = "CATGGAACCTTGGACCAT"
_BB_MID = "AGCTTAGGTCAACGGTTCAGATCCGA"
_TU_MID = "GATCCAGTTGGACTAGCTTGACCAGT"
_MIN_ANNEAL = 15
_MIN_OVERLAP = 15


def _install_fake_concat(monkeypatch: pytest.MonkeyPatch) -> None:
    def _concat(parts: list[SeqRecord], *, new_id: str) -> SeqRecord:
        seq = "".join(str(p.seq) for p in parts)
        rec = SeqRecord(Seq(seq), id=new_id, name=new_id)
        rec.annotations["molecule_type"] = "DNA"
        rec.features = []
        return rec

    monkeypatch.setattr(pcr_mod, "_concat_records_with_features", _concat, raising=True)


def _install_fake_golden_gate(
    monkeypatch: pytest.MonkeyPatch, *, raise_exc: Exception | None = None
) -> None:
    def _fake_generate(*, category_dirs, designs, category_order, output_dir):
        if raise_exc is not None:
            raise raise_exc
        return []

    def _fake_organize(*, report_dir, reports, delete_zip, final_only):
        return None

    monkeypatch.setattr(
        pcr_mod, "generate_recipe_assemblies", _fake_generate, raising=True
    )
    monkeypatch.setattr(
        pcr_mod, "organize_assembly_reports", _fake_organize, raising=True
    )


def _write_tu1_final(reports_dir: Path, seq: str) -> None:
    assembly_dir = reports_dir / "Assembly"
    assembly_dir.mkdir(parents=True, exist_ok=True)
    _write_gb(
        assembly_dir / "TU1_final.gb",
        _mk_record(seq, rid="TU1_final", rname="TU1_final"),
    )


def _three_g_dataframes(
    *, include_overhangs_backbone: bool, tu1_fwd: str, tu1_rev: str
):
    pcr_plan_df = pd.DataFrame(
        [
            {
                "ConstructID": "C01",
                "Role": "Backbone",
                "FwdPrimer": _J2 + _BB_MID[:_MIN_ANNEAL],
                "RevPrimer": str(
                    Seq(_BB_MID[-_MIN_ANNEAL:] + _J1).reverse_complement()
                ),
                "Circular": False,
                "IncludeOverhangs": include_overhangs_backbone,
                "MinAnneal": _MIN_ANNEAL,
            },
            {
                "ConstructID": "C01",
                "Role": "TU1",
                "FwdPrimer": tu1_fwd,
                "RevPrimer": tu1_rev,
                "Circular": False,
                "IncludeOverhangs": False,
                "MinAnneal": _MIN_ANNEAL,
            },
        ]
    )
    assembly_df = pd.DataFrame(
        [
            {"ConstructID": "C01", "Order": 1, "FragmentRole": "Backbone"},
            {"ConstructID": "C01", "Order": 2, "FragmentRole": "TU1"},
        ]
    )
    constructs_df = pd.DataFrame(
        [{"ConstructID": "C01", "MinOverlap": _MIN_OVERLAP, "Circularize": True}]
    )
    return pcr_plan_df, assembly_df, constructs_df


def _template_lookup_factory(reports_dir: Path, bb_template: SeqRecord):
    def _lookup(cid: str, role: str) -> SeqRecord:
        if role.lower() == "backbone":
            return bb_template
        path = reports_dir / "Assembly" / "TU1_final.gb"
        return SeqIO.read(str(path), "genbank")

    return _lookup


@pytest.mark.skipif(run_3g_batch_safe is None, reason="run_3g_batch_safe not exported")
def test_run_3g_batch_safe_happy_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Golden Gate (mocked) -> real PCR (with + without overhangs) -> real Gibson."""
    _install_fake_concat(monkeypatch)
    _install_fake_golden_gate(monkeypatch)

    reports_dir = tmp_path / "reports"
    tu1_frag = _J1 + _TU_MID + _J2
    _write_tu1_final(reports_dir, tu1_frag)

    bb_template = _mk_record(_BB_MID, rid="pVector", rname="pVector")
    tu1_fwd = tu1_frag[:_MIN_ANNEAL]
    tu1_rev = str(Seq(tu1_frag[-_MIN_ANNEAL:]).reverse_complement())
    pcr_plan_df, assembly_df, constructs_df = _three_g_dataframes(
        include_overhangs_backbone=True, tu1_fwd=tu1_fwd, tu1_rev=tu1_rev
    )

    result = run_3g_batch_safe(
        category_dirs={},
        designs=[],
        category_order=[],
        reports_dir=reports_dir,
        pcr_plan_df=pcr_plan_df,
        assembly_df=assembly_df,
        tus_df=pd.DataFrame(),
        constructs_df=constructs_df,
        template_lookup=_template_lookup_factory(reports_dir, bb_template),
    )

    status = result.status_df.loc["C01"]
    assert status["gg_ok"] and status["pcr_ok"] and status["gibson_ok"]
    assert status["final_bp"] == len(_J2 + _BB_MID + _J1 + _TU_MID)

    assert ("C01", "Backbone") in result.products
    assert ("C01", "TU1") in result.products

    final = result.finals["C01"]
    assert str(final.seq) == _J2 + _BB_MID + _J1 + _TU_MID

    final_path = reports_dir / "Finals" / "C01_final.gb"
    assert final_path.exists()


@pytest.mark.skipif(run_3g_batch_safe is None, reason="run_3g_batch_safe not exported")
def test_run_3g_batch_safe_golden_gate_failure_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A Golden Gate exception is wrapped in a RuntimeError, not swallowed."""
    _install_fake_golden_gate(monkeypatch, raise_exc=RuntimeError("DNAcauldron boom"))

    reports_dir = tmp_path / "reports"
    pcr_plan_df, assembly_df, constructs_df = _three_g_dataframes(
        include_overhangs_backbone=True, tu1_fwd="AAAA", tu1_rev="TTTT"
    )

    with pytest.raises(RuntimeError, match="Golden Gate stage failed"):
        run_3g_batch_safe(
            category_dirs={},
            designs=[],
            category_order=[],
            reports_dir=reports_dir,
            pcr_plan_df=pcr_plan_df,
            assembly_df=assembly_df,
            tus_df=pd.DataFrame(),
            constructs_df=constructs_df,
            template_lookup=lambda cid, role: None,
        )


@pytest.mark.skipif(run_3g_batch_safe is None, reason="run_3g_batch_safe not exported")
def test_run_3g_batch_safe_no_tu_finals_marks_construct_failed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No TU finals in reports/Assembly -> gg_ok False, PCR/Gibson skipped entirely."""
    _install_fake_concat(monkeypatch)
    _install_fake_golden_gate(monkeypatch)  # writes nothing into reports/Assembly

    reports_dir = tmp_path / "reports"
    pcr_plan_df, assembly_df, constructs_df = _three_g_dataframes(
        include_overhangs_backbone=True, tu1_fwd="AAAA", tu1_rev="TTTT"
    )

    def _lookup(cid: str, role: str) -> SeqRecord:
        raise AssertionError("template_lookup must not be called when gg_ok is False")

    result = run_3g_batch_safe(
        category_dirs={},
        designs=[],
        category_order=[],
        reports_dir=reports_dir,
        pcr_plan_df=pcr_plan_df,
        assembly_df=assembly_df,
        tus_df=pd.DataFrame(),
        constructs_df=constructs_df,
        template_lookup=_lookup,
    )

    status = result.status_df.loc["C01"]
    assert not status["gg_ok"]
    assert not status["pcr_ok"]
    assert not status["gibson_ok"]
    assert "No TU finals found" in status["gg_msg"]
    assert result.products == {}
    assert result.finals == {}


@pytest.mark.skipif(run_3g_batch_safe is None, reason="run_3g_batch_safe not exported")
def test_run_3g_batch_safe_pcr_failure_skips_gibson(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A bad primer for one role fails PCR for the whole construct; Gibson is skipped."""
    _install_fake_concat(monkeypatch)
    _install_fake_golden_gate(monkeypatch)

    reports_dir = tmp_path / "reports"
    tu1_frag = _J1 + _TU_MID + _J2
    _write_tu1_final(reports_dir, tu1_frag)

    bb_template = _mk_record(_BB_MID, rid="pVector", rname="pVector")
    # Wrong primers: won't be found in the TU1 template -> simulate_pcr raises.
    pcr_plan_df, assembly_df, constructs_df = _three_g_dataframes(
        include_overhangs_backbone=True,
        tu1_fwd="GGGGGGGGGGGGGGG",
        tu1_rev="CCCCCCCCCCCCCCC",
    )

    result = run_3g_batch_safe(
        category_dirs={},
        designs=[],
        category_order=[],
        reports_dir=reports_dir,
        pcr_plan_df=pcr_plan_df,
        assembly_df=assembly_df,
        tus_df=pd.DataFrame(),
        constructs_df=constructs_df,
        template_lookup=_template_lookup_factory(reports_dir, bb_template),
    )

    status = result.status_df.loc["C01"]
    assert status["gg_ok"]
    assert not status["pcr_ok"]
    assert "TU1" in status["pcr_msg"]
    assert not status["gibson_ok"]

    assert ("C01", "Backbone") in result.products  # backbone PCR still succeeded
    assert ("C01", "TU1") not in result.products
    assert result.finals == {}
