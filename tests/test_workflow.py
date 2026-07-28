import sys
import types
from pathlib import Path

import pytest

NETWORKX_AVAILABLE = True
try:
    import networkx as _nx  # noqa: F401
except Exception:
    NETWORKX_AVAILABLE = False
    sys.modules["networkx"] = types.SimpleNamespace()

from assembly_designer.workflow.workflow import LiquidHandlingWorkflow  # noqa: E402


def test_add_and_list_steps():
    wf = LiquidHandlingWorkflow()
    assert wf.list_steps() == []

    wf.add_step(
        process_label="Worklist_MM",
        source_plate="src",
        destination_plate="dst",
        liquid_handling="dispense",
        sample_handling="sample_distribution",
        start="Yes",
        additional_liquid_factor=1.2,
    )

    steps = wf.list_steps()
    assert len(steps) == 1

    step = steps[0]
    assert step["process_label"] == "Worklist_MM"
    assert step["source_plate"] == "src"
    assert step["destination_plate"] == "dst"
    assert step["liquid_handling"] == "dispense"
    assert step["sample_handling"] == "sample_distribution"
    assert step["start"] == "Yes"
    assert step["additional_liquid_factor"] == 1.2


def test_add_step_defaults():
    wf = LiquidHandlingWorkflow()
    wf.add_step(process_label="Step1")

    step = wf.list_steps()[0]
    assert step["source_plate"] == "hier bitte Info"
    assert step["destination_plate"] == "hier bitte Info"
    assert step["liquid_handling"] == "hier bitte Info"
    assert step["sample_handling"] == "transfer"
    assert step["fake_source_plate"] is None
    assert step["start"] == "No"
    assert step["additional_liquid_factor"] is None


def test_remove_step_ok_and_fail():
    wf = LiquidHandlingWorkflow()
    wf.add_step("A")
    wf.add_step("B")

    wf.remove_step("A")
    labels = [step["process_label"] for step in wf.list_steps()]
    assert labels == ["B"]

    with pytest.raises(ValueError):
        wf.remove_step("not-there")


@pytest.mark.skipif(not NETWORKX_AVAILABLE, reason="networkx not installed")
def test_generate_flowchart_saves_png_when_suffix_is_given(tmp_path):
    wf = LiquidHandlingWorkflow()
    wf.add_step("Mix1", "src1", "dst1")
    wf.add_step("Mix2", "src1", "dst1")
    wf.add_step("QC", "dst1", "dst1")

    out = tmp_path / "flow.png"
    wf.generate_flowchart(
        save_path=str(out),
        fig_size=(6, 4),
        dpi=72,
        show_plot=False,
    )

    expected = Path("results") / "flow.png"
    assert expected.exists()
    assert expected.suffix == ".png"
    assert expected.stat().st_size > 0

    expected.unlink()


@pytest.mark.skipif(not NETWORKX_AVAILABLE, reason="networkx not installed")
def test_generate_flowchart_uses_export_format_when_no_suffix_is_given(tmp_path):
    wf = LiquidHandlingWorkflow()
    wf.add_step("Mix1", "src1", "dst1")

    out_base = tmp_path / "flowchart_output"
    wf.generate_flowchart(
        save_path=str(out_base),
        export_format="pdf",
        fig_size=(6, 4),
        dpi=72,
        show_plot=False,
    )

    expected = Path("results") / "flowchart_output.pdf"
    assert expected.exists()
    assert expected.suffix == ".pdf"
    assert expected.stat().st_size > 0

    expected.unlink()


@pytest.mark.skipif(not NETWORKX_AVAILABLE, reason="networkx not installed")
def test_generate_flowchart_empty_workflow_saves_file(tmp_path):
    wf = LiquidHandlingWorkflow()
    out = tmp_path / "empty.png"

    wf.generate_flowchart(
        save_path=str(out),
        fig_size=(4, 3),
        dpi=72,
        show_plot=False,
    )

    expected = Path("results") / "empty.png"
    assert expected.exists()
    assert expected.stat().st_size > 0

    expected.unlink()


@pytest.mark.skipif(not NETWORKX_AVAILABLE, reason="networkx not installed")
def test_generate_flowchart_self_loop_only_saves_file(tmp_path):
    wf = LiquidHandlingWorkflow()
    wf.add_step(
        process_label="Worklist_PM_parts",
        source_plate="mtp_source",
        destination_plate="mtp_source",
    )

    out = tmp_path / "self_loop.png"
    wf.generate_flowchart(
        save_path=str(out),
        fig_size=(6, 4),
        dpi=72,
        show_plot=False,
    )

    expected = Path("results") / "self_loop.png"
    assert expected.exists()
    assert expected.stat().st_size > 0

    expected.unlink()


@pytest.mark.skipif(not NETWORKX_AVAILABLE, reason="networkx not installed")
def test_generate_flowchart_show_plot_true_does_not_fail_in_test_backend(tmp_path):
    wf = LiquidHandlingWorkflow()
    wf.add_step("Mix1", "src1", "dst1")

    out = tmp_path / "show_plot.png"
    wf.generate_flowchart(
        save_path=str(out),
        fig_size=(6, 4),
        dpi=72,
        show_plot=True,
    )

    expected = Path("results") / "show_plot.png"
    assert expected.exists()
    assert expected.stat().st_size > 0

    expected.unlink()
