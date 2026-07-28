import pytest

# Robust import: try the top-level, then nested variants
from assembly_designer.mtp_manager import MTPManager  # standard path


def test_initial_templates():
    """By default, the PCR template should be active and available in the template list."""
    manager = MTPManager()
    assert manager.current_template == "PCR"
    templates = manager.list_templates()
    assert "PCR" in templates
    assert "Golden Gate" in templates
    assert "Gibson" in templates


def test_switch_template_valid():
    """Switching to a valid template should update current_template."""
    manager = MTPManager()
    manager.switch_template("Golden Gate")
    assert manager.current_template == "Golden Gate"


def test_switch_template_invalid():
    """Switching to an invalid template should raise ValueError."""
    manager = MTPManager()
    with pytest.raises(ValueError):
        manager.switch_template("InvalidTemplate")


def test_list_mtps_and_get_mtp():
    """Listing MTPs should return a non-empty list, and get_mtp should return an object."""
    manager = MTPManager()
    mtps = manager.list_mtps()
    assert isinstance(mtps, list)
    # There should be at least one MTP in the PCR template
    assert "mtp_source" in mtps or "mtp_source[001]" in mtps

    # Retrieve an MTP object
    mtp = manager.get_mtp(mtps[0])
    assert mtp is not None


def test_get_mtp_no_active_template():
    """If no template is active, get_mtp should raise ValueError."""
    manager = MTPManager(default_template=None)
    with pytest.raises(ValueError):
        manager.get_mtp("mtp_source")
