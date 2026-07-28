from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

import robotools
from IPython.display import Image  # for GIF return in notebooks

from assembly_designer.utils import fn_plot, fn_plot2, plot_gif

LOGGER = logging.getLogger(__name__)


class MTPManager:
    """Manage multiple MTP (labware) objects with workflow-specific templates.

    The manager provides:
    - predefined templates (PCR / Golden Gate / Gibson)
    - switching between templates
    - lookup of MTP objects by name
    - convenience plotting and history animation
    """

    def __init__(self, default_template: Optional[str] = "PCR") -> None:
        """Initialize the manager with optional default template.

        Parameters
        ----------
        default_template : str or None, optional
            Initial template to activate (``"PCR"``, ``"Golden Gate"``,
            or ``"Gibson"``). If ``None``, the manager starts empty with
            no active template.
        """
        self.templates: dict[str, dict[str, Any]] = {
            "PCR": self._create_pcr_template(),
            "Golden Gate": self._create_golden_gate_template(),
            "Gibson": self._create_gibson_template(),
        }
        if default_template is None:
            self.templates = {}
        self.current_template: Optional[str] = (
            default_template if default_template else None
        )

    # --------------------------------------------------------------------- #
    # Template management
    # --------------------------------------------------------------------- #
    def switch_template(self, template_name: str) -> None:
        """Activate a specific template by name.

        Parameters
        ----------
        template_name : str
            Name of the template (e.g., ``"PCR"``, ``"Golden Gate"``,
            ``"Gibson"``).

        Raises
        ------
        ValueError
            If the requested template does not exist.
        """
        if template_name not in self.templates:
            available = list(self.templates.keys())
            raise ValueError(
                f"Template {template_name!r} not found. Available: {available}"
            )
        self.current_template = template_name
        LOGGER.info("Switched to template: %s", template_name)

    def list_templates(self) -> list[str]:
        """Return the list of available template names."""
        return list(self.templates.keys())

    def list_mtps(self) -> list[str]:
        """List names of all MTPs in the active template.

        Returns
        -------
        list of str
            MTP names. Returns an empty list if no valid template is active.
        """
        if not self.current_template or self.current_template not in self.templates:
            return []
        return list(self.templates[self.current_template].keys())

    def get_mtp(self, name: str) -> Optional[Any]:
        """Retrieve an MTP object from the active template by name.

        Parameters
        ----------
        name : str
            Name of the MTP to retrieve.

        Returns
        -------
        object or None
            The MTP object if found; otherwise ``None``.

        Raises
        ------
        ValueError
            If no valid template is active.
        """
        if not self.current_template or self.current_template not in self.templates:
            raise ValueError("No valid template is active.")
        return self.templates[self.current_template].get(name)

    # --------------------------------------------------------------------- #
    # Predefined templates
    # --------------------------------------------------------------------- #
    def _create_pcr_template(self) -> dict[str, Any]:
        """Create the default PCR template."""
        return {
            "mtp_source": robotools.Labware(
                "mtp_source[001]",
                8,
                12,
                min_volume=5,
                max_volume=450,
                initial_volumes=0,
            ),
            "mtp_destination": robotools.Labware(
                "Destination_Plate_1[001]",
                8,
                12,
                min_volume=5,
                max_volume=150,
            ),
            "master_mix": robotools.Trough(
                "mastermix",
                8,
                12,
                min_volume=20,
                max_volume=4300,
                initial_volumes=2000,
            ),
            "primer_mix": robotools.Trough(
                "primer_mix",
                8,
                12,
                min_volume=20,
                max_volume=4300,
                initial_volumes=4300,
            ),
            "mtp_primer_stocks": robotools.Labware(
                "Primer_stocks[001]",
                8,
                12,
                min_volume=5,
                max_volume=600,
                initial_volumes=0,
            ),
            "water": robotools.Trough(
                "water",
                8,
                1,
                min_volume=20,
                max_volume=100000,
                initial_volumes=80000,
            ),
            "mtp_template": robotools.Labware(
                "mtp_template[001]",
                8,
                12,
                min_volume=5,
                max_volume=200,
                initial_volumes=0,
            ),
        }

    def _create_golden_gate_template(self) -> dict[str, Any]:
        """Create the default Golden Gate template."""
        return {
            "mtp_source": robotools.Labware(
                "mtp_source",
                8,
                12,
                min_volume=5,
                max_volume=450,
                initial_volumes=0,
            ),
            "mtp_destination": robotools.Labware(
                "Destination_Plate_1",
                8,
                12,
                min_volume=5,
                max_volume=150,
            ),
            "master_mix": robotools.Trough(
                "mastermix",
                9,
                9,
                min_volume=20,
                max_volume=4300,
                initial_volumes=2000,
            ),
            "mtp_mastermix": robotools.Labware(
                "mtp_mastermix",
                8,
                12,
                min_volume=5,
                max_volume=400,
                initial_volumes=0,
            ),
            "mtp_dna_stocks": robotools.Labware(
                "mtp_dna_stocks",
                8,
                12,
                min_volume=5,
                max_volume=400,
                initial_volumes=0,
            ),
            "water": robotools.Trough(
                "water",
                8,
                1,
                min_volume=20,
                max_volume=100000,
                initial_volumes=80000,
            ),
        }

    def _create_gibson_template(self) -> dict[str, Any]:
        """Create the default Gibson Assembly template."""
        return {
            "mtp_source": robotools.Labware(
                "mtp_source",
                8,
                12,
                min_volume=5,
                max_volume=450,
                initial_volumes=0,
            ),
            "mtp_destination": robotools.Labware(
                "Destination_Plate_1",
                8,
                12,
                min_volume=5,
                max_volume=150,
            ),
            "master_mix": robotools.Trough(
                "mastermix",
                9,
                9,
                min_volume=20,
                max_volume=4300,
                initial_volumes=2000,
            ),
            "mtp_mastermix": robotools.Labware(
                "mtp_mastermix",
                8,
                12,
                min_volume=5,
                max_volume=400,
                initial_volumes=0,
            ),
            "mtp_dna_stocks": robotools.Labware(
                "mtp_dna_stocks",
                8,
                12,
                min_volume=5,
                max_volume=400,
                initial_volumes=0,
            ),
            "water": robotools.Trough(
                "water",
                8,
                1,
                min_volume=20,
                max_volume=100000,
                initial_volumes=80000,
            ),
        }

    # --------------------------------------------------------------------- #
    # Visualization helpers
    # --------------------------------------------------------------------- #
    def display_mtp_volumes(self, list_of_mtps: Optional[Sequence[str]] = None) -> None:
        """Display volume heatmaps for MTPs in the active template.

        Parameters
        ----------
        list_of_mtps : sequence of str or None, optional
            If provided, only these MTP names are plotted. If ``None``,
            all MTPs from the active template are plotted.
        """
        try:
            if list_of_mtps is None:
                if (
                    not self.current_template
                    or self.current_template not in self.templates
                ):
                    LOGGER.warning(
                        "No active template. Switch to a valid template first."
                    )
                    return

                for name, mtp in self.templates[self.current_template].items():
                    if self._validate_mtp(mtp):
                        fn_plot(mtp.volumes, f"{mtp.name} Volumes")
                    else:
                        LOGGER.error("%s lacks required attributes.", name)
            else:
                for name in list_of_mtps:
                    mtp = self.get_mtp(name)
                    if mtp is None:
                        LOGGER.warning("MTP %r not found.", name)
                        continue
                    if self._validate_mtp(mtp):
                        fn_plot(mtp.volumes, f"{mtp.name} Volumes")
                    else:
                        LOGGER.error("%s lacks required attributes.", name)
        except Exception as err:  # pylint: disable=broad-except
            LOGGER.exception("Error displaying MTP volumes: %s", err)

    def simulate_plate_history(
        self,
        mtp: Any,
        fps: float = 0.1,
        delay_frames: int = 1,
    ) -> Image:
        """Render a GIF from an MTP's history and return it as an IPython Image.

        Parameters
        ----------
        mtp : object
            An MTP-like object with ``name`` and ``history`` attributes.
        fps : float, optional
            Frames per second of the GIF. Fractional values are supported.
        delay_frames : int, optional
            Number of duplicate frames between history steps. Rounded; minimum is 0.

        Returns
        -------
        IPython.display.Image
            The generated GIF wrapped as an Image for display in notebooks.
        """
        if fps <= 0:
            raise ValueError("fps must be > 0")

        delay_i = max(0, int(round(delay_frames)))

        # Always save to the results folder
        filename = f"{getattr(mtp, 'name', 'plate')}.gif"
        out_path = Path("results") / filename

        fp = plot_gif(
            fn_plot=fn_plot2,
            fp_out=out_path,
            data=getattr(mtp, "history", []),
            fps=fps,
            delay_frames=delay_i,
        )

        LOGGER.info("Plate history GIF written to %s", fp)
        return Image(filename=str(fp))

    # --------------------------------------------------------------------- #
    # Private helpers
    # --------------------------------------------------------------------- #
    def _validate_mtp(self, mtp: Any) -> bool:
        """Return True if `mtp` exposes minimal attributes used for plotting."""
        return hasattr(mtp, "volumes") and hasattr(mtp, "name")
