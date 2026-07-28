"""Restriction digest & agarose gel simulation."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .io_ import load_dna_file


def get_ladder(ladder_type: str) -> list[int]:
    """
    Return fragment sizes (in base pairs) for standard DNA ladders.

    Parameters
    ----------
    ladder_type : str
        Type of DNA ladder. Supported types: "10kb", "5kb", "2kb", "1kb", "100bp".

    Returns
    -------
    list of int
        List of fragment sizes in base pairs (bp), sorted from largest to smallest.

    Raises
    ------
    ValueError
        If an unsupported ladder type is provided.
    """
    ladders: dict[str, list[int]] = {
        "10kb": [10000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 500],
        "5kb": [5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250],
        "2kb": [2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        "1kb": [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        "100bp": [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
    }

    if ladder_type not in ladders:
        raise ValueError(
            "Invalid ladder type. Choose from '10kb', '5kb', '2kb', '1kb', or '100bp'."
        )

    return ladders[ladder_type]


def perform_restriction_digest(
    sequence: Union[str, Seq], enzymes: list[str]
) -> list[int]:
    """
    Simulate a restriction digest and return the resulting fragment sizes.

    Parameters
    ----------
    sequence : str
        DNA sequence to digest.
    enzymes : list of str
        List of enzyme names as strings (e.g., ["EcoRI", "BamHI"]).

    Returns
    -------
    list of int
        Sorted list of fragment sizes in base pairs (largest to smallest).
    """
    batch = RestrictionBatch(enzymes)
    cuts = batch.search(sequence)
    cut_positions = sorted(
        {0, len(sequence)} | {pos for positions in cuts.values() for pos in positions}
    )
    fragments = [
        cut_positions[i + 1] - cut_positions[i] for i in range(len(cut_positions) - 1)
    ]
    return sorted(fragments, reverse=True)


def read_plasmid_file(file_path: str) -> SeqRecord:
    """Read a plasmid file (GenBank or SnapGene) as a normalized SeqRecord.

    Convenience wrapper around ``load_dna_file`` for callers that pass a string
    path. Applies the same normalization rules (fixes empty/placeholder id/name
    and enforces ``molecule_type='DNA'``).

    Parameters
    ----------
    file_path : str
        Path to the input file. Supports ``.gb``, ``.gbk``, ``.genbank``, and ``.dna``.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        A Biopython ``SeqRecord`` with normalized identifiers.

    Raises
    ------
    FileNotFoundError
        If ``file_path`` does not exist.
    ValueError
        If the file extension is unsupported.
    RuntimeError
        If reading a SnapGene (``.dna``) file without ``snapgene_reader`` installed.
    """
    return load_dna_file(Path(file_path))


def plot_gel(
    fragments_by_lane: dict[str, list[int]],
    lane_labels: list[str],
    ladder_type: str = "2kb",
    save_path: Optional[str] = None,
) -> None:
    """
    Plots a simulated agarose gel for multiple lanes of DNA fragments.

    Parameters
    ----------
    fragments_by_lane : dict[str, list[int]]
        Dictionary with lane names as keys and lists of fragment sizes (bp) as values.
    lane_labels : list[str]
        Labels for each lane including the ladder.
    ladder_type : str, optional
        Type of DNA ladder to use ("10kb", "5kb", "2kb", "1kb", "100bp"), by default "2kb".
    save_path : Optional[str], optional
        File path to save the plot as an image (e.g. "gel.png"). If None, plot is not saved, by default None.

    Returns
    -------
    None
        Shows the gel plot and optionally saves it to a file.
    """
    gel_height = 12
    gel_width = len(fragments_by_lane) + 1
    band_height = 0.2  # Smaller bands for better visualization

    fig, ax = plt.subplots(figsize=(gel_width * 2, 8))
    ax.imshow(
        np.ones((gel_height * 20, gel_width * 50)),
        cmap="Greys",
        aspect="auto",
        extent=(0, gel_width, 0, gel_height),
    )

    lane_idx = 0
    ladder_fragments = get_ladder(ladder_type)
    ladder_positions = [gel_height - (np.log10(size) * 3) for size in ladder_fragments]

    # Plot ladder (Lane 1)
    for position in ladder_positions:
        ax.hlines(
            y=position,
            xmin=lane_idx + 0.3,
            xmax=lane_idx + 0.7,
            color="cornflowerblue",
            linewidth=band_height * 10,
        )
    for position, size in zip(ladder_positions, ladder_fragments, strict=False):
        ax.text(lane_idx + 0.8, position, f"{size} bp", fontsize=8, va="center")
    lane_idx += 1

    # Plot other lanes (Lane 2+)
    for _lane_name, fragments in fragments_by_lane.items():
        sorted_fragments = sorted(fragments, reverse=True)
        positions = [gel_height - (np.log10(size) * 3) for size in sorted_fragments]
        for position, size in zip(positions, sorted_fragments, strict=False):
            ax.hlines(
                y=position,
                xmin=lane_idx + 0.3,
                xmax=lane_idx + 0.7,
                color="cornflowerblue",
                linewidth=band_height * 10,
            )
            ax.text(lane_idx + 0.8, position, f"{size} bp", fontsize=8, va="center")
        lane_idx += 1

    # Add lane numbers with white bars at the bottom
    for i, label in enumerate(range(1, len(lane_labels) + 1)):
        ax.add_patch(
            plt.Rectangle((i + 0.1, 1), 1, 2, color="white")
        )  # White background bar
        ax.text(
            i + 0.5,
            10,
            str(label),
            fontsize=8,
            color="black",
            ha="center",
            va="center",
            weight="bold",
        )

    ax.set_xlim(0, gel_width)
    ax.set_ylim(-1.5, gel_height)  # Extra space for labels
    ax.invert_yaxis()  # Simulate gel migration from top to bottom
    ax.set_title("Simulated Agarose Gel", fontsize=16)
    ax.axis("off")

    # Add lane legend
    plt.legend(
        [f"{i + 1}. {label}" for i, label in enumerate(lane_labels)],
        loc="upper right",
        fontsize=8,
        title="Lanes",
    )

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Gel image saved to {save_path}")

    plt.show()


def digest_and_plot_plasmids(
    plasmid_files_with_enzymes: list[tuple[str, list[str]]],
    ladder_type: str = "2kb",
    save_path: Optional[str] = None,
) -> None:
    """
    Reads plasmids, performs restriction digests, and plots the gel.

    Parameters
    ----------
    plasmid_files_with_enzymes : list[tuple[str, list[str]]]
        List of tuples with (file_path, list of enzyme names).
    ladder_type : str, optional
        DNA ladder type, by default "2kb".
    save_path : Optional[str], optional
        File path to save the gel plot image, by default None.

    Returns
    -------
    None
    """
    fragments_by_lane = {}

    for plasmid_file, enzymes in plasmid_files_with_enzymes:
        record = read_plasmid_file(plasmid_file)
        fragments = perform_restriction_digest(record.seq, enzymes)
        fragments_by_lane[record.name] = fragments

    lane_labels = ["Ladder"] + [
        read_plasmid_file(f).name for f, _ in plasmid_files_with_enzymes
    ]

    plot_gel(
        fragments_by_lane, lane_labels, ladder_type=ladder_type, save_path=save_path
    )
