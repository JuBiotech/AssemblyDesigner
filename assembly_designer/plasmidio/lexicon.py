"""Flexible lexicon for parts (promoter/RBS/gene/term) and small string helpers."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from re import Pattern
from typing import Optional

_ONLY_DIGITS_RE: Pattern[str] = re.compile(r"[0-9]+")


def normalize_ref_id(text: str) -> str:
    """
    Extract only the digits from a sequencing/trace identifier.

    Examples
    --------
    >>> normalize_ref_id("EF73802034")
    '73802034'
    >>> normalize_ref_id("  run42_read007  ")
    '42'

    Parameters
    ----------
    text : str
        Original ID string (e.g., from the sequencing provider).

    Returns
    -------
    str
        Digit-only substring, or empty string if no digits are found.
    """
    match = _ONLY_DIGITS_RE.search(text)
    return match.group(0) if match else ""


def revcomp(seq: str) -> str:
    """
    Reverse-complement a DNA string (A/C/G/T/N, case preserved).
    """
    tbl = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(tbl)[::-1]


@dataclass(frozen=True)
class Lexicon:
    """
    Configurable parts lexicon with compiled regex patterns and alias maps.

    Attributes
    ----------
    patterns : dict[str, list[Pattern[str]]]
        Compiled regular expressions per slot: "promoter", "rbs", "gene", "term".
        Patterns should usually include `(?i)` for case-insensitive matching.
    aliases : dict[str, dict[str, str]]
        Canonicalization maps per slot: raw token → canonical token.
        Example: {"rbs": {"b0033": "B0033m", "BCD12m": "BCD12"}}
    accept_unknown_label : bool
        If True and a slot is strongly implied by context (upstream code),
        unknown labels may be accepted verbatim as tokens.
    """

    patterns: dict[str, list[Pattern[str]]] = field(default_factory=dict)
    aliases: dict[str, dict[str, str]] = field(default_factory=dict)
    accept_unknown_label: bool = False


def compile_lexicon(
    config: Mapping[str, Sequence[str]] | None = None,
    aliases: Mapping[str, Mapping[str, str]] | None = None,
    *,
    accept_unknown_label: bool = False,
) -> Lexicon:
    """
    Build a compiled `Lexicon` from regex strings and alias mappings.

    Parameters
    ----------
    config : mapping or None
        Keys must be among {"promoter", "rbs", "gene", "term"}.
        Values are lists of regex strings (ideally with `(?i)` to ignore case).
        If None, a sensible default is used.
    aliases : mapping or None
        Optional alias maps per slot: raw → canonical (case-insensitive entries recommended).
    accept_unknown_label : bool, optional
        See `Lexicon.accept_unknown_label`.

    Returns
    -------
    Lexicon
        Compiled lexicon.

    Examples
    --------
    >>> my_lex = compile_lexicon(
    ...     config={
    ...         "promoter": [r"(?i)J23\\d{3}", r"(?i)pLac"],
    ...         "rbs":      [r"(?i)BCD\\d+", r"(?i)B003\\d+m?"],
    ...         "gene":     [r"(?i)ecpand", r"(?i)cgpand", r"(?i)gfp"],
    ...         "term":     [r"(?i)B001\\d+", r"(?i)terminator"],
    ...     },
    ...     aliases={"gene": {"ecpand": "ecPanD", "cgpand": "cgPanD"}},
    ... )
    """
    base: dict[str, list[str]] = {
        "promoter": [r"(?i)J23\d{3}", r"(?i)pLac"],
        "rbs": [r"(?i)BCD\d+", r"(?i)B003\d+m?"],
        "gene": [r"(?i)ecpand", r"(?i)cgpand", r"(?i)pand"],
        "term": [r"(?i)B001\d+", r"(?i)terminator"],
    }
    if config is not None:
        # Make a copy so we don't mutate caller state.
        for slot, pats in config.items():
            base[slot] = list(pats)

    compiled: dict[str, list[Pattern[str]]] = {
        slot: [re.compile(pat) for pat in pats] for slot, pats in base.items()
    }
    alias_maps: dict[str, dict[str, str]] = {
        slot: {k: v for k, v in (aliases.get(slot, {}) if aliases else {}).items()}
        for slot in ("promoter", "rbs", "gene", "term")
    }
    return Lexicon(
        patterns=compiled,
        aliases=alias_maps,
        accept_unknown_label=accept_unknown_label,
    )


def default_lexicon() -> Lexicon:
    """
    Return the default `Lexicon` compatible with the typical MoClo setup.

    Returns
    -------
    Lexicon
        Default lexicon with common promoter/RBS/gene/term patterns and aliases.
    """
    return compile_lexicon(
        config={
            "promoter": [r"(?i)J23\d{3}", r"(?i)pLac"],
            "rbs": [r"(?i)BCD\d+", r"(?i)B003\d+m?"],
            "gene": [r"(?i)ecpand", r"(?i)cgpand", r"(?i)pand"],
            "term": [r"(?i)B001\d+", r"(?i)terminator"],
        },
        aliases={
            "gene": {"ecpand": "ecPanD", "cgpand": "cgPanD", "pand": "PanD"},
            # Example for RBS canonicalization:
            # "rbs": {"b0033": "B0033m", "b0034": "B0034m"},
        },
        accept_unknown_label=False,
    )


def _apply_alias(slot: str, token: str, lex: Lexicon) -> str:
    """
    Apply alias canonicalization for a token.

    Parameters
    ----------
    slot : {"promoter", "rbs", "gene", "term"}
        Parts slot.
    token : str
        Raw token matched by a regex.
    lex : Lexicon
        Lexicon with alias maps.

    Returns
    -------
    str
        Canonicalized token if an alias exists; otherwise the original token.
    """
    amap = lex.aliases.get(slot, {})
    # Try exact key, then lower-cased key for convenience.
    return amap.get(token, amap.get(token.lower(), token))


def match_token(name: str, slot: str, lex: Optional[Lexicon] = None) -> Optional[str]:
    """
    Match the first token for a given slot within a string.

    Parameters
    ----------
    name : str
        Free-form string, e.g., feature label, filename, or concatenated name.
    slot : {"promoter", "rbs", "gene", "term"}
        Slot to match against.
    lex : Lexicon or None
        Patterns & aliases; defaults to `default_lexicon()`.

    Returns
    -------
    str or None
        Canonical token if a pattern matches, otherwise None.

    Examples
    --------
    >>> match_token("misc_feature:J23100", "promoter")
    'J23100'
    >>> match_token("5'UTR:BCD2", "rbs")
    'BCD2'
    >>> match_token("CDS:EcPanD", "gene")
    'ecPanD'
    """
    use_lex = lex or default_lexicon()
    for pat in use_lex.patterns.get(slot, []):
        found = pat.search(name)
        if found:
            return _apply_alias(slot, found.group(0), use_lex)
    return None
