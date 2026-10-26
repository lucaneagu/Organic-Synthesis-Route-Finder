"""Turn what a user types into a :class:`Species`.

Accepts systematic names ("propan-1-ol", "1-bromobutane"), a few old/common
names ("acetone", "toluene"), condensed formulae ("CH3CH2OH") and whole classes
of compound ("a primary alcohol").
"""

from __future__ import annotations

import difflib
import re

from .species import GROUPS, STEMS, Species, formula_of, name_of


class ParseError(ValueError):
    """Raised when the input cannot be understood."""

    def __init__(self, text: str, message: str, suggestions: list[str] | None = None):
        self.text = text
        self.suggestions = suggestions or []
        super().__init__(message)


# --------------------------------------------------------------------------
# Lookup tables
# --------------------------------------------------------------------------

#: Whole classes of compound, for when no particular molecule is meant.
CLASS_ALIASES: dict[str, str] = {
    "alkane": "alkane",
    "alkene": "alkene",
    "haloalkane": "haloalkane",
    "halogenoalkane": "haloalkane",
    "haloalkanes": "haloalkane",
    "halogenoalkanes": "haloalkane",
    "alkyl halide": "haloalkane",
    "primary halogenoalkane": "haloalkane",
    "primary haloalkane": "haloalkane",
    "secondary halogenoalkane": "haloalkane_secondary",
    "secondary haloalkane": "haloalkane_secondary",
    "dihaloalkane": "dihaloalkane",
    "dihalogenoalkane": "dihaloalkane",
    "alcohol": "alcohol_primary",
    "primary alcohol": "alcohol_primary",
    "secondary alcohol": "alcohol_secondary",
    "tertiary alcohol": "alcohol_tertiary",
    "aldehyde": "aldehyde",
    "ketone": "ketone",
    "carboxylic acid": "carboxylic_acid",
    "acid": "carboxylic_acid",
    "carboxylate": "carboxylate_salt",
    "carboxylate salt": "carboxylate_salt",
    "ester": "ester",
    "acyl chloride": "acyl_chloride",
    "acid chloride": "acyl_chloride",
    "acid anhydride": "acid_anhydride",
    "anhydride": "acid_anhydride",
    "amide": "amide",
    "primary amide": "amide",
    "n-substituted amide": "n_substituted_amide",
    "secondary amide": "n_substituted_amide",
    "nitrile": "nitrile",
    "hydroxynitrile": "hydroxynitrile",
    "hydroxycarboxylic acid": "hydroxy_acid",
    "hydroxy acid": "hydroxy_acid",
    "amine": "amine_primary",
    "primary amine": "amine_primary",
    "secondary amine": "amine_secondary",
    "tertiary amine": "amine_tertiary",
    "quaternary ammonium salt": "quaternary_ammonium",
    "quaternary ammonium ion": "quaternary_ammonium",
    "amino acid": "amino_acid",
    "2-amino acid": "amino_acid",
    "addition polymer": "addition_polymer",
    "polymer": "addition_polymer",
    "arene": "benzene",
    "alkylbenzene": "alkylbenzene",
    "aromatic ketone": "aromatic_ketone",
    "aryl ketone": "aromatic_ketone",
}

#: Old, trivial or trade names that exams still occasionally use.
NAME_ALIASES: dict[str, str] = {
    "acetone": "propanone",
    "acetaldehyde": "ethanal",
    "formaldehyde": "methanal",
    "acetic acid": "ethanoic acid",
    "formic acid": "methanoic acid",
    "acetyl chloride": "ethanoyl chloride",
    "acetic anhydride": "ethanoic anhydride",
    "acetamide": "ethanamide",
    "ethylene": "ethene",
    "propylene": "propene",
    "ethyl alcohol": "ethanol",
    "methyl alcohol": "methanol",
    "isopropanol": "propan-2-ol",
    "isopropyl alcohol": "propan-2-ol",
    "toluene": "methylbenzene",
    "aniline": "phenylamine",
    "aminobenzene": "phenylamine",
    "phenylethanone": "phenylethanone",
    "acetophenone": "phenylethanone",
    "glycine": "2-aminoethanoic acid",
    "alanine": "2-aminopropanoic acid",
    "lactic acid": "2-hydroxypropanoic acid",
}

#: Compounds people ask for that this reaction map deliberately does not cover.
UNSUPPORTED: dict[str, str] = {
    "phenol": "Phenol is on the AQA spec, but none of its reactions connect it "
              "to the aliphatic synthesis map used here.",
    "benzoic acid": "Oxidation of a benzene side chain is not on the AQA "
                    "A-level specification, so there is no route to benzoic acid.",
    "chlorobenzene": "Halogenation of benzene needs a halogen carrier and is "
                     "not part of the AQA synthesis routes.",
    "glucose": "Carbohydrate chemistry is not part of AQA organic synthesis.",
    "polyester": "Condensation polymers need two difunctional monomers, which "
                 "this one-molecule-at-a-time map cannot represent.",
    "polyamide": "Condensation polymers need two difunctional monomers, which "
                 "this one-molecule-at-a-time map cannot represent.",
    "nylon": "Condensation polymers need two difunctional monomers, which this "
             "one-molecule-at-a-time map cannot represent.",
    "terylene": "Condensation polymers need two difunctional monomers, which "
                "this one-molecule-at-a-time map cannot represent.",
}

STEM_TO_CARBONS = {s: i + 1 for i, s in enumerate(STEMS)}
_STEM_RE = "|".join(STEMS)
_HALO_RE = "chloro|bromo|iodo|fluoro"


def _stem_carbons(s: str) -> int:
    return STEM_TO_CARBONS[s]


# (pattern, group id, function from match -> carbon count, or -> a
# (total carbons, acid-half carbons) pair for esters)
_PATTERNS: list[tuple[str, str, object]] = [
    # esters:  "ethyl propanoate"
    (rf"^({_STEM_RE})yl ({_STEM_RE})anoate$", "ester",
     lambda m: (_stem_carbons(m.group(1)) + _stem_carbons(m.group(2)),
                _stem_carbons(m.group(2)))),
    (rf"^(?:sodium|potassium) ({_STEM_RE})anoate$", "carboxylate_salt",
     lambda m: _stem_carbons(m.group(1))),

    # aromatics
    (r"^benzene$", "benzene", lambda m: 6),
    (r"^nitrobenzene$", "nitrobenzene", lambda m: 6),
    (r"^phenylamine$", "phenylamine", lambda m: 6),
    (rf"^({_STEM_RE})ylbenzene$", "alkylbenzene",
     lambda m: 6 + _stem_carbons(m.group(1))),
    (rf"^phenyl({_STEM_RE})anone$", "aromatic_ketone",
     lambda m: 6 + _stem_carbons(m.group(1))),

    # substituted acids and nitriles (check before the plain versions)
    (rf"^2-amino({_STEM_RE})anoic acid$", "amino_acid",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^2-hydroxy({_STEM_RE})anenitrile$", "hydroxynitrile",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^2-hydroxy({_STEM_RE})anoic acid$", "hydroxy_acid",
     lambda m: _stem_carbons(m.group(1))),

    # acid derivatives
    (rf"^({_STEM_RE})anoic acid$", "carboxylic_acid",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})ano(?:y|i)l chloride$", "acyl_chloride",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})anoic (?:acid )?anhydride$", "acid_anhydride",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})anamide$", "amide", lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})anenitrile$", "nitrile", lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})anonitrile$", "nitrile", lambda m: _stem_carbons(m.group(1))),

    # carbonyls
    (rf"^({_STEM_RE})anal$", "aldehyde", lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})an(?:-\d+)?-?one$", "ketone",
     lambda m: _stem_carbons(m.group(1))),

    # alcohols
    (rf"^2-methyl({_STEM_RE})an-2-ol$", "alcohol_tertiary",
     lambda m: _stem_carbons(m.group(1)) + 1),
    (rf"^({_STEM_RE})an-1-ol$", "alcohol_primary",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})an-2-ol$", "alcohol_secondary",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})an-3-ol$", "alcohol_secondary",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})anol$", "alcohol_primary", lambda m: _stem_carbons(m.group(1))),

    # amines
    (rf"^di({_STEM_RE})ylamine$", "amine_secondary",
     lambda m: 2 * _stem_carbons(m.group(1))),
    (rf"^tri({_STEM_RE})ylamine$", "amine_tertiary",
     lambda m: 3 * _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})an(?:-\d+)?-?amine$", "amine_primary",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})ylamine$", "amine_primary", lambda m: _stem_carbons(m.group(1))),
    (rf"^\d+-amino({_STEM_RE})ane$", "amine_primary",
     lambda m: _stem_carbons(m.group(1))),

    # halogenoalkanes
    (rf"^(?:\d+,\d+-)?di(?:{_HALO_RE})({_STEM_RE})ane$", "dihaloalkane",
     lambda m: _stem_carbons(m.group(1))),
    # A halogen on carbon 2 or 3 is secondary; parse() turns "2-bromoethane"
    # (really bromoethane) back into the primary compound.
    (rf"^[23]-(?:{_HALO_RE})({_STEM_RE})ane$", "haloalkane_secondary",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^(?:1-)?(?:{_HALO_RE})({_STEM_RE})ane$", "haloalkane",
     lambda m: _stem_carbons(m.group(1))),

    # addition polymers: poly(ethene), polyethene
    (rf"^poly\(({_STEM_RE})(?:-\d+)?-?ene\)$", "addition_polymer",
     lambda m: _stem_carbons(m.group(1))),
    (rf"^poly({_STEM_RE})ene$", "addition_polymer",
     lambda m: _stem_carbons(m.group(1))),

    # hydrocarbons
    (rf"^({_STEM_RE})(?:-\d+)?-?ene$", "alkene", lambda m: _stem_carbons(m.group(1))),
    (rf"^({_STEM_RE})ane$", "alkane", lambda m: _stem_carbons(m.group(1))),
]

_COMPILED = [(re.compile(p), g, f) for p, g, f in _PATTERNS]


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------


def normalise(text: str) -> str:
    """Lower-case, tidy dashes and drop leading articles."""
    s = text.strip().lower()
    s = s.replace("‐", "-").replace("‑", "-")
    s = s.replace("‒", "-").replace("–", "-").replace("—", "-")
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"\s*-\s*", "-", s)
    s = re.sub(r"^(a|an|the|some)\s+", "", s)
    return s.strip()


def _formula_lookup(text: str) -> Species | None:
    """Match a condensed formula such as CH3CH2OH or CH3COOH."""
    wanted = re.sub(r"[\s()]", "", text).upper()
    if not wanted or not re.fullmatch(r"[A-Z0-9=+\-]+", wanted):
        return None
    for group in GROUPS:
        for carbons in range(1, 11):
            f = formula_of(group, carbons)
            if not f:
                continue
            if re.sub(r"[\s()]", "", f).upper() == wanted:
                return Species(group, carbons)
    return None


def _known_names() -> list[str]:
    """Every name this parser can produce, for "did you mean" suggestions."""
    names = set(CLASS_ALIASES) | set(NAME_ALIASES)
    for group in GROUPS:
        names.add(GROUPS[group].label)
        for carbons in range(1, 9):
            name = name_of(group, carbons)
            if not name.startswith("a ") and "(" not in name:
                names.add(name)
    return sorted(names)


def parse(text: str) -> Species:
    """Parse ``text`` into a :class:`Species`, or raise :class:`ParseError`."""
    if not text or not text.strip():
        raise ParseError(text, "Please give the name of a chemical.")

    s = normalise(text)

    if s in UNSUPPORTED:
        raise ParseError(text, UNSUPPORTED[s])

    # A whole class of compound, e.g. "a primary alcohol".
    if s in CLASS_ALIASES:
        return Species(CLASS_ALIASES[s], None)
    if s.endswith("s") and s[:-1] in CLASS_ALIASES:
        return Species(CLASS_ALIASES[s[:-1]], None)

    # Trivial names.
    s = NAME_ALIASES.get(s, s)
    if s in UNSUPPORTED:
        raise ParseError(text, UNSUPPORTED[s])
    if s in CLASS_ALIASES:
        return Species(CLASS_ALIASES[s], None)

    for pattern, group, carbons_of in _COMPILED:
        m = pattern.match(s)
        if m:
            result = carbons_of(m)
            if isinstance(result, tuple):
                return Species(group, result[0], result[1])
            if group == "haloalkane_secondary" and result < 3:
                group = "haloalkane"
            return Species(group, result)

    by_formula = _formula_lookup(text)
    if by_formula is not None:
        return by_formula

    suggestions = difflib.get_close_matches(s, _known_names(), n=4, cutoff=0.6)
    raise ParseError(
        text,
        f"I do not recognise '{text.strip()}'.",
        suggestions,
    )
