"""Functional group classes, plus naming and formula generation for them.

Every species handled by this program is modelled as a pair:

    (functional group class, number of carbon atoms)

That is enough to name and draw the straight-chain compounds that come up in
AQA A-level organic synthesis questions, and it is what the reaction map in
``reactions.py`` connects together.
"""

from __future__ import annotations

from dataclasses import dataclass

# --------------------------------------------------------------------------
# Functional group classes
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Group:
    """A class of organic compound (the node type in the reaction map)."""

    id: str
    label: str
    #: Short reminder of what the group looks like.
    display: str
    #: AQA specification section this group is introduced in.
    spec: str
    aromatic: bool = False


GROUPS: dict[str, Group] = {
    g.id: g
    for g in [
        Group("alkane", "alkane", "C-C single bonds only", "3.3.2"),
        Group("alkene", "alkene", "C=C", "3.3.4"),
        Group("haloalkane", "halogenoalkane", "C-X (X = Cl, Br, I)", "3.3.3"),
        Group("dihaloalkane", "dihalogenoalkane", "X-C-C-X on adjacent carbons", "3.3.4"),
        Group("alcohol_primary", "primary alcohol", "R-CH2-OH", "3.3.5"),
        Group("alcohol_secondary", "secondary alcohol", "R2CH-OH", "3.3.5"),
        Group("alcohol_tertiary", "tertiary alcohol", "R3C-OH", "3.3.5"),
        Group("aldehyde", "aldehyde", "R-CHO", "3.3.6"),
        Group("ketone", "ketone", "R-CO-R", "3.3.6"),
        Group("carboxylic_acid", "carboxylic acid", "R-COOH", "3.3.10"),
        Group("carboxylate_salt", "carboxylate salt", "R-COO- Na+", "3.3.10"),
        Group("ester", "ester", "R-COO-R", "3.3.10"),
        Group("acyl_chloride", "acyl chloride", "R-COCl", "3.3.10"),
        Group("acid_anhydride", "acid anhydride", "R-CO-O-CO-R", "3.3.10"),
        Group("amide", "amide", "R-CONH2", "3.3.10"),
        Group("n_substituted_amide", "N-substituted amide", "R-CONHR", "3.3.10"),
        Group("nitrile", "nitrile", "R-C(triple bond)N", "3.3.9"),
        Group("hydroxynitrile", "hydroxynitrile", "R-CH(OH)-CN", "3.3.9"),
        Group("hydroxy_acid", "hydroxycarboxylic acid", "R-CH(OH)-COOH", "3.3.9"),
        Group("amine_primary", "primary amine", "R-NH2", "3.3.11"),
        Group("amine_secondary", "secondary amine", "R2NH", "3.3.11"),
        Group("amine_tertiary", "tertiary amine", "R3N", "3.3.11"),
        Group("quaternary_ammonium", "quaternary ammonium salt", "R4N+ X-", "3.3.11"),
        Group("amino_acid", "2-amino acid", "R-CH(NH2)-COOH", "3.3.12"),
        Group("addition_polymer", "addition polymer", "repeating -CH2-CHR- units", "3.3.4"),
        Group("benzene", "benzene", "C6H6", "3.3.7", aromatic=True),
        Group("nitrobenzene", "nitrobenzene", "C6H5NO2", "3.3.8", aromatic=True),
        Group("phenylamine", "phenylamine", "C6H5NH2", "3.3.11", aromatic=True),
        Group("alkylbenzene", "alkylbenzene", "C6H5-R", "3.3.7", aromatic=True),
        Group("aromatic_ketone", "aromatic ketone", "C6H5-CO-R", "3.3.7", aromatic=True),
    ]
}


@dataclass(frozen=True)
class Species:
    """A specific compound: a functional group class on a carbon skeleton.

    ``carbons`` is the total number of carbon atoms, with two exceptions:

    * an acid anhydride counts one of its two identical acyl groups, so
      ethanoic anhydride is 2, matching the name;
    * ``None`` means the user named a whole class of compound ("a primary
      alcohol") rather than one molecule.
    """

    group: str
    carbons: int | None = None
    #: Esters only: how many of those carbons came from the acid.  Ethyl
    #: propanoate is ``carbons=5, acyl=3``.  ``None`` means "either half".
    acyl: int | None = None

    @property
    def info(self) -> Group:
        return GROUPS[self.group]

    @property
    def name(self) -> str:
        return name_of(self.group, self.carbons, self.acyl)

    @property
    def formula(self) -> str:
        return formula_of(self.group, self.carbons, self.acyl)

    def __str__(self) -> str:  # pragma: no cover - convenience only
        return self.name


# --------------------------------------------------------------------------
# Naming
# --------------------------------------------------------------------------

STEMS = [
    "meth", "eth", "prop", "but", "pent", "hex",
    "hept", "oct", "non", "dec", "undec", "dodec",
]


def stem(carbons: int) -> str:
    """IUPAC stem for a chain of ``carbons`` carbon atoms."""
    if 1 <= carbons <= len(STEMS):
        return STEMS[carbons - 1]
    return f"C{carbons}-alk"


def _generic(group: str) -> str:
    return f"a {GROUPS[group].label}"


def name_of(group: str, carbons: int | None, acyl: int | None = None) -> str:
    """Systematic name for ``(group, carbons)``.

    Falls back to the name of the class when the carbon count is unknown, or
    when an ester's two halves are not known separately.
    """
    if carbons is None:
        return _generic(group)

    if group == "ester":
        if acyl is None or not 1 <= acyl < carbons:
            return f"an ester (C{carbons} in total)"
        return f"{stem(carbons - acyl)}yl {stem(acyl)}anoate"

    s = stem(carbons)

    if group == "alkane":
        return f"{s}ane"
    if group == "alkene":
        if carbons < 2:
            return _generic(group)
        return f"{s}ene" if carbons <= 3 else f"{s}-1-ene"
    if group == "haloalkane":
        return f"bromo{s}ane" if carbons == 1 else f"1-bromo{s}ane"
    if group == "dihaloalkane":
        return f"1,2-dibromo{s}ane"
    if group == "alcohol_primary":
        return f"{s}anol" if carbons <= 2 else f"{s}an-1-ol"
    if group == "alcohol_secondary":
        if carbons < 3:
            return _generic(group)
        return f"{s}an-2-ol"
    if group == "alcohol_tertiary":
        if carbons == 4:
            return "2-methylpropan-2-ol"
        if carbons < 4:
            return _generic(group)
        return f"2-methyl{stem(carbons - 1)}an-2-ol"
    if group == "aldehyde":
        return f"{s}anal"
    if group == "ketone":
        if carbons < 3:
            return _generic(group)
        return "propanone" if carbons == 3 else f"{s}an-2-one"
    if group == "carboxylic_acid":
        return f"{s}anoic acid"
    if group == "carboxylate_salt":
        return f"sodium {s}anoate"
    if group == "acyl_chloride":
        return f"{s}anoyl chloride"
    if group == "acid_anhydride":
        return f"{s}anoic anhydride"
    if group == "amide":
        return f"{s}anamide"
    if group == "n_substituted_amide":
        return f"an N-substituted {s}anamide"
    if group == "nitrile":
        return "hydrogen cyanide" if carbons == 1 else f"{s}anenitrile"
    if group == "hydroxynitrile":
        return f"2-hydroxy{s}anenitrile"
    if group == "hydroxy_acid":
        return f"2-hydroxy{s}anoic acid"
    if group == "amine_primary":
        return f"{s}anamine" if carbons <= 2 else f"{s}an-1-amine"
    if group in ("amine_secondary", "amine_tertiary", "quaternary_ammonium"):
        return f"{_generic(group)} (C{carbons})"
    if group == "amino_acid":
        return f"2-amino{s}anoic acid"
    if group == "addition_polymer":
        return f"poly({name_of('alkene', carbons)})"
    if group == "benzene":
        return "benzene"
    if group == "nitrobenzene":
        return "nitrobenzene"
    if group == "phenylamine":
        return "phenylamine"
    if group == "alkylbenzene":
        side = carbons - 6
        return f"{stem(side)}ylbenzene" if side >= 1 else "benzene"
    if group == "aromatic_ketone":
        side = carbons - 6
        return f"phenyl{stem(side)}anone" if side >= 1 else _generic(group)
    return _generic(group)


# --------------------------------------------------------------------------
# Condensed structural formulae
# --------------------------------------------------------------------------


def alkyl(carbons: int) -> str:
    """Condensed formula for an alkyl group with ``carbons`` carbon atoms.

    ``alkyl(0)`` is a hydrogen atom, so that ``alkyl(n - 1) + "COOH"`` gives
    HCOOH for methanoic acid and CH3COOH for ethanoic acid.
    """
    if carbons <= 0:
        return "H"
    if carbons == 1:
        return "CH3"
    middle = carbons - 1
    if middle >= 4:
        return f"CH3(CH2){middle}"
    return "CH3" + "CH2" * middle


def alkyl_rev(carbons: int) -> str:
    """An alkyl group written for the right-hand side of a bond.

    ``alkyl(2)`` is "CH3CH2" (as in CH3CH2Br); after an oxygen the same group
    is written the other way round, "OCH2CH3", which is what this returns.
    """
    if carbons <= 0:
        return "H"
    if carbons == 1:
        return "CH3"
    middle = carbons - 1
    if middle >= 4:
        return f"(CH2){middle}CH3"
    return "CH2" * middle + "CH3"


def _tidy(formula: str) -> str:
    """Fix the leading carbon when the alkyl group is just a hydrogen atom.

    ``alkyl(0)`` is "H", so building a formula from it gives HCH2Br or HCH=CH2;
    written properly those are CH3Br and CH2=CH2.  HCHO (methanal) is already
    correct and is left alone.
    """
    if formula.startswith("HCH2"):
        return "CH3" + formula[4:]
    if formula.startswith("HCH") and formula != "HCHO":
        return "CH2" + formula[3:]
    return formula


def formula_of(group: str, carbons: int | None, acyl: int | None = None) -> str:
    """Condensed structural formula, or "" when one cannot be written."""
    if carbons is None:
        return ""

    if group == "ester":
        if acyl is None or not 1 <= acyl < carbons:
            return ""
        return _tidy(alkyl(acyl - 1) + "COO" + alkyl_rev(carbons - acyl))
    if group == "acid_anhydride":
        return _tidy(alkyl(carbons - 1) + "COOCO" + alkyl_rev(carbons - 1))

    r = alkyl(carbons - 1)

    if group == "alkane":
        return "CH4" if carbons == 1 else _tidy(alkyl(carbons - 1) + "CH3")
    if group == "alkene":
        if carbons < 2:
            return ""
        return _tidy(alkyl(carbons - 2) + "CH=CH2")
    if group == "haloalkane":
        return _tidy(r + "CH2Br")
    if group == "dihaloalkane":
        if carbons < 2:
            return ""
        return _tidy(alkyl(carbons - 2) + "CHBrCH2Br")
    if group == "alcohol_primary":
        return _tidy(r + "CH2OH")
    if group == "alcohol_secondary":
        if carbons < 3:
            return ""
        return _tidy(alkyl(carbons - 2) + "CH(OH)CH3")
    if group == "alcohol_tertiary":
        if carbons < 4:
            return ""
        return _tidy(alkyl(carbons - 3) + "C(CH3)2OH")
    if group == "aldehyde":
        return _tidy(r + "CHO")
    if group == "ketone":
        if carbons < 3:
            return ""
        return _tidy(alkyl(carbons - 2) + "COCH3")
    if group == "carboxylic_acid":
        return _tidy(r + "COOH")
    if group == "carboxylate_salt":
        return _tidy(r + "COO-Na+")
    if group == "acyl_chloride":
        return _tidy(r + "COCl")
    if group == "amide":
        return _tidy(r + "CONH2")
    if group == "nitrile":
        return _tidy(r + "CN")
    if group == "hydroxynitrile":
        if carbons < 2:
            return ""
        return _tidy(alkyl(carbons - 2) + "CH(OH)CN")
    if group == "hydroxy_acid":
        if carbons < 2:
            return ""
        return _tidy(alkyl(carbons - 2) + "CH(OH)COOH")
    if group == "amine_primary":
        return _tidy(r + "CH2NH2")
    if group == "amino_acid":
        if carbons < 2:
            return ""
        return _tidy(alkyl(carbons - 2) + "CH(NH2)COOH")
    if group == "benzene":
        return "C6H6"
    if group == "nitrobenzene":
        return "C6H5NO2"
    if group == "phenylamine":
        return "C6H5NH2"
    if group == "alkylbenzene":
        side = carbons - 6
        return f"C6H5{alkyl(side)}" if side >= 1 else "C6H6"
    if group == "aromatic_ketone":
        side = carbons - 6
        return f"C6H5CO{alkyl(side - 1)}" if side >= 1 else ""
    return ""
