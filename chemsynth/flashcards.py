"""Flash cards: synthesis routes with gaps to fill in.

Two decks:

* **Past papers** - multi-step routes that appeared in AQA A-level Chemistry
  (7405) Papers 2 and 3, June 2019 to June 2024.  Each card records the paper
  and question it came from.  Only the route itself is used here, not the
  question wording.  ``from_mark_scheme`` marks the steps whose reagents and
  conditions the mark scheme itself gives; the other steps use the standard
  AQA answer for that reaction.

* **Classics** - the routes examiners keep coming back to.  These are written
  as a list of compound names and every step is looked up in the reaction map
  (:mod:`chemsynth.reactions`), so the flash cards and the route finder can
  never disagree.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .parser import parse
from .pathfinder import find_route


@dataclass(frozen=True)
class Compound:
    name: str
    formula: str = ""
    #: Other names that count as correct when typed in.
    also: tuple[str, ...] = ()


@dataclass(frozen=True)
class CardStep:
    reagents: str
    conditions: str
    rxn_type: str
    from_mark_scheme: bool = False


@dataclass(frozen=True)
class Source:
    paper: str      # "7405/2" or "7405/3"
    series: str     # as printed on the cover, e.g. "June 2024"
    question: str   # e.g. "08"


@dataclass(frozen=True)
class Card:
    id: str
    compounds: tuple[Compound, ...]
    steps: tuple[CardStep, ...]
    source: Source | None = None
    note: str = ""
    #: Classic cards only: the route finder's species for each compound, so
    #: any name the parser understands is accepted as an answer.
    species: tuple = field(default=(), compare=False)


def _c(name: str, formula: str = "", *also: str) -> Compound:
    return Compound(name, formula, tuple(also))


def _s(reagents: str, conditions: str, rxn_type: str, ms: bool = False) -> CardStep:
    return CardStep(reagents, conditions, rxn_type, ms)


def _p(paper: int, year: int, question: str) -> Source:
    return Source(f"7405/{paper}", f"June {year}", question)


# Reagent wordings used on several cards.
_NITRATION = ("Conc. HNO3 and conc. H2SO4", "About 50 °C (below 55 °C)",
              "Electrophilic substitution (nitration)")
_SN_HCL = ("Sn and conc. HCl, then NaOH(aq)", "Heat under reflux", "Reduction")
_KCN = ("KCN (or NaCN)", "In aqueous ethanol, heat under reflux",
        "Nucleophilic substitution")
_H2_NI = ("H2 (or LiAlH4 in dry ether)", "Nickel catalyst", "Reduction")


PAST_PAPER_CARDS: list[Card] = [
    Card(
        id="2019-2-01",
        source=_p(2, 2019, "01"),
        compounds=(
            _c("1,4-dibromobutane", "Br(CH2)4Br"),
            _c("hexanedinitrile", "NC(CH2)4CN", "hexane-1,6-dinitrile",
               "1,4-dicyanobutane", "hexanedinitrile"),
            _c("hexane-1,6-diamine", "H2N(CH2)6NH2", "1,6-diaminohexane"),
        ),
        steps=(
            _s(*_KCN, ms=True),
            _s("H2", "Nickel (or Pt or Pd) catalyst", "Reduction", ms=True),
        ),
        note="Two CN groups go on, one at each end, so the chain grows from "
             "4 carbons to 6.",
    ),
    Card(
        id="2019-2-08",
        source=_p(2, 2019, "08"),
        compounds=(
            _c("phenol", "C6H5OH"),
            _c("4-nitrophenol", "HOC6H4NO2"),
            _c("4-aminophenol", "HOC6H4NH2"),
            _c("paracetamol", "HOC6H4NHCOCH3",
               "n-(4-hydroxyphenyl)ethanamide", "4-acetamidophenol"),
        ),
        steps=(
            _s("Dilute HNO3", "Room temperature (phenol is more reactive than "
               "benzene)", "Electrophilic substitution (nitration)"),
            _s(*_SN_HCL),
            _s("Ethanoic anhydride (or ethanoyl chloride)", "Warm",
               "Nucleophilic addition-elimination"),
        ),
        note="The paper asked for the mechanism of step 1 and the equation for "
             "step 2. Industry uses ethanoic anhydride in step 3 because "
             "ethanoyl chloride is corrosive and gives off HCl fumes.",
    ),
    Card(
        id="2019-3-04",
        source=_p(3, 2019, "04"),
        compounds=(
            _c("ethanal", "CH3CHO", "acetaldehyde"),
            _c("2-hydroxypropanenitrile", "CH3CH(OH)CN"),
            _c("propenenitrile", "CH2=CHCN", "acrylonitrile",
               "prop-2-enenitrile"),
        ),
        steps=(
            _s("KCN, then dilute acid", "Room temperature",
               "Nucleophilic addition", ms=True),
            _s("Conc. H2SO4 (or conc. H3PO4)", "Heat (about 170 °C)",
               "Elimination (dehydration)", ms=True),
        ),
        note="Step 1 makes a racemic mixture: the planar C=O is attacked "
             "equally from either side. Passing the vapour over hot Al2O3 "
             "also works for step 2.",
    ),
    Card(
        id="2019-3-21",
        source=_p(3, 2019, "21"),
        compounds=(
            _c("methylbenzene", "C6H5CH3", "toluene"),
            _c("(chloromethyl)benzene", "C6H5CH2Cl", "chloromethylbenzene",
               "benzyl chloride", "chlorophenylmethane"),
            _c("phenylmethanol", "C6H5CH2OH", "benzyl alcohol"),
            _c("phenylmethyl ethanoate", "CH3COOCH2C6H5", "benzyl ethanoate",
               "benzyl acetate"),
        ),
        steps=(
            _s("Cl2", "UV light", "Free-radical substitution"),
            _s("NaOH(aq)", "Heat under reflux", "Nucleophilic substitution"),
            _s("Ethanoyl chloride (or ethanoic anhydride)", "Room temperature",
               "Nucleophilic addition-elimination"),
        ),
        note="A multiple-choice question: the answer was that electrophilic "
             "substitution is not used. The chlorine goes on the side chain, "
             "not the ring, because the conditions are free-radical ones.",
    ),
    Card(
        id="2020-2-06",
        source=_p(2, 2020, "06"),
        compounds=(
            _c("benzene", "C6H6"),
            _c("phenylethanone", "C6H5COCH3", "acetophenone",
               "1-phenylethanone"),
            _c("1-phenylethanol", "C6H5CH(OH)CH3", "1-phenylethan-1-ol"),
            _c("phenylethene", "C6H5CH=CH2", "styrene", "ethenylbenzene"),
        ),
        steps=(
            _s("Ethanoyl chloride with AlCl3", "Anhydrous, heat under reflux",
               "Electrophilic substitution (Friedel-Crafts acylation)", ms=True),
            _s("NaBH4", "Room temperature, aqueous methanol",
               "Nucleophilic addition (reduction)", ms=True),
            _s("Conc. H2SO4", "Heat", "Elimination (dehydration)", ms=True),
        ),
        note="In the paper, phenylethene was then polymerised to make "
             "polystyrene.",
    ),
    Card(
        id="2021-2-08",
        source=_p(2, 2021, "08"),
        compounds=(
            _c("cyclohexanol", "C6H11OH"),
            _c("cyclohexene", "C6H10"),
            _c("1,2-dibromocyclohexane", "C6H10Br2"),
            _c("cyclohexane-1,2-diol", "C6H10(OH)2"),
            _c("the cyclic diester", "C6H10(OCO)2",
               "cyclohexane-1,2-diyl ethanedioate"),
        ),
        steps=(
            _s("Conc. H2SO4 (or conc. H3PO4)", "Heat",
               "Elimination (dehydration)", ms=True),
            _s("Br2", "Room temperature", "Electrophilic addition", ms=True),
            _s("NaOH(aq)", "Heat under reflux", "Nucleophilic substitution",
               ms=True),
            _s("Ethanedioyl dichloride, (COCl)2", "Room temperature",
               "Nucleophilic addition-elimination", ms=True),
        ),
        note="Both OH groups react with the two acyl chloride groups of "
             "ethanedioyl dichloride, closing a six-membered ring.",
    ),
    Card(
        id="2021-3-30",
        source=_p(3, 2021, "30"),
        compounds=(
            _c("bromoethane", "CH3CH2Br"),
            _c("propanenitrile", "CH3CH2CN"),
            _c("propan-1-amine", "CH3CH2CH2NH2", "propylamine",
               "1-aminopropane"),
        ),
        steps=(
            _s(*_KCN),
            _s(*_H2_NI),
        ),
        note="A multiple-choice question asking for the intermediate. Making "
             "the amine via the nitrile avoids the mixture of amines you get "
             "with ammonia, and adds one carbon.",
    ),
    Card(
        id="2022-2-05",
        source=_p(2, 2022, "05"),
        compounds=(
            _c("hex-1-ene", "CH3(CH2)3CH=CH2"),
            _c("hexan-2-ol", "CH3(CH2)3CH(OH)CH3"),
        ),
        steps=(
            _s("Conc. H3PO4, then water", "Cool in an ice bath, then add water",
               "Electrophilic addition (hydration)", ms=True),
        ),
        note="Markovnikov: the OH ends up on carbon 2, so the product is the "
             "secondary alcohol, not hexan-1-ol.",
    ),
    Card(
        id="2022-2-10a",
        source=_p(2, 2022, "10"),
        compounds=(
            _c("methylbenzene", "C6H5CH3", "toluene"),
            _c("1-methyl-4-nitrobenzene", "CH3C6H4NO2", "4-nitromethylbenzene",
               "4-nitrotoluene", "methyl-4-nitrobenzene"),
            _c("4-methylphenylamine", "CH3C6H4NH2", "4-methylaniline",
               "4-aminomethylbenzene"),
            _c("N-(4-methylphenyl)ethanamide", "CH3C6H4NHCOCH3"),
        ),
        steps=(
            _s(*_NITRATION, ms=True),
            _s("Sn and conc. HCl", "Heat under reflux", "Reduction", ms=True),
            _s("Ethanoyl chloride (or ethanoic anhydride)", "Room temperature",
               "Nucleophilic addition-elimination"),
        ),
        note="The mark scheme gives steps 1 and 2 and names the mechanism for "
             "step 3.",
    ),
    Card(
        id="2022-2-10b",
        source=_p(2, 2022, "10"),
        compounds=(
            _c("methylbenzene", "C6H5CH3", "toluene"),
            _c("(chloromethyl)benzene", "C6H5CH2Cl", "chloromethylbenzene",
               "benzyl chloride", "chlorophenylmethane"),
            _c("phenylmethanamine", "C6H5CH2NH2", "benzylamine",
               "phenylmethylamine"),
        ),
        steps=(
            _s("Cl2", "UV light", "Free-radical substitution", ms=True),
            _s("Excess NH3", "Heat in a sealed tube", "Nucleophilic substitution"),
        ),
        note="The yield of this amine is lower than the nitro-reduction route, "
             "because it reacts further to give secondary and tertiary amines. "
             "It is a stronger base than an aromatic amine: its lone pair is "
             "not delocalised into the ring.",
    ),
    Card(
        id="2022-3-30",
        source=_p(3, 2022, "30"),
        compounds=(
            _c("phenylethanone", "C6H5COCH3", "acetophenone",
               "1-phenylethanone"),
            _c("1-(3-nitrophenyl)ethanone", "O2NC6H4COCH3",
               "3-nitrophenylethanone"),
            _c("1-(3-nitrophenyl)ethanol", "O2NC6H4CH(OH)CH3",
               "1-(3-nitrophenyl)ethan-1-ol"),
            _c("1-(3-nitrophenyl)ethyl ethanoate", "O2NC6H4CH(OCOCH3)CH3"),
        ),
        steps=(
            _s(*_NITRATION),
            _s("NaBH4", "Room temperature, aqueous methanol",
               "Nucleophilic addition (reduction)"),
            _s("Ethanoic anhydride (or ethanoic acid with conc. H2SO4)",
               "Warm", "Esterification"),
        ),
        note="A multiple-choice question: the answer was that hydrolysis is "
             "not involved. NaBH4 reduces the C=O but leaves the NO2 group "
             "alone.",
    ),
    Card(
        id="2023-2-04",
        source=_p(2, 2023, "04"),
        compounds=(
            _c("propenenitrile", "CH2=CHCN", "acrylonitrile",
               "prop-2-enenitrile"),
            _c("3-bromopropanenitrile", "BrCH2CH2CN"),
            _c("butanedinitrile", "NCCH2CH2CN", "butane-1,4-dinitrile",
               "1,2-dicyanoethane"),
            _c("butane-1,4-diamine", "H2N(CH2)4NH2", "1,4-diaminobutane"),
        ),
        steps=(
            _s("HBr", "Room temperature", "Electrophilic addition", ms=True),
            _s("KCN (or NaCN)", "Warm, in aqueous ethanol",
               "Nucleophilic substitution", ms=True),
            _s(*_H2_NI, ms=True),
        ),
        note="Step 1 gives a mixture: 3-bromopropanenitrile and the two "
             "optical isomers of 2-bromopropanenitrile.",
    ),
    Card(
        id="2023-3-26",
        source=_p(3, 2023, "26"),
        compounds=(
            _c("propanone", "CH3COCH3", "acetone"),
            _c("2-hydroxy-2-methylpropanenitrile", "(CH3)2C(OH)CN"),
            _c("2-methylpropenenitrile", "CH2=C(CH3)CN",
               "2-methylprop-2-enenitrile", "methacrylonitrile"),
            _c("2-methylpropenoic acid", "CH2=C(CH3)COOH",
               "2-methylprop-2-enoic acid", "methacrylic acid"),
            _c("methyl 2-methylpropenoate", "CH2=C(CH3)COOCH3",
               "methyl 2-methylprop-2-enoate", "methyl methacrylate"),
            _c("methyl 2-methylpropanoate", "(CH3)2CHCOOCH3"),
        ),
        steps=(
            _s("KCN, then dilute H2SO4", "Room temperature",
               "Nucleophilic addition"),
            _s("Conc. H2SO4", "Heat", "Elimination (dehydration)"),
            _s("Dilute HCl(aq)", "Heat under reflux", "Hydrolysis"),
            _s("Methanol, with conc. H2SO4 as catalyst", "Heat under reflux",
               "Esterification"),
            _s("H2", "Nickel catalyst", "Addition (reduction)", ms=True),
        ),
        note="A multiple-choice question asking which step had a correct "
             "reagent and condition. The answer was step 5. The wrong options "
             "were HCN in water, KOH in warm water for the dehydration, and an "
             "alkaline catalyst for the esterification.",
    ),
    Card(
        id="2023-3-31",
        source=_p(3, 2023, "31"),
        compounds=(
            _c("propene", "CH3CH=CH2", "propylene"),
            _c("propan-2-yl hydrogensulfate", "CH3CH(OSO3H)CH3",
               "2-propyl hydrogensulfate", "propyl hydrogensulfate",
               "isopropyl hydrogensulfate", "1-methylethyl hydrogensulfate",
               "an alkyl hydrogensulfate"),
            _c("propan-2-ol", "CH3CH(OH)CH3"),
        ),
        steps=(
            _s("Conc. H2SO4", "Cold", "Electrophilic addition"),
            _s("Water", "Warm", "Hydrolysis"),
        ),
        note="A multiple-choice question: propene was the compound that turns "
             "into an alcohol in two stages. Markovnikov again: the OH ends up "
             "on carbon 2.",
    ),
    Card(
        id="2024-2-08",
        source=_p(2, 2024, "08"),
        compounds=(
            _c("but-1-ene", "CH3CH2CH=CH2"),
            _c("2-bromobutane", "CH3CH2CHBrCH3"),
            _c("butan-2-ol", "CH3CH2CH(OH)CH3"),
            _c("butanone", "CH3CH2COCH3", "butan-2-one"),
            _c("2-hydroxy-2-methylbutanenitrile", "CH3CH2C(OH)(CN)CH3"),
        ),
        steps=(
            _s("HBr", "Room temperature", "Electrophilic addition", ms=True),
            _s("NaOH(aq)", "Warm", "Nucleophilic substitution", ms=True),
            _s("Acidified K2Cr2O7", "Heat under reflux", "Oxidation"),
            _s("HCN (or KCN, then dilute H2SO4)", "Room temperature",
               "Nucleophilic addition", ms=True),
        ),
        note="Step 4 gives a racemic mixture, because the planar C=O is "
             "attacked equally from above and below.",
    ),
    Card(
        id="2024-2-09",
        source=_p(2, 2024, "09"),
        compounds=(
            _c("2-bromopropane", "CH3CHBrCH3"),
            _c("2-methylpropanenitrile", "(CH3)2CHCN"),
            _c("2-methylpropan-1-amine", "(CH3)2CHCH2NH2",
               "2-methylpropylamine", "1-amino-2-methylpropane"),
        ),
        steps=(
            _s(*_KCN),
            _s("H2", "Nickel catalyst", "Reduction", ms=True),
        ),
        note="The CN replaces the Br on the middle carbon, so the chain "
             "becomes branched.",
    ),
    Card(
        id="2024-2-10",
        source=_p(2, 2024, "10"),
        compounds=(
            _c("1,2-dimethylbenzene", "C6H4(CH3)2"),
            _c("1,2-dimethyl-4-nitrobenzene", "(CH3)2C6H3NO2",
               "4-nitro-1,2-dimethylbenzene", "3,4-dimethylnitrobenzene"),
            _c("3,4-dimethylphenylamine", "(CH3)2C6H3NH2",
               "3,4-dimethylaniline"),
        ),
        steps=(
            _s(*_NITRATION, ms=True),
            _s("Sn and conc. HCl", "Heat under reflux", "Reduction", ms=True),
        ),
        note="The electrophile is NO2+, made when H2SO4 protonates HNO3: "
             "HNO3 + 2H2SO4 -> NO2+ + H3O+ + 2HSO4-.",
    ),
]


#: Classic routes, as names the parser understands.  Every neighbouring pair
#: must be exactly one reaction on the map.
CLASSIC_ROUTES: list[tuple[str, ...]] = [
    ("ethene", "bromoethane", "ethanol", "ethanal"),
    ("ethanol", "ethanoic acid", "ethyl ethanoate"),
    ("propene", "propan-2-ol", "propanone"),
    ("propan-1-ol", "propanal", "2-hydroxybutanenitrile", "2-hydroxybutanoic acid"),
    ("1-bromobutane", "pentanenitrile", "pentanoic acid", "methyl pentanoate"),
    ("propanal", "propan-1-ol", "1-bromopropane", "butanenitrile"),
    ("propanoic acid", "propan-1-ol", "propene", "1,2-dibromopropane"),
    ("propan-2-ol", "2-bromopropane", "propene", "propane"),
    ("ethyl propanoate", "sodium propanoate", "propanoic acid"),
    ("ethanal", "2-hydroxypropanenitrile", "2-hydroxypropanoic acid"),
    ("benzene", "nitrobenzene", "phenylamine"),
]


def classic_card(names: tuple[str, ...]) -> Card:
    """Build a card from the reaction map, one reaction per pair of names."""
    species = [parse(name) for name in names]
    steps = []
    for start, end in zip(species, species[1:]):
        route = find_route(start, end)
        if route.length != 1:
            raise ValueError(f"{start.name} -> {end.name} is {route.length} "
                             f"steps on the map, not 1")
        r = route.steps[0].reaction
        steps.append(CardStep(r.reagents, r.conditions, r.rxn_type))
    compounds = tuple(Compound(s.name, s.formula) for s in species)
    return Card(
        id="classic-" + "-".join(n.replace(" ", "_") for n in (names[0], names[-1])),
        compounds=compounds,
        steps=tuple(steps),
        species=tuple(species),
    )


def all_cards() -> list[Card]:
    return PAST_PAPER_CARDS + [classic_card(names) for names in CLASSIC_ROUTES]
