"""The AQA A-level organic reaction map.

Each :class:`Reaction` is one arrow on the synthesis map: it turns one
functional group class into another, using the reagents and conditions AQA
expect you to quote in an answer.

Carbon bookkeeping
------------------
``delta_min`` / ``delta_max`` give the change in the number of carbon atoms.
Most steps leave the chain alone (0, 0).  Steps that use a co-reactant whose
size you choose - esterification, Friedel-Crafts - have a range instead, and
``None`` means "unbounded in that direction".
"""

from __future__ import annotations

from dataclasses import dataclass

CORE = "core"          # explicitly named in the AQA specification
EXTENSION = "extension"  # useful, but check with your teacher before quoting


@dataclass(frozen=True)
class Reaction:
    id: str
    start: str
    end: str
    reagents: str
    conditions: str
    rxn_type: str
    spec: str
    delta_min: int = 0
    delta_max: int = 0
    note: str = ""
    #: Search cost.  1.0 is a clean, high-yield, exam-standard step; higher
    #: values push a route down the rankings without forbidding it.
    cost: float = 1.0
    status: str = CORE
    #: Restrictions on the carbon count of the *starting* species.
    carbons_min: int | None = None
    carbons_max: int | None = None
    #: Name of a special carbon-bookkeeping rule (see ``pathfinder``).  Used by
    #: the ester reactions, where the two halves of the molecule have to be
    #: tracked separately, and ``delta_min`` / ``delta_max`` are ignored.
    carbon_rule: str | None = None

    @property
    def fixed_delta(self) -> bool:
        return self.delta_min == self.delta_max

    def allows(self, carbons: int | None) -> bool:
        """Is this reaction available to a starting species of this size?"""
        if carbons is None:
            return True
        if self.carbons_min is not None and carbons < self.carbons_min:
            return False
        if self.carbons_max is not None and carbons > self.carbons_max:
            return False
        return True


REACTIONS: list[Reaction] = [
    # ---------------------------------------------------------------- alkanes
    Reaction(
        id="alkane_halogenation",
        start="alkane", end="haloalkane",
        reagents="Cl2 (or Br2)",
        conditions="UV light (sunlight), gas phase",
        rxn_type="Free-radical substitution",
        spec="3.3.2", cost=1.8,
        note="Gives a mixture of mono- and poly-substituted products, so the "
             "yield of any one halogenoalkane is poor. Avoid it in a synthesis "
             "if another route exists.",
    ),

    # ---------------------------------------------------------------- alkenes
    Reaction(
        id="alkene_hydrogenation",
        start="alkene", end="alkane",
        reagents="H2",
        conditions="Nickel catalyst, 150 °C",
        rxn_type="Addition (reduction)",
        spec="3.3.4",
    ),
    Reaction(
        id="alkene_hydrogen_halide",
        start="alkene", end="haloalkane",
        reagents="HBr (or HCl) gas",
        conditions="Room temperature",
        rxn_type="Electrophilic addition",
        spec="3.3.4",
        note="Markovnikov: the H adds to the carbon that already has more "
             "hydrogens, so an unsymmetrical alkene gives mainly the secondary "
             "halogenoalkane (via the more stable carbocation).",
    ),
    Reaction(
        id="alkene_bromine",
        start="alkene", end="dihaloalkane",
        reagents="Br2 (bromine water, or in an organic solvent)",
        conditions="Room temperature",
        rxn_type="Electrophilic addition",
        spec="3.3.4",
        note="Orange bromine water is decolourised - the test for a C=C bond.",
    ),
    Reaction(
        id="alkene_hydration_primary",
        start="alkene", end="alcohol_primary",
        reagents="Steam (H2O)",
        conditions="Conc. H3PO4 catalyst, 300 °C, 60 atm",
        rxn_type="Electrophilic addition (hydration)",
        spec="3.3.5", carbons_max=2,
        note="Only ethene gives a primary alcohol this way; longer alkenes "
             "follow Markovnikov and give a secondary alcohol.",
    ),
    Reaction(
        id="alkene_hydration_secondary",
        start="alkene", end="alcohol_secondary",
        reagents="Steam (H2O)",
        conditions="Conc. H3PO4 catalyst, 300 °C, 60 atm",
        rxn_type="Electrophilic addition (hydration)",
        spec="3.3.5", carbons_min=3,
        note="Markovnikov addition, so the OH ends up on the middle carbon.",
    ),
    Reaction(
        id="alkene_sulfuric_hydration",
        start="alkene", end="alcohol_secondary",
        reagents="Conc. H2SO4, then water",
        conditions="Cold for the first stage, then warm with water to hydrolyse "
                   "the alkyl hydrogensulfate",
        rxn_type="Electrophilic addition (two-stage hydration)",
        spec="3.3.4", cost=1.2, carbons_min=3,
        note="The laboratory version of hydration; the H2SO4 is regenerated.",
    ),
    Reaction(
        id="alkene_polymerisation",
        start="alkene", end="addition_polymer",
        reagents="Many alkene molecules (the monomer)",
        conditions="High pressure and a Ziegler-Natta catalyst",
        rxn_type="Addition polymerisation",
        spec="3.3.4",
        note="The C=C opens up; every atom of the monomer ends up in the polymer.",
    ),

    # --------------------------------------------------------- halogenoalkanes
    Reaction(
        id="haloalkane_hydrolysis",
        start="haloalkane", end="alcohol_primary",
        reagents="Aqueous NaOH (or KOH)",
        conditions="Heat under reflux, in aqueous solution",
        rxn_type="Nucleophilic substitution",
        spec="3.3.3",
        note="Water is a poor nucleophile, so use hydroxide ions. Aqueous "
             "conditions favour substitution over elimination.",
    ),
    Reaction(
        id="haloalkane_elimination",
        start="haloalkane", end="alkene",
        reagents="Ethanolic NaOH (or KOH)",
        conditions="Heat under reflux in ethanol - no water present",
        rxn_type="Elimination",
        spec="3.3.3",
        note="Hot ethanolic conditions make OH- act as a base rather than a "
             "nucleophile. Substitution is the competing reaction.",
    ),
    Reaction(
        id="haloalkane_cyanide",
        start="haloalkane", end="nitrile",
        reagents="KCN (or NaCN) dissolved in ethanol",
        conditions="Heat under reflux",
        rxn_type="Nucleophilic substitution",
        spec="3.3.3", delta_min=1, delta_max=1,
        note="The standard way to LENGTHEN a carbon chain by one: the CN carbon "
             "becomes part of the molecule.",
    ),
    Reaction(
        id="haloalkane_ammonia",
        start="haloalkane", end="amine_primary",
        reagents="Excess ethanolic ammonia (NH3)",
        conditions="Heat in a sealed tube under pressure",
        rxn_type="Nucleophilic substitution",
        spec="3.3.11",
        note="Use a large excess of ammonia, otherwise the primary amine reacts "
             "again to give secondary and tertiary amines and quaternary salts.",
    ),

    # ---------------------------------------------------------------- alcohols
    Reaction(
        id="alcohol1_to_aldehyde",
        start="alcohol_primary", end="aldehyde",
        reagents="Acidified potassium dichromate(VI) (K2Cr2O7 / dilute H2SO4)",
        conditions="Warm gently and DISTIL the aldehyde off as it forms",
        rxn_type="Oxidation",
        spec="3.3.5",
        note="Distilling immediately removes the aldehyde before it can be "
             "oxidised further. Orange dichromate turns green.",
    ),
    Reaction(
        id="alcohol1_to_acid",
        start="alcohol_primary", end="carboxylic_acid",
        reagents="Excess acidified potassium dichromate(VI)",
        conditions="Heat under reflux",
        rxn_type="Oxidation",
        spec="3.3.5",
        note="Reflux with excess oxidising agent takes the oxidation all the "
             "way to the acid.",
    ),
    Reaction(
        id="alcohol2_to_ketone",
        start="alcohol_secondary", end="ketone",
        reagents="Acidified potassium dichromate(VI)",
        conditions="Heat under reflux (or distil the ketone off)",
        rxn_type="Oxidation",
        spec="3.3.5",
        note="Ketones resist further oxidation, so the conditions are not "
             "critical here.",
    ),
    Reaction(
        id="alcohol1_dehydration",
        start="alcohol_primary", end="alkene",
        reagents="Conc. H2SO4 (or conc. H3PO4) catalyst",
        conditions="Heat under reflux at about 170 °C",
        rxn_type="Elimination (dehydration)",
        spec="3.3.5",
        note="Also works by passing the alcohol vapour over hot Al2O3.",
    ),
    Reaction(
        id="alcohol2_dehydration",
        start="alcohol_secondary", end="alkene",
        reagents="Conc. H2SO4 (or conc. H3PO4) catalyst",
        conditions="Heat under reflux at about 170 °C",
        rxn_type="Elimination (dehydration)",
        spec="3.3.5",
        note="A secondary alcohol can give more than one alkene - watch for "
             "E/Z isomers in the products.",
    ),
    Reaction(
        id="alcohol1_to_haloalkane",
        start="alcohol_primary", end="haloalkane",
        reagents="NaBr with conc. H2SO4 (makes HBr in situ)",
        conditions="Heat under reflux",
        rxn_type="Nucleophilic substitution",
        spec="3.3.3",
        note="For a chloroalkane use PCl5, PCl3 or SOCl2; for an iodoalkane use "
             "red phosphorus with iodine.",
    ),
    Reaction(
        id="alcohol2_to_haloalkane",
        start="alcohol_secondary", end="haloalkane",
        reagents="NaBr with conc. H2SO4 (makes HBr in situ)",
        conditions="Heat under reflux",
        rxn_type="Nucleophilic substitution",
        spec="3.3.3",
    ),
    Reaction(
        id="alcohol_esterification",
        start="alcohol_primary", end="ester",
        reagents="A carboxylic acid, with a few drops of conc. H2SO4 as catalyst",
        conditions="Warm / heat under reflux",
        rxn_type="Esterification (condensation)",
        spec="3.3.10", carbon_rule="ester_from_alcohol", cost=1.1,
        note="Reversible and slow, so the yield is limited. Sweet-smelling "
             "product.",
    ),
    Reaction(
        id="alcohol_plus_acyl_chloride",
        start="alcohol_primary", end="ester",
        reagents="An acyl chloride (e.g. ethanoyl chloride)",
        conditions="Room temperature",
        rxn_type="Addition-elimination",
        spec="3.3.10", carbon_rule="ester_from_alcohol",
        note="Much faster and higher-yielding than esterification, and it goes "
             "to completion. Misty fumes of HCl are given off.",
    ),

    # --------------------------------------------------- aldehydes and ketones
    Reaction(
        id="aldehyde_reduction",
        start="aldehyde", end="alcohol_primary",
        reagents="NaBH4 in aqueous methanol, then dilute acid",
        conditions="Room temperature",
        rxn_type="Nucleophilic addition (reduction)",
        spec="3.3.9",
        note="NaBH4 supplies the hydride ion :H- , which attacks the delta+ "
             "carbon of C=O.",
    ),
    Reaction(
        id="ketone_reduction",
        start="ketone", end="alcohol_secondary",
        reagents="NaBH4 in aqueous methanol, then dilute acid",
        conditions="Room temperature",
        rxn_type="Nucleophilic addition (reduction)",
        spec="3.3.9",
    ),
    Reaction(
        id="aldehyde_oxidation",
        start="aldehyde", end="carboxylic_acid",
        reagents="Acidified potassium dichromate(VI)",
        conditions="Heat under reflux",
        rxn_type="Oxidation",
        spec="3.3.6",
        note="Tollens' reagent (silver mirror) and Fehling's solution (brick-red "
             "precipitate) also oxidise aldehydes - they are the tests that "
             "distinguish an aldehyde from a ketone.",
    ),
    Reaction(
        id="aldehyde_hcn",
        start="aldehyde", end="hydroxynitrile",
        reagents="KCN followed by dilute H2SO4 (a source of HCN)",
        conditions="Room temperature",
        rxn_type="Nucleophilic addition",
        spec="3.3.9", delta_min=1, delta_max=1,
        note="Lengthens the chain by one carbon. The planar C=O is attacked "
             "from either side equally, so the product is a racemic mixture of "
             "the two optical isomers.",
    ),
    Reaction(
        id="ketone_hcn",
        start="ketone", end="hydroxynitrile",
        reagents="KCN followed by dilute H2SO4 (a source of HCN)",
        conditions="Room temperature",
        rxn_type="Nucleophilic addition",
        spec="3.3.9", delta_min=1, delta_max=1,
        note="HCN is highly toxic, so this is done in a fume cupboard - a "
             "common 'suggest a hazard' answer.",
    ),

    # --------------------------------------------------------------- nitriles
    Reaction(
        id="nitrile_hydrolysis",
        start="nitrile", end="carboxylic_acid",
        reagents="Dilute HCl (aq)",
        conditions="Heat under reflux",
        rxn_type="Hydrolysis",
        spec="3.3.9",
        note="Ammonium ions are the other product. The chain keeps the extra "
             "carbon gained from the cyanide.",
    ),
    Reaction(
        id="nitrile_reduction",
        start="nitrile", end="amine_primary",
        reagents="LiAlH4 in dry ether, then dilute acid (or H2 with a Ni catalyst)",
        conditions="Room temperature for LiAlH4; high temperature and pressure "
                   "for H2 / Ni",
        rxn_type="Reduction",
        spec="3.3.11",
        note="LiAlH4 is too expensive for industry, which uses H2 / Ni instead.",
    ),
    Reaction(
        id="hydroxynitrile_hydrolysis",
        start="hydroxynitrile", end="hydroxy_acid",
        reagents="Dilute HCl (aq)",
        conditions="Heat under reflux",
        rxn_type="Hydrolysis",
        spec="3.3.9",
        note="Gives a 2-hydroxycarboxylic acid, which is still a racemic mixture.",
    ),

    # ------------------------------------------------------- carboxylic acids
    Reaction(
        id="acid_esterification",
        start="carboxylic_acid", end="ester",
        reagents="An alcohol, with a few drops of conc. H2SO4 as catalyst",
        conditions="Warm / heat under reflux",
        rxn_type="Esterification (condensation)",
        spec="3.3.10", carbon_rule="ester_from_acid", cost=1.1,
        note="Reversible - remove the ester by distillation to improve the yield.",
    ),
    Reaction(
        id="acid_neutralisation",
        start="carboxylic_acid", end="carboxylate_salt",
        reagents="NaOH(aq), Na2CO3(aq) or NaHCO3(aq)",
        conditions="Room temperature",
        rxn_type="Neutralisation (acid-base)",
        spec="3.3.10",
        note="Fizzing with a carbonate (CO2 released) is the test for a "
             "carboxylic acid.",
    ),
    Reaction(
        id="salt_acidification",
        start="carboxylate_salt", end="carboxylic_acid",
        reagents="Dilute HCl (or dilute H2SO4)",
        conditions="Room temperature",
        rxn_type="Acidification",
        spec="3.3.10",
    ),
    Reaction(
        id="acid_to_acyl_chloride",
        start="carboxylic_acid", end="acyl_chloride",
        reagents="SOCl2 (thionyl chloride)",
        conditions="Room temperature; SO2 and HCl are released",
        rxn_type="Substitution",
        spec="3.3.10", cost=1.4, status=EXTENSION,
        note="AQA ask about the reactions OF acyl chlorides rather than how to "
             "make them - check whether your teacher wants this step.",
    ),
    Reaction(
        id="acid_reduction",
        start="carboxylic_acid", end="alcohol_primary",
        reagents="LiAlH4 in dry ether, then dilute acid",
        conditions="Room temperature",
        rxn_type="Reduction",
        spec="3.3.10", cost=1.2,
        note="LiAlH4 is a much stronger reducing agent than NaBH4, which will "
             "not reduce a carboxylic acid.",
    ),

    # ----------------------------------------------------------------- esters
    Reaction(
        id="ester_acid_hydrolysis_acid",
        start="ester", end="carboxylic_acid",
        reagents="Dilute HCl or dilute H2SO4, with water",
        conditions="Heat under reflux",
        rxn_type="Acid hydrolysis (reversible)",
        spec="3.3.10", carbon_rule="ester_to_acid",
        note="Reversible, so use an excess of water to push it to the right. "
             "An alcohol is the other product.",
    ),
    Reaction(
        id="ester_acid_hydrolysis_alcohol",
        start="ester", end="alcohol_primary",
        reagents="Dilute HCl or dilute H2SO4, with water",
        conditions="Heat under reflux",
        rxn_type="Acid hydrolysis (reversible)",
        spec="3.3.10", carbon_rule="ester_to_alcohol",
        note="The alcohol half of the ester; a carboxylic acid is formed too.",
    ),
    Reaction(
        id="ester_alkaline_hydrolysis",
        start="ester", end="carboxylate_salt",
        reagents="NaOH(aq)",
        conditions="Heat under reflux",
        rxn_type="Alkaline hydrolysis (saponification)",
        spec="3.3.10", carbon_rule="ester_to_salt",
        note="Goes to completion, unlike acid hydrolysis. Add dilute acid "
             "afterwards if you want the carboxylic acid itself.",
    ),
    Reaction(
        id="ester_alkaline_hydrolysis_alcohol",
        start="ester", end="alcohol_primary",
        reagents="NaOH(aq)",
        conditions="Heat under reflux",
        rxn_type="Alkaline hydrolysis (saponification)",
        spec="3.3.10", carbon_rule="ester_to_alcohol",
        note="The alcohol is the other product; the acid is left as its sodium "
             "salt.",
    ),

    # --------------------------------------------------------- acyl chlorides
    Reaction(
        id="acyl_water",
        start="acyl_chloride", end="carboxylic_acid",
        reagents="Water",
        conditions="Room temperature - vigorous, with misty fumes of HCl",
        rxn_type="Addition-elimination",
        spec="3.3.10",
    ),
    Reaction(
        id="acyl_alcohol",
        start="acyl_chloride", end="ester",
        reagents="An alcohol",
        conditions="Room temperature",
        rxn_type="Addition-elimination",
        spec="3.3.10", carbon_rule="ester_from_acyl",
        note="A better way to make an ester than esterification: fast, "
             "high-yielding and irreversible.",
    ),
    Reaction(
        id="acyl_ammonia",
        start="acyl_chloride", end="amide",
        reagents="Excess concentrated ammonia",
        conditions="Room temperature",
        rxn_type="Addition-elimination",
        spec="3.3.10",
        note="Ammonium chloride is also formed, from the HCl released.",
    ),
    Reaction(
        id="acyl_amine",
        start="acyl_chloride", end="n_substituted_amide",
        reagents="A primary amine",
        conditions="Room temperature",
        rxn_type="Addition-elimination",
        spec="3.3.10", delta_min=1, delta_max=None,
    ),

    # -------------------------------------------------------- acid anhydrides
    Reaction(
        id="anhydride_water",
        start="acid_anhydride", end="carboxylic_acid",
        reagents="Water",
        conditions="Warm",
        rxn_type="Addition-elimination",
        spec="3.3.10",
        note="Anhydrides react in the same way as acyl chlorides but less "
             "violently, and a carboxylic acid is the by-product instead of HCl.",
    ),
    Reaction(
        id="anhydride_alcohol",
        start="acid_anhydride", end="ester",
        reagents="An alcohol",
        conditions="Room temperature / warm",
        rxn_type="Addition-elimination",
        spec="3.3.10", carbon_rule="ester_from_acyl",
        note="Cheaper, safer and less corrosive than the acyl chloride - this "
             "is how aspirin is made industrially.",
    ),
    Reaction(
        id="anhydride_ammonia",
        start="acid_anhydride", end="amide",
        reagents="Ammonia",
        conditions="Room temperature",
        rxn_type="Addition-elimination",
        spec="3.3.10",
    ),
    Reaction(
        id="anhydride_amine",
        start="acid_anhydride", end="n_substituted_amide",
        reagents="A primary amine",
        conditions="Room temperature",
        rxn_type="Addition-elimination",
        spec="3.3.10", delta_min=1, delta_max=None,
    ),

    # ----------------------------------------------------------------- amines
    Reaction(
        id="amine_acylation",
        start="amine_primary", end="n_substituted_amide",
        reagents="An acyl chloride (or an acid anhydride)",
        conditions="Room temperature",
        rxn_type="Addition-elimination (acylation)",
        spec="3.3.10", delta_min=1, delta_max=None,
        note="The amine acts as a nucleophile through the lone pair on nitrogen.",
    ),
    Reaction(
        id="amine1_to_amine2",
        start="amine_primary", end="amine_secondary",
        reagents="A halogenoalkane",
        conditions="Heat in a sealed tube; use excess amine",
        rxn_type="Nucleophilic substitution",
        spec="3.3.11", delta_min=1, delta_max=None,
    ),
    Reaction(
        id="amine2_to_amine3",
        start="amine_secondary", end="amine_tertiary",
        reagents="A halogenoalkane",
        conditions="Heat in a sealed tube",
        rxn_type="Nucleophilic substitution",
        spec="3.3.11", delta_min=1, delta_max=None,
    ),
    Reaction(
        id="amine3_to_quaternary",
        start="amine_tertiary", end="quaternary_ammonium",
        reagents="Excess halogenoalkane",
        conditions="Heat in a sealed tube",
        rxn_type="Nucleophilic substitution",
        spec="3.3.11", delta_min=1, delta_max=None,
        note="Quaternary ammonium salts are used as cationic surfactants in "
             "fabric conditioners.",
    ),

    # ------------------------------------------------------ aromatic chemistry
    Reaction(
        id="benzene_nitration",
        start="benzene", end="nitrobenzene",
        reagents="Conc. HNO3 with conc. H2SO4 as catalyst",
        conditions="Warm in a water bath at 50 °C (keep below 55 °C to avoid "
                   "further substitution)",
        rxn_type="Electrophilic substitution (nitration)",
        spec="3.3.7",
        note="The electrophile is NO2+ , generated by H2SO4 protonating HNO3.",
    ),
    Reaction(
        id="nitrobenzene_reduction",
        start="nitrobenzene", end="phenylamine",
        reagents="Tin (Sn) with conc. HCl, then excess NaOH(aq)",
        conditions="Heat under reflux",
        rxn_type="Reduction",
        spec="3.3.11",
        note="The acid gives the phenylammonium salt; the NaOH then liberates "
             "the free phenylamine. Used to make azo dyes.",
    ),
    Reaction(
        id="friedel_crafts_alkylation",
        start="benzene", end="alkylbenzene",
        reagents="A halogenoalkane with anhydrous AlCl3 (the halogen carrier)",
        conditions="Heat under reflux in dry conditions",
        rxn_type="Electrophilic substitution (Friedel-Crafts alkylation)",
        spec="3.3.7", delta_min=1, delta_max=None, cost=1.3,
        note="Further substitution is hard to control because the alkyl group "
             "activates the ring.",
    ),
    Reaction(
        id="friedel_crafts_acylation",
        start="benzene", end="aromatic_ketone",
        reagents="An acyl chloride with anhydrous AlCl3",
        conditions="Heat under reflux in dry ether",
        rxn_type="Electrophilic substitution (Friedel-Crafts acylation)",
        spec="3.3.7", delta_min=2, delta_max=None,
        note="Preferred to alkylation because it stops cleanly after one "
             "substitution.",
    ),
    Reaction(
        id="phenylamine_acylation",
        start="phenylamine", end="n_substituted_amide",
        reagents="An acyl chloride (or an acid anhydride)",
        conditions="Room temperature",
        rxn_type="Addition-elimination (acylation)",
        spec="3.3.10", delta_min=1, delta_max=None,
    ),

    # ------------------------------------------------------------ amino acids
    Reaction(
        id="acid_to_amino_acid",
        start="carboxylic_acid", end="amino_acid",
        reagents="Cl2 with UV light (to make the 2-chloro acid), then excess "
                 "concentrated NH3",
        conditions="Reflux with excess ammonia",
        rxn_type="Free-radical substitution, then nucleophilic substitution",
        spec="3.3.12", cost=2.0, status=EXTENSION,
        note="AQA expect you to recognise and use amino acids rather than "
             "synthesise them - treat this route as background.",
    ),
]


REACTIONS_BY_ID = {r.id: r for r in REACTIONS}


def outgoing(group: str) -> list[Reaction]:
    """Every reaction that starts from ``group``."""
    return [r for r in REACTIONS if r.start == group]


def incoming(group: str) -> list[Reaction]:
    """Every reaction that produces ``group``."""
    return [r for r in REACTIONS if r.end == group]


#: Groups that cannot be turned into anything else on this map.
DEAD_ENDS = {
    "alcohol_tertiary",  # resists oxidation
}
