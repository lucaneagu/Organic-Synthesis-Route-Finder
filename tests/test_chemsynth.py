"""Tests for the AQA organic synthesis route finder.

Run with:  python3 -m unittest discover -s tests
"""

import unittest

from chemsynth import find_route, find_routes, parse
from chemsynth.parser import ParseError
from chemsynth.pathfinder import NoRouteFound
from chemsynth.reactions import REACTIONS
from chemsynth.species import GROUPS, Species, formula_of, name_of


def route(start: str, target: str):
    return find_route(parse(start), parse(target))


def reaction_ids(start: str, target: str) -> list[str]:
    return [step.reaction.id for step in route(start, target).steps]


class TestNaming(unittest.TestCase):

    def test_names(self):
        cases = [
            ("alkane", 3, None, "propane"),
            ("alkene", 2, None, "ethene"),
            ("alkene", 4, None, "but-1-ene"),
            ("haloalkane", 1, None, "bromomethane"),
            ("haloalkane", 2, None, "bromoethane"),
            ("haloalkane", 3, None, "1-bromopropane"),
            ("haloalkane_secondary", 3, None, "2-bromopropane"),
            ("haloalkane_secondary", 4, None, "2-bromobutane"),
            ("alcohol_primary", 2, None, "ethanol"),
            ("alcohol_primary", 3, None, "propan-1-ol"),
            ("alcohol_secondary", 3, None, "propan-2-ol"),
            ("alcohol_tertiary", 4, None, "2-methylpropan-2-ol"),
            ("aldehyde", 2, None, "ethanal"),
            ("ketone", 3, None, "propanone"),
            ("carboxylic_acid", 4, None, "butanoic acid"),
            ("nitrile", 3, None, "propanenitrile"),
            ("amine_primary", 3, None, "propan-1-amine"),
            ("amide", 2, None, "ethanamide"),
            ("acyl_chloride", 2, None, "ethanoyl chloride"),
            ("hydroxynitrile", 3, None, "2-hydroxypropanenitrile"),
            ("ester", 4, 2, "ethyl ethanoate"),
            ("ester", 5, 3, "ethyl propanoate"),
            ("ester", 5, 2, "propyl ethanoate"),
            ("alkylbenzene", 7, None, "methylbenzene"),
            ("aromatic_ketone", 8, None, "phenylethanone"),
        ]
        for group, carbons, acyl, expected in cases:
            with self.subTest(expected=expected):
                self.assertEqual(name_of(group, carbons, acyl), expected)

    def test_formulae(self):
        cases = [
            ("alkane", 1, None, "CH4"),
            ("alkene", 2, None, "CH2=CH2"),
            ("alkene", 3, None, "CH3CH=CH2"),
            ("haloalkane", 1, None, "CH3Br"),
            ("haloalkane_secondary", 3, None, "CH3CHBrCH3"),
            ("haloalkane_secondary", 4, None, "CH3CH2CHBrCH3"),
            ("alcohol_primary", 1, None, "CH3OH"),
            ("alcohol_primary", 2, None, "CH3CH2OH"),
            ("alcohol_secondary", 3, None, "CH3CH(OH)CH3"),
            ("aldehyde", 1, None, "HCHO"),
            ("aldehyde", 2, None, "CH3CHO"),
            ("carboxylic_acid", 1, None, "HCOOH"),
            ("carboxylic_acid", 2, None, "CH3COOH"),
            ("nitrile", 2, None, "CH3CN"),
            ("ester", 4, 2, "CH3COOCH2CH3"),
            ("acid_anhydride", 2, None, "CH3COOCOCH3"),
        ]
        for group, carbons, acyl, expected in cases:
            with self.subTest(expected=expected):
                self.assertEqual(formula_of(group, carbons, acyl), expected)

    def test_every_group_names_without_crashing(self):
        for group in GROUPS:
            for carbons in range(1, 11):
                self.assertTrue(name_of(group, carbons))
                formula_of(group, carbons)  # must not raise


class TestParser(unittest.TestCase):

    def test_systematic_names(self):
        cases = {
            "ethanol": Species("alcohol_primary", 2),
            "propan-1-ol": Species("alcohol_primary", 3),
            "propan-2-ol": Species("alcohol_secondary", 3),
            "1-bromopropane": Species("haloalkane", 3),
            "bromoethane": Species("haloalkane", 2),
            "2-bromopropane": Species("haloalkane_secondary", 3),
            "2-chlorobutane": Species("haloalkane_secondary", 4),
            "3-bromopentane": Species("haloalkane_secondary", 5),
            "2-bromoethane": Species("haloalkane", 2),
            "butanoic acid": Species("carboxylic_acid", 4),
            "propanenitrile": Species("nitrile", 3),
            "ethanal": Species("aldehyde", 2),
            "propanone": Species("ketone", 3),
            "but-1-ene": Species("alkene", 4),
            "benzene": Species("benzene", 6),
            "nitrobenzene": Species("nitrobenzene", 6),
            "phenylamine": Species("phenylamine", 6),
            "methylbenzene": Species("alkylbenzene", 7),
            "ethyl propanoate": Species("ester", 5, 3),
            "sodium ethanoate": Species("carboxylate_salt", 2),
            "ethanoic anhydride": Species("acid_anhydride", 2),
        }
        for text, expected in cases.items():
            with self.subTest(text=text):
                self.assertEqual(parse(text), expected)

    def test_common_names_and_case(self):
        self.assertEqual(parse("Acetone"), Species("ketone", 3))
        self.assertEqual(parse("  ACETIC ACID "), Species("carboxylic_acid", 2))
        self.assertEqual(parse("aniline"), Species("phenylamine", 6))
        self.assertEqual(parse("toluene"), Species("alkylbenzene", 7))

    def test_classes_of_compound(self):
        self.assertEqual(parse("a primary alcohol"), Species("alcohol_primary", None))
        self.assertEqual(parse("carboxylic acid"), Species("carboxylic_acid", None))
        self.assertEqual(parse("aldehydes"), Species("aldehyde", None))

    def test_formulae(self):
        self.assertEqual(parse("CH3CH2OH"), Species("alcohol_primary", 2))
        self.assertEqual(parse("CH3COOH"), Species("carboxylic_acid", 2))
        self.assertEqual(parse("CH3CHBrCH3"), Species("haloalkane_secondary", 3))

    def test_round_trip_through_names(self):
        for group in ("alkane", "alkene", "haloalkane", "haloalkane_secondary",
                      "alcohol_primary",
                      "alcohol_secondary", "aldehyde", "ketone",
                      "carboxylic_acid", "nitrile", "amine_primary", "amide",
                      "acyl_chloride", "hydroxynitrile"):
            for carbons in range(3, 7):
                name = name_of(group, carbons)
                with self.subTest(name=name):
                    self.assertEqual(parse(name), Species(group, carbons))

    def test_unknown_name_suggests_alternatives(self):
        with self.assertRaises(ParseError) as caught:
            parse("ethanoll")
        self.assertIn("ethanol", caught.exception.suggestions)

    def test_explicitly_unsupported(self):
        with self.assertRaises(ParseError) as caught:
            parse("benzoic acid")
        self.assertIn("not on the AQA", str(caught.exception))


class TestRoutes(unittest.TestCase):

    def test_single_step_oxidation(self):
        self.assertEqual(reaction_ids("ethanol", "ethanoic acid"),
                         ["alcohol1_to_acid"])
        self.assertEqual(reaction_ids("propan-2-ol", "propanone"),
                         ["alcohol2_to_ketone"])

    def test_distil_for_the_aldehyde(self):
        steps = route("ethanol", "ethanal").steps
        self.assertEqual(len(steps), 1)
        self.assertIn("DISTIL", steps[0].reaction.conditions)

    def test_chain_extension_uses_cyanide(self):
        self.assertEqual(reaction_ids("1-bromopropane", "butanoic acid"),
                         ["haloalkane_cyanide", "nitrile_hydrolysis"])

    def test_chain_extension_via_hydroxynitrile(self):
        self.assertEqual(reaction_ids("ethanal", "2-hydroxypropanoic acid"),
                         ["aldehyde_hcn", "hydroxynitrile_hydrolysis"])

    def test_nitrile_to_amine_keeps_chain_length(self):
        found = route("bromoethane", "propan-1-amine")
        self.assertEqual([s.reaction.id for s in found.steps],
                         ["haloalkane_cyanide", "nitrile_reduction"])
        self.assertEqual(found.end.carbons, 3)

    def test_readme_example(self):
        self.assertEqual(reaction_ids("ethene", "propanoic acid"),
                         ["alkene_hydrogen_halide", "haloalkane_cyanide",
                          "nitrile_hydrolysis"])

    def test_markovnikov_gives_the_2_halogenoalkane(self):
        self.assertEqual(reaction_ids("propene", "2-bromopropane"),
                         ["alkene_hydrogen_halide_markovnikov"])
        self.assertEqual(reaction_ids("ethene", "bromoethane"),
                         ["alkene_hydrogen_halide"])

    def test_hbr_never_puts_the_halogen_on_carbon_1(self):
        """Propene + HBr gives 2-bromopropane, so no route may use it for 1-."""
        for start, target in [("propene", "1-bromopropane"),
                              ("propene", "propan-1-ol"),
                              ("propene", "butanoic acid"),
                              ("but-1-ene", "pentanenitrile")]:
            with self.subTest(start=start, target=target):
                self.assertNotIn("alkene_hydrogen_halide",
                                 reaction_ids(start, target))

    def test_secondary_halogenoalkane_keeps_its_position(self):
        self.assertEqual(reaction_ids("2-bromopropane", "propan-2-ol"),
                         ["haloalkane2_hydrolysis"])
        self.assertEqual(reaction_ids("propan-2-ol", "2-bromopropane"),
                         ["alcohol2_to_haloalkane"])

    def test_cannot_move_a_group_from_carbon_2_to_carbon_1(self):
        """The only way across is the messy free-radical step via the alkane."""
        ids = reaction_ids("propan-2-ol", "propan-1-ol")
        self.assertIn("alkane_halogenation", ids)
        self.assertNotIn("haloalkane_hydrolysis", ids[:ids.index("alkane_halogenation")])

    def test_aromatic_route(self):
        self.assertEqual(reaction_ids("benzene", "phenylamine"),
                         ["benzene_nitration", "nitrobenzene_reduction"])

    def test_ester_keeps_the_two_halves_straight(self):
        """propan-1-ol supplies the acid half of ethyl propanoate."""
        found = route("propan-1-ol", "ethyl propanoate")
        self.assertEqual([s.reaction.id for s in found.steps],
                         ["alcohol1_to_acid", "acid_esterification"])
        self.assertEqual(found.end.name, "ethyl propanoate")

    def test_ester_from_the_alcohol_half(self):
        found = route("propan-1-ol", "propyl ethanoate")
        self.assertEqual(found.length, 1)
        self.assertEqual(found.end.name, "propyl ethanoate")

    def test_ester_hydrolysis_returns_both_halves(self):
        self.assertEqual(route("ethyl propanoate", "propanoic acid").length, 1)
        self.assertEqual(route("ethyl propanoate", "ethanol").length, 1)

    def test_co_reactants_cannot_smuggle_in_a_carbon_skeleton(self):
        """You cannot 'make' decanoic acid by hydrolysing someone else's ester.

        The only AQA ways to lengthen a chain add one carbon at a time, so a
        C1 to C10 synthesis has to be long.
        """
        found = route("methanol", "decanoic acid")
        self.assertGreater(found.length, 10)

    def test_start_equals_target(self):
        self.assertEqual(route("ethanol", "ethanol").length, 0)

    def test_generic_target_keeps_the_starting_chain(self):
        found = route("ethanol", "a carboxylic acid")
        self.assertEqual(found.end.name, "ethanoic acid")

    def test_no_route_from_a_dead_end(self):
        with self.assertRaises(NoRouteFound):
            route("an addition polymer", "ethanol")

    def test_tertiary_alcohols_are_not_oxidised(self):
        with self.assertRaises(NoRouteFound):
            route("2-methylpropan-2-ol", "propanone")

    def test_alternatives_are_distinct_and_ordered(self):
        routes = find_routes(parse("ethanol"), parse("ethanoic acid"), limit=3)
        self.assertGreater(len(routes), 1)
        lengths = [r.length for r in routes]
        self.assertEqual(lengths, sorted(lengths))
        signatures = {r.signature for r in routes}
        self.assertEqual(len(signatures), len(routes))

    def test_every_step_is_chemically_connected(self):
        """The product of each step must be the starting material of the next."""
        found = route("propene", "butanoic acid")
        for earlier, later in zip(found.steps, found.steps[1:]):
            self.assertEqual(earlier.end, later.start)
            self.assertEqual(earlier.reaction.end, later.reaction.start)


class TestReactionData(unittest.TestCase):

    def test_groups_exist(self):
        for reaction in REACTIONS:
            with self.subTest(reaction=reaction.id):
                self.assertIn(reaction.start, GROUPS)
                self.assertIn(reaction.end, GROUPS)

    def test_ids_are_unique(self):
        ids = [r.id for r in REACTIONS]
        self.assertEqual(len(ids), len(set(ids)))

    def test_every_reaction_is_fully_described(self):
        for reaction in REACTIONS:
            with self.subTest(reaction=reaction.id):
                self.assertTrue(reaction.reagents)
                self.assertTrue(reaction.conditions)
                self.assertTrue(reaction.rxn_type)
                self.assertRegex(reaction.spec, r"^3\.3\.\d+$")

    def test_secondary_halogenoalkanes_make_no_branched_products(self):
        """KCN or NH3 on a 2-halogenoalkane gives a branched chain, which the
        straight-chain names used here cannot describe."""
        for reaction in REACTIONS:
            if reaction.start == "haloalkane_secondary":
                with self.subTest(reaction=reaction.id):
                    self.assertIn(reaction.end, ("alcohol_secondary", "alkene"))

    def test_carbon_ranges_are_sensible(self):
        for reaction in REACTIONS:
            with self.subTest(reaction=reaction.id):
                if reaction.carbon_rule:
                    continue
                if reaction.delta_min is not None and reaction.delta_max is not None:
                    self.assertLessEqual(reaction.delta_min, reaction.delta_max)


class TestFlashcards(unittest.TestCase):

    def setUp(self):
        from chemsynth.flashcards import PAST_PAPER_CARDS, all_cards
        self.cards = all_cards()
        self.past = PAST_PAPER_CARDS

    def test_ids_are_unique(self):
        ids = [card.id for card in self.cards]
        self.assertEqual(len(ids), len(set(ids)))

    def test_every_card_is_a_complete_route(self):
        for card in self.cards:
            with self.subTest(card=card.id):
                self.assertGreaterEqual(len(card.steps), 1)
                self.assertEqual(len(card.compounds), len(card.steps) + 1)
                for compound in card.compounds:
                    self.assertTrue(compound.name)
                for step in card.steps:
                    self.assertTrue(step.reagents)
                    self.assertTrue(step.conditions)
                    self.assertTrue(step.rxn_type)

    def test_past_paper_sources_are_real_papers(self):
        for card in self.past:
            with self.subTest(card=card.id):
                self.assertIn(card.source.paper, ("7405/2", "7405/3"))
                self.assertRegex(card.source.series, r"^June 20(19|2[0-4])$")
                self.assertRegex(card.source.question, r"^\d\d$")

    def test_every_year_from_2019_to_2024_is_covered(self):
        years = {card.source.series[-4:] for card in self.past}
        self.assertEqual(years, {str(y) for y in range(2019, 2025)})

    def test_classic_cards_come_from_the_reaction_map(self):
        classics = [card for card in self.cards if card.source is None]
        self.assertGreater(len(classics), 5)
        for card in classics:
            with self.subTest(card=card.id):
                self.assertEqual(len(card.species), len(card.compounds))

    def test_markovnikov_cards_agree_with_the_route_finder(self):
        """The 2024 paper's but-1-ene + HBr step must match the map."""
        self.assertEqual(reaction_ids("but-1-ene", "2-bromobutane"),
                         ["alkene_hydrogen_halide_markovnikov"])
        self.assertEqual(reaction_ids("2-bromobutane", "butan-2-ol"),
                         ["haloalkane2_hydrolysis"])

    def test_web_build_includes_every_card(self):
        import build_web
        cards = build_web.build_cards()
        self.assertEqual(len(cards), len(self.cards))
        for card in cards:
            if card["source"] is None:
                for compound in card["compounds"]:
                    self.assertIn("key", compound)


if __name__ == "__main__":
    unittest.main()
