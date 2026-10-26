#!/usr/bin/env python3
"""Build the self-contained web version from the Python reaction map.

``web/index.html`` is generated, not hand-edited: the reaction data, the
compound names and the list of names the search box accepts all come from
:mod:`chemsynth`, so the browser and the command line can never disagree.

    python3 build_web.py
"""

from __future__ import annotations

import json
from pathlib import Path

from chemsynth.flashcards import all_cards
from chemsynth.parser import (CLASS_ALIASES, NAME_ALIASES, UNSUPPORTED,
                              ParseError, normalise, parse)
from chemsynth.reactions import REACTIONS
from chemsynth.species import GROUPS, STEMS, formula_of, name_of

ROOT = Path(__file__).parent
TEMPLATE = ROOT / "web" / "template.html"
OUTPUT = ROOT / "web" / "index.html"

#: Chain lengths to generate names for.  The search allows two carbons of
#: head-room above the larger molecule, so C12 covers a C10 target.
MAX_CARBONS = 12

EXAMPLES = [
    ("ethene", "propanoic acid"),
    ("propene", "propan-2-ol"),
    ("ethanol", "ethanal"),
    ("1-bromopropane", "propan-1-amine"),
    ("benzene", "phenylamine"),
    ("propan-1-ol", "ethyl propanoate"),
    ("ethanal", "2-hydroxypropanoic acid"),
]


def species_key(group: str, carbons: int | None, acyl: int | None) -> str:
    return f"{group}|{'' if carbons is None else carbons}|{'' if acyl is None else acyl}"


def build_names() -> dict[str, dict[str, str]]:
    """Name and formula for every species the search can reach."""
    names: dict[str, dict[str, str]] = {}

    for group in GROUPS:
        for carbons in range(1, MAX_CARBONS + 1):
            if group == "ester":
                continue
            names[species_key(group, carbons, None)] = {
                "name": name_of(group, carbons),
                "formula": formula_of(group, carbons),
            }

    for total in range(2, MAX_CARBONS + 1):
        names[species_key("ester", total, None)] = {
            "name": name_of("ester", total),
            "formula": "",
        }
        for acyl in range(1, total):
            names[species_key("ester", total, acyl)] = {
                "name": name_of("ester", total, acyl),
                "formula": formula_of("ester", total, acyl),
            }

    return names


def candidate_names() -> set[str]:
    """Everything worth offering the parser, including spelling variants.

    Each candidate is run through the real parser, so the web page accepts
    exactly what the command line accepts.
    """
    candidates: set[str] = set(CLASS_ALIASES) | set(NAME_ALIASES)

    for group in GROUPS:
        for carbons in range(1, MAX_CARBONS + 1):
            candidates.add(name_of(group, carbons))
            formula = formula_of(group, carbons)
            if formula:
                candidates.add(formula)

    for total in range(2, MAX_CARBONS + 1):
        for acyl in range(1, total):
            candidates.add(name_of("ester", total, acyl))

    # Spellings a student is just as likely to type.
    for index, stem in enumerate(STEMS):
        carbons = index + 1
        for halogen in ("chloro", "bromo", "iodo"):
            candidates.add(f"{halogen}{stem}ane")
            for locant in (1, 2):
                candidates.add(f"{locant}-{halogen}{stem}ane")
        candidates.add(f"{stem}ylamine")
        candidates.add(f"1-amino{stem}ane")
        candidates.add(f"{stem}-1-ene")
        candidates.add(f"{stem}-2-ene")
        candidates.add(f"{stem}an-2-one")
        candidates.add(f"{stem}an-2-ol")
        candidates.add(f"{stem}an-1-ol")
        candidates.add(f"sodium {stem}anoate")
        candidates.add(f"poly({stem}ene)")

    return candidates


def build_lookup() -> tuple[dict[str, str], list[str]]:
    """Map typed text to a species key, and collect autocomplete entries."""
    lookup: dict[str, str] = {}
    suggestions: set[str] = set()

    for candidate in candidate_names():
        try:
            species = parse(candidate)
        except ParseError:
            continue
        lookup[normalise(candidate)] = species_key(
            species.group, species.carbons, species.acyl)
        if species.carbons is not None:
            suggestions.add(candidate)

    for text, group in CLASS_ALIASES.items():
        lookup[normalise(text)] = species_key(group, None, None)

    return lookup, sorted(suggestions, key=lambda s: (len(s), s))


def build_cards() -> list[dict]:
    """Flash cards.  Classic cards carry species keys, so any name the parser
    understands counts as a correct answer."""
    cards = []
    for card in all_cards():
        compounds = []
        for i, compound in enumerate(card.compounds):
            entry = {"name": compound.name, "formula": compound.formula,
                     "also": list(compound.also)}
            if card.species:
                s = card.species[i]
                entry["key"] = species_key(s.group, s.carbons, s.acyl)
            compounds.append(entry)
        cards.append({
            "id": card.id,
            "source": None if card.source is None else {
                "paper": card.source.paper, "series": card.source.series,
                "question": card.source.question,
            },
            "compounds": compounds,
            "steps": [
                {"reagents": s.reagents, "conditions": s.conditions,
                 "type": s.rxn_type, "ms": s.from_mark_scheme}
                for s in card.steps
            ],
            "note": card.note,
        })
    return cards


def build_data() -> dict:
    return {
        "groups": {
            group.id: {"label": group.label, "display": group.display,
                       "spec": group.spec}
            for group in GROUPS.values()
        },
        "reactions": [
            {
                "id": r.id, "start": r.start, "end": r.end,
                "reagents": r.reagents, "conditions": r.conditions,
                "type": r.rxn_type, "spec": r.spec, "note": r.note,
                "cost": r.cost, "status": r.status,
                "dmin": r.delta_min, "dmax": r.delta_max,
                "cmin": r.carbons_min, "cmax": r.carbons_max,
                "rule": r.carbon_rule,
            }
            for r in REACTIONS
        ],
        "names": build_names(),
        "lookup": build_lookup()[0],
        "suggestions": build_lookup()[1],
        "unsupported": {normalise(k): v for k, v in UNSUPPORTED.items()},
        "examples": EXAMPLES,
        "cards": build_cards(),
    }


def main() -> None:
    template = TEMPLATE.read_text(encoding="utf-8")
    marker = "/*__DATA__*/{}"
    if marker not in template:
        raise SystemExit(f"{TEMPLATE} no longer contains the data marker")

    payload = json.dumps(build_data(), separators=(",", ":"), ensure_ascii=False)
    OUTPUT.write_text(template.replace(marker, payload), encoding="utf-8")

    data = build_data()
    print(f"wrote {OUTPUT.relative_to(ROOT)}  "
          f"({OUTPUT.stat().st_size / 1024:.0f} kB, "
          f"{len(data['reactions'])} reactions, "
          f"{len(data['cards'])} flash cards, "
          f"{len(data['lookup'])} names understood)")


if __name__ == "__main__":
    main()
