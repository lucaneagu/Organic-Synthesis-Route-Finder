"""Command line interface: chemsynth <start chemical> <target chemical>."""

from __future__ import annotations

import argparse
import json
import sys

from . import formatter
from .formatter import Style
from .parser import CLASS_ALIASES, ParseError, parse
from .pathfinder import NoRouteFound, find_routes, reachable_groups
from .reactions import REACTIONS
from .species import GROUPS, name_of

EXAMPLES = """\
examples:
  chemsynth ethene "ethanoic acid"
  chemsynth 1-bromopropane "butanoic acid"
  chemsynth propan-2-ol propanone --routes 3
  chemsynth benzene phenylamine
  chemsynth ethanol "a carboxylic acid"
  chemsynth --list compounds
"""


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="chemsynth",
        description="Work out the reagents and conditions needed to get from "
                    "one organic compound to another, using the reactions on "
                    "the AQA A-level Chemistry specification.",
        epilog=EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("start", nargs="?", help="the chemical you are starting from")
    p.add_argument("target", nargs="?", help="the chemical you want to make")
    p.add_argument("-n", "--routes", type=int, default=1, metavar="N",
                   help="show up to N routes, shortest first (default 1)")
    p.add_argument("--json", action="store_true",
                   help="print the route as JSON instead of formatted text")
    p.add_argument("--list", dest="listing",
                   choices=["compounds", "groups", "reactions"],
                   help="list what the program knows about, then exit")
    p.add_argument("--width", type=int, default=78, metavar="COLS",
                   help="output width in characters (default 78)")
    p.add_argument("--no-color", dest="color", action="store_false", default=None,
                   help="turn off colour output")
    return p


# --------------------------------------------------------------------------
# Listings
# --------------------------------------------------------------------------


def show_groups(style: Style) -> None:
    print(style.bold("Functional group classes on the map:"))
    print()
    for group in GROUPS.values():
        print(f"  {group.label:<28}{style.dim(group.display):<34}"
              f"{style.dim('AQA ' + group.spec)}")


def show_compounds(style: Style) -> None:
    print(style.bold("Names you can type (any chain length from C1 to C10):"))
    print()
    interesting = [
        "alkane", "alkene", "haloalkane", "alcohol_primary", "alcohol_secondary",
        "aldehyde", "ketone", "carboxylic_acid", "ester", "acyl_chloride",
        "amide", "nitrile", "hydroxynitrile", "amine_primary", "benzene",
        "nitrobenzene", "phenylamine", "alkylbenzene",
    ]
    for group in interesting:
        samples = []
        for carbons in (2, 3, 7):
            name = name_of(group, carbons)
            if not name.startswith("a ") and name not in samples:
                samples.append(name)
        print(f"  {GROUPS[group].label:<24}e.g. {style.cyan(', '.join(samples[:2]))}")
    print()
    print(style.bold("You can also name a whole class of compound, such as:"))
    print()
    classes = ["a primary alcohol", "an aldehyde", "a carboxylic acid",
               "a nitrile", "an ester", "a primary amine"]
    print("  " + ", ".join(style.cyan(c) for c in classes))
    print()
    print(style.dim("  Condensed formulae work too, e.g. CH3CH2OH or CH3COOH."))


def show_reactions(style: Style) -> None:
    print(style.bold(f"All {len(REACTIONS)} reactions on the map:"))
    print()
    for reaction in sorted(REACTIONS, key=lambda r: (r.spec, r.id)):
        start = GROUPS[reaction.start].label
        end = GROUPS[reaction.end].label
        print(f"  {style.cyan(start)} {formatter.ARROW} {style.green(end)}")
        print(f"      {reaction.reagents}; {reaction.conditions}")
        print(f"      {style.dim(reaction.rxn_type + '  |  AQA ' + reaction.spec)}")
        print()


# --------------------------------------------------------------------------
# Running a search
# --------------------------------------------------------------------------


def read_species(prompt: str, style: Style):
    """Ask until the user types something the parser understands."""
    while True:
        try:
            text = input(style.bold(prompt))
        except (EOFError, KeyboardInterrupt):
            print()
            raise SystemExit(0)
        if not text.strip():
            continue
        try:
            return parse(text)
        except ParseError as err:
            print(style.red(f"  {err}"))
            if err.suggestions:
                print(style.dim(f"  Did you mean: {', '.join(err.suggestions)}?"))
            print()


def run_search(start, target, args, style: Style) -> int:
    try:
        routes = find_routes(start, target, limit=max(1, args.routes))
    except NoRouteFound:
        if args.json:
            print(json.dumps({"error": "no route found",
                              "start": start.name, "target": target.name}, indent=2))
        else:
            print(formatter.format_no_route(
                start, target, reachable_groups(start.group), style, args.width))
        return 1

    if args.json:
        payload = {
            "start": start.name,
            "target": target.name,
            "routes": [formatter.route_to_dict(r) for r in routes],
        }
        print(json.dumps(payload, indent=2))
        return 0

    print()
    print(formatter.format_header(start, target, routes[0], style, args.width))
    print()

    for index, route in enumerate(routes, start=1):
        title = None
        if len(routes) > 1:
            label = "BEST ROUTE" if index == 1 else f"ALTERNATIVE {index - 1}"
            plural = "" if route.length == 1 else "s"
            title = f"  {label}  ({route.length} step{plural})"
        print(formatter.format_route(route, style, args.width, title))
        print()

    print(formatter.format_summary(routes[0], style, args.width))
    print()
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    style = Style(formatter.colour_enabled(args.color))

    if args.listing:
        {"groups": show_groups, "compounds": show_compounds,
         "reactions": show_reactions}[args.listing](style)
        return 0

    if args.start and args.target:
        try:
            start = parse(args.start)
            target = parse(args.target)
        except ParseError as err:
            print(style.red(str(err)), file=sys.stderr)
            if err.suggestions:
                print(f"Did you mean: {', '.join(err.suggestions)}?", file=sys.stderr)
            else:
                print("Try 'chemsynth --list compounds' to see what is supported.",
                      file=sys.stderr)
            return 2
        return run_search(start, target, args, style)

    if args.start and not args.target:
        print("Give both a starting chemical and a target chemical.", file=sys.stderr)
        print('  e.g. chemsynth ethene "ethanoic acid"', file=sys.stderr)
        return 2

    # No arguments: interactive mode.
    print()
    print(style.bold("AQA A-level organic synthesis route finder"))
    print(style.dim("Type the chemical you have, then the one you want. "
                    "Ctrl-C to quit."))
    print()
    while True:
        start = read_species("  Start chemical:  ", style)
        target = read_species("  Target chemical: ", style)
        run_search(start, target, args, style)
        print(style.dim("  " + "-" * 40))
        print()


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
