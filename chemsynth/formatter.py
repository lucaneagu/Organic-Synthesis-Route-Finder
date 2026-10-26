"""Lay a route out neatly for the terminal."""

from __future__ import annotations

import os
import sys
import textwrap

from .pathfinder import Route
from .species import GROUPS, Species

ARROW = "->"
LABEL_WIDTH = 12


class Style:
    """ANSI styling that quietly turns itself off when it is not wanted."""

    def __init__(self, enabled: bool = True):
        self.enabled = enabled

    def _wrap(self, code: str, text: str) -> str:
        return f"\033[{code}m{text}\033[0m" if self.enabled else text

    def bold(self, text: str) -> str:
        return self._wrap("1", text)

    def dim(self, text: str) -> str:
        return self._wrap("2", text)

    def cyan(self, text: str) -> str:
        return self._wrap("36", text)

    def green(self, text: str) -> str:
        return self._wrap("32", text)

    def yellow(self, text: str) -> str:
        return self._wrap("33", text)

    def red(self, text: str) -> str:
        return self._wrap("31", text)


def colour_enabled(force: bool | None = None) -> bool:
    if force is not None:
        return force
    if os.environ.get("NO_COLOR"):
        return False
    return sys.stdout.isatty()


def _field(style: Style, label: str, value: str, width: int, indent: int = 4) -> list[str]:
    """A wrapped 'Label   value' block."""
    pad = " " * indent
    body_width = max(20, width - indent - LABEL_WIDTH)
    wrapped = textwrap.wrap(value, body_width) or [""]
    head = f"{pad}{style.dim(label.ljust(LABEL_WIDTH))}{wrapped[0]}"
    tail = [f"{pad}{' ' * LABEL_WIDTH}{line}" for line in wrapped[1:]]
    return [head, *tail]


def format_species(species: Species, style: Style) -> str:
    name = species.name
    formula = species.formula
    return f"{name}  ({formula})" if formula else name


def format_route(route: Route, style: Style, width: int = 78,
                 title: str | None = None) -> str:
    lines: list[str] = []

    if title:
        lines.append(style.bold(title))
        lines.append("")

    if not route.steps:
        lines.append("  Start and target are the same compound - no reaction needed.")
        return "\n".join(lines)

    for number, step in enumerate(route.steps, start=1):
        reaction = step.reaction
        heading = (
            f"  {style.bold(f'STEP {number}')}  "
            f"{style.cyan(step.start.name)} {ARROW} {style.green(step.end.name)}"
        )
        lines.append(heading)

        start_formula, end_formula = step.start.formula, step.end.formula
        if start_formula and end_formula:
            lines.append(
                f"    {' ' * LABEL_WIDTH}"
                f"{style.dim(f'{start_formula} {ARROW} {end_formula}')}"
            )

        lines.extend(_field(style, "Reagents", reaction.reagents, width))
        lines.extend(_field(style, "Conditions", reaction.conditions, width))
        lines.extend(_field(style, "Type", reaction.rxn_type, width))

        delta = step.carbon_delta
        if delta:
            plural = "" if abs(delta) == 1 else "s"
            if reaction.carbon_rule:
                change = (f"the other reactant contributes {abs(delta)} "
                          f"carbon{plural}")
            else:
                verb = "lengthened" if delta > 0 else "shortened"
                change = f"chain {verb} by {abs(delta)} carbon{plural}"
            lines.extend(_field(style, "Carbons", style.yellow(change), width))

        if reaction.note:
            lines.extend(_field(style, "Note", reaction.note, width))

        spec = f"AQA {reaction.spec}"
        if reaction.status != "core":
            spec += "  (beyond the core specification - check with your teacher)"
        lines.extend(_field(style, "Spec", style.dim(spec), width))
        lines.append("")

    return "\n".join(lines).rstrip()


def format_header(start: Species, target: Species, route: Route, style: Style,
                  width: int = 78) -> str:
    rule = "=" * width
    steps = route.length
    step_word = "step" if steps == 1 else "steps"

    target_name = target.name
    if target.carbons is None and route.steps:
        target_name = route.end.name
    summary = f"{start.name} {ARROW} {target_name}"
    detail = f"{steps} {step_word}"
    if start.carbons is not None and target.carbons is not None:
        if start.carbons == target.carbons:
            detail += f" | carbon chain stays at {start.carbons}"
        else:
            detail += f" | carbon chain {start.carbons} {ARROW} {target.carbons}"

    return "\n".join([
        style.cyan(rule),
        f" {style.bold(summary.upper())}",
        f" {style.dim(detail)}",
        style.cyan(rule),
    ])


def format_summary(route: Route, style: Style, width: int = 78) -> str:
    """One-line reminder of the whole sequence."""
    if not route.steps:
        return ""
    names = [route.steps[0].start.name] + [s.end.name for s in route.steps]
    chain = f"  {ARROW}  ".join(names)
    lines = [style.bold("  OVERALL")]
    for line in textwrap.wrap(chain, max(20, width - 4)):
        lines.append(f"    {line}")
    if route.uses_extension():
        lines.append("")
        lines.append(style.yellow(
            "  ! This route uses a step beyond the core AQA specification."
        ))
    return "\n".join(lines)


def format_no_route(start: Species, target: Species, reachable: set[str],
                    style: Style, width: int = 78) -> str:
    lines = [
        style.red(f"  No AQA route found from {start.name} to {target.name}."),
        "",
    ]
    if target.group not in reachable:
        lines.append(
            f"  Nothing on the map turns {GROUPS[start.group].label}s into "
            f"{GROUPS[target.group].label}s."
        )
    else:
        lines.append(
            "  The functional groups can be linked, but not with these carbon "
            "chain lengths - AQA reactions change the chain by at most one "
            "carbon at a time (via KCN or HCN)."
        )
    lines.append("")
    labels = sorted(GROUPS[g].label for g in reachable if g != start.group)
    if labels:
        wrapped = textwrap.wrap(", ".join(labels), max(20, width - 4))
        lines.append(f"  From {start.name} you can reach:")
        lines.extend(f"    {line}" for line in wrapped)
    return "\n".join(lines)


def route_to_dict(route: Route) -> dict:
    """JSON-friendly view of a route."""
    return {
        "steps": len(route.steps),
        "start": {
            "name": route.start.name,
            "formula": route.start.formula,
            "group": route.start.group,
            "carbons": route.start.carbons,
        } if route.steps else None,
        "target": {
            "name": route.end.name,
            "formula": route.end.formula,
            "group": route.end.group,
            "carbons": route.end.carbons,
        } if route.steps else None,
        "uses_extension": route.uses_extension(),
        "sequence": [
            {
                "step": i,
                "from": {"name": s.start.name, "formula": s.start.formula},
                "to": {"name": s.end.name, "formula": s.end.formula},
                "reagents": s.reaction.reagents,
                "conditions": s.reaction.conditions,
                "type": s.reaction.rxn_type,
                "carbon_change": s.carbon_delta,
                "note": s.reaction.note,
                "spec": s.reaction.spec,
                "status": s.reaction.status,
            }
            for i, s in enumerate(route.steps, start=1)
        ],
    }
