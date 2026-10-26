"""Search the reaction map for routes from one compound to another.

The search space is a graph whose nodes are ``(functional group, number of
carbons)`` states and whose edges are the reactions in :mod:`reactions`.  A
reaction that can add carbons (esterification, Friedel-Crafts, ...) becomes
several edges, one per possible chain length.
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass, field

from .reactions import REACTIONS, Reaction
from .species import Species

#: Head-room, in carbon atoms, above the larger of the start and target
#: molecules.  Two is enough for the "lengthen then shorten" routes that come
#: up at A-level while keeping the search small.
CARBON_HEADROOM = 2

#: Safety valve for the alternative-route search.
NODE_BUDGET = 300_000


#: A node in the search graph: functional group, total carbons, the ester's
#: acid-derived carbons, and which half of an ester descends from the starting
#: material ("acyl", "alkyl", or None when both do).
State = tuple[str, int | None, int | None, str | None]


@dataclass(frozen=True)
class Step:
    """One arrow of a synthesis."""

    reaction: Reaction
    start: Species
    end: Species

    @property
    def carbon_delta(self) -> int | None:
        if self.start.carbons is None or self.end.carbons is None:
            return None
        return self.end.carbons - self.start.carbons


@dataclass
class Route:
    steps: list[Step] = field(default_factory=list)

    @property
    def cost(self) -> float:
        return sum(s.reaction.cost for s in self.steps)

    @property
    def length(self) -> int:
        return len(self.steps)

    @property
    def start(self) -> Species:
        return self.steps[0].start

    @property
    def end(self) -> Species:
        return self.steps[-1].end

    @property
    def signature(self) -> tuple:
        return tuple((s.reaction.id, s.end.carbons, s.end.acyl) for s in self.steps)

    def uses_extension(self) -> bool:
        return any(s.reaction.status != "core" for s in self.steps)


class NoRouteFound(Exception):
    """Raised when the map contains no path between the two compounds."""


# --------------------------------------------------------------------------
# Graph expansion
# --------------------------------------------------------------------------

_BY_START: dict[str, list[Reaction]] = {}
for _r in REACTIONS:
    _BY_START.setdefault(_r.start, []).append(_r)


def _ester_successors(reaction: Reaction, state: State,
                      max_carbons: int) -> list[State]:
    """States produced by a reaction that makes or breaks an ester.

    An ester has two carbon skeletons and only one of them can have come from
    the compound the student started with.  ``lineage`` records which, and
    hydrolysis is then only allowed to hand back the half they actually own -
    otherwise the search would happily "make" any molecule at all by picking a
    convenient co-reactant and hydrolysing it straight back off.
    """
    _, carbons, acyl, lineage = state
    if carbons is None:
        return [(reaction.end, None, None, None)]

    rule = reaction.carbon_rule

    if rule == "ester_from_alcohol":
        # The starting alcohol supplies the alkyl half; the acid half is the
        # co-reactant, so its size is ours to choose.
        return [("ester", carbons + acid, acid, "alkyl")
                for acid in range(1, max_carbons - carbons + 1)]

    if rule in ("ester_from_acid", "ester_from_acyl"):
        # Here the starting material supplies the acyl half instead.
        return [("ester", carbons + alcohol, carbons, "acyl")
                for alcohol in range(1, max_carbons - carbons + 1)]

    # Breaking an ester apart: work out which halves we are entitled to keep.
    if acyl is None:
        splits = [(a, carbons - a) for a in range(1, carbons)]
    else:
        splits = [(acyl, carbons - acyl)]

    if rule == "ester_to_alcohol":
        if lineage == "acyl":
            return []
        return [(reaction.end, alkyl, None, None) for _, alkyl in splits if alkyl >= 1]

    # ester_to_acid / ester_to_salt
    if lineage == "alkyl":
        return []
    return [(reaction.end, acid, None, None) for acid, _ in splits if acid >= 1]


def _successors(state: State, max_carbons: int) -> list[tuple[State, Reaction]]:
    """Every state reachable from ``state`` in one reaction."""
    group, carbons, _, _ = state
    out: list[tuple[State, Reaction]] = []

    for reaction in _BY_START.get(group, []):
        if not reaction.allows(carbons):
            continue

        if reaction.carbon_rule:
            for next_state in _ester_successors(reaction, state, max_carbons):
                if next_state[1] is None or 1 <= next_state[1] <= max_carbons:
                    out.append((next_state, reaction))
            continue

        if carbons is None:
            out.append(((reaction.end, None, None, None), reaction))
            continue

        # An open-ended range is clamped to what a real molecule allows: never
        # fewer than one carbon, never more than the head-room limit.
        low = reaction.delta_min if reaction.delta_min is not None else 1 - carbons
        high = reaction.delta_max if reaction.delta_max is not None else max_carbons - carbons

        for delta in range(low, high + 1):
            new_carbons = carbons + delta
            if 1 <= new_carbons <= max_carbons:
                out.append(((reaction.end, new_carbons, None, None), reaction))

    return out


def carbon_limit(start: Species, target: Species) -> int:
    counts = [c for c in (start.carbons, target.carbons) if c is not None]
    return (max(counts) if counts else 6) + CARBON_HEADROOM


def _route_from_trail(
    trail: list[tuple[State, Reaction]], start_state: State
) -> Route:
    steps: list[Step] = []
    current = start_state
    for state, reaction in trail:
        steps.append(
            Step(reaction,
                 Species(current[0], current[1], current[2]),
                 Species(state[0], state[1], state[2]))
        )
        current = state
    return Route(steps)


# --------------------------------------------------------------------------
# Shortest route (Dijkstra on reaction cost)
# --------------------------------------------------------------------------


def _setup(start: Species, target: Species):
    """Start state and goal test for a search.

    Carbon counts are followed whenever the starting compound has one, so the
    steps can be named properly, but they are only *checked* when the target is
    a specific molecule rather than a whole class of compound.
    """
    start_state: State = (start.group, start.carbons, start.acyl, None)
    check_carbons = start.carbons is not None and target.carbons is not None

    def is_goal(state: State) -> bool:
        if state[0] != target.group:
            return False
        if not check_carbons:
            return True
        if state[1] != target.carbons:
            return False
        if target.acyl is not None and state[2] != target.acyl:
            return False
        return True

    return start_state, is_goal


def find_route(start: Species, target: Species) -> Route:
    """Cheapest route from ``start`` to ``target``.

    "Cheapest" is the number of steps, nudged by a small penalty on messy or
    off-specification reactions so the textbook route comes out on top.
    """
    if (start.group == target.group and start.carbons == target.carbons
            and (target.acyl is None or start.acyl == target.acyl)):
        return Route([])

    max_carbons = carbon_limit(start, target)
    start_state, is_goal = _setup(start, target)

    queue: list[tuple[float, int, State]] = [(0.0, 0, start_state)]
    best: dict[State, float] = {start_state: 0.0}
    came_from: dict[State, tuple[State, Reaction]] = {}
    counter = 1

    while queue:
        cost, _, state = heapq.heappop(queue)
        if cost > best.get(state, float("inf")):
            continue
        if is_goal(state) and state != start_state:
            trail: list[tuple[State, Reaction]] = []
            node = state
            while node != start_state:
                previous, reaction = came_from[node]
                trail.append((node, reaction))
                node = previous
            trail.reverse()
            return _route_from_trail(trail, start_state)

        for next_state, reaction in _successors(state, max_carbons):
            new_cost = cost + reaction.cost
            if new_cost < best.get(next_state, float("inf")):
                best[next_state] = new_cost
                came_from[next_state] = (state, reaction)
                heapq.heappush(queue, (new_cost, counter, next_state))
                counter += 1

    raise NoRouteFound(
        f"No route from {start.name} to {target.name} using AQA A-level reactions."
    )


# --------------------------------------------------------------------------
# Alternative routes
# --------------------------------------------------------------------------


def find_routes(start: Species, target: Species, limit: int = 3,
                extra_steps: int = 1) -> list[Route]:
    """Up to ``limit`` routes, shortest first.

    Alternatives may be up to ``extra_steps`` reactions longer than the best
    route found.
    """
    best = find_route(start, target)
    if not best.steps or limit <= 1:
        return [best]

    max_carbons = carbon_limit(start, target)
    start_state, is_goal = _setup(start, target)
    max_depth = best.length + extra_steps

    found: dict[tuple, Route] = {best.signature: best}
    budget = NODE_BUDGET

    def walk(state: State, trail: list[tuple[State, Reaction]], seen: set[State]) -> None:
        nonlocal budget
        if budget <= 0 or len(trail) >= max_depth:
            return
        for next_state, reaction in _successors(state, max_carbons):
            budget -= 1
            if budget <= 0:
                return
            if next_state in seen:
                continue
            trail.append((next_state, reaction))
            if is_goal(next_state):
                route = _route_from_trail(trail, start_state)
                found.setdefault(route.signature, route)
            else:
                seen.add(next_state)
                walk(next_state, trail, seen)
                seen.discard(next_state)
            trail.pop()

    walk(start_state, [], {start_state})

    routes = sorted(found.values(), key=lambda r: (r.length, r.cost, r.signature))
    return routes[:limit]


# --------------------------------------------------------------------------
# Reachability, for better error messages
# --------------------------------------------------------------------------


def reachable_groups(start_group: str) -> set[str]:
    """Functional groups reachable from ``start_group``, ignoring chain length."""
    seen = {start_group}
    stack = [start_group]
    while stack:
        group = stack.pop()
        for reaction in _BY_START.get(group, []):
            if reaction.end not in seen:
                seen.add(reaction.end)
                stack.append(reaction.end)
    return seen
