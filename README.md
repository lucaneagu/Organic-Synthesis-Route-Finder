# Organic Synthesis Route Finder

Type in the chemical you have and the chemical you want, and this tells you the
reagents and conditions for every step in between — using only the reactions on
the **AQA A-level Chemistry** specification (3.3 Organic chemistry).

```
$ chemsynth ethene "propanoic acid"

==============================================================================
 ETHENE -> PROPANOIC ACID
 3 steps | carbon chain 2 -> 3
==============================================================================

  STEP 1  ethene -> bromoethane
                CH2=CH2 -> CH3CH2Br
    Reagents    HBr (or HCl) gas
    Conditions  Room temperature
    Type        Electrophilic addition
    Note        Ethene is symmetrical, so there is only one possible product.
    Spec        AQA 3.3.4

  STEP 2  bromoethane -> propanenitrile
                CH3CH2Br -> CH3CH2CN
    Reagents    KCN (or NaCN) dissolved in ethanol
    Conditions  Heat under reflux
    Type        Nucleophilic substitution
    Carbons     chain lengthened by 1 carbon
    Note        The standard way to LENGTHEN a carbon chain by one: the CN
                carbon becomes part of the molecule.
    Spec        AQA 3.3.3

  STEP 3  propanenitrile -> propanoic acid
                CH3CH2CN -> CH3CH2COOH
    Reagents    Dilute HCl (aq)
    Conditions  Heat under reflux
    Type        Hydrolysis
    Note        Ammonium ions are the other product. The chain keeps the extra
                carbon gained from the cyanide.
    Spec        AQA 3.3.9

  OVERALL
    ethene  ->  bromoethane  ->  propanenitrile  ->  propanoic acid
```

## Running it

No installation and no dependencies — just Python 3.10 or newer:

```bash
python3 -m chemsynth ethene "propanoic acid"
```

Run it with no arguments to get an interactive prompt that keeps asking for
pairs of chemicals:

```bash
python3 -m chemsynth
```

Or install it so that `chemsynth` works as a command anywhere:

```bash
pip install -e .
chemsynth benzene phenylamine
```

There is also a **web version**: open [`web/index.html`](web/index.html) in any
browser. It is a single self-contained file with the same reaction data, so it
works offline and on a phone.

That file is generated, not hand-edited — the reaction map, the compound names
and the list of names the search box accepts are all taken from the Python
package, so the two versions cannot drift apart. After changing any chemistry,
rebuild it:

```bash
python3 build_web.py
```

## What you can type

| You can type | Examples |
| --- | --- |
| Systematic names | `propan-1-ol`, `2-bromobutane`, `ethyl propanoate`, `2-hydroxypropanenitrile` |
| Older/common names | `acetone`, `acetic acid`, `toluene`, `aniline` |
| Condensed formulae | `CH3CH2OH`, `CH3COOH`, `CH2=CH2` |
| A whole class of compound | `a primary alcohol`, `an aldehyde`, `a carboxylic acid` |

Names are case-insensitive, and a misspelling gets a "did you mean" suggestion.
`chemsynth --list compounds` prints the full range.

## Options

| Option | What it does |
| --- | --- |
| `-n N`, `--routes N` | Show up to N routes, shortest first, so you can compare alternatives |
| `--json` | Print the route as JSON (for revision cards, spreadsheets, other programs) |
| `--list groups\|compounds\|reactions` | Show what the program knows about |
| `--width COLS` | Change the output width |
| `--no-color` | Turn off colour |

## How it works

The 60 reactions in [`chemsynth/reactions.py`](chemsynth/reactions.py) form a
map. Each one is an arrow from one class of compound to another, carrying the
reagents, conditions, reaction type and AQA specification reference:

```python
Reaction(
    id="haloalkane_cyanide",
    start="haloalkane", end="nitrile",
    reagents="KCN (or NaCN) dissolved in ethanol",
    conditions="Heat under reflux",
    rxn_type="Nucleophilic substitution",
    spec="3.3.3", delta_min=1, delta_max=1,
)
```

Finding a synthesis is then a shortest-path search (Dijkstra) across that map.
The search tracks the **number of carbon atoms** as well as the functional
group, so it knows that KCN and HCN are the two ways to lengthen a chain by one
carbon, and it will not offer you a route that quietly changes the size of the
molecule. Step costs are all roughly 1, nudged upwards for messy reactions like
free-radical substitution, so the route an examiner would expect comes out on
top.

Esters are tracked as two separate halves, because `propyl ethanoate` and
`ethyl propanoate` both have five carbons but are made in completely different
ways. This also stops the search cheating: it will not "make" a molecule by
esterifying with a convenient co-reactant and then hydrolysing the other half
back off.

Halogenoalkanes are tracked as primary (1-bromo) or secondary (2-bromo),
because the position of the halogen decides what you can make. HBr adds to an
unsymmetrical alkene by **Markovnikov's rule**, so propene gives
2-bromopropane, not 1-bromopropane, and that carries through to every later
step: 2-bromopropane hydrolyses to propan-2-ol, not propan-1-ol. Ask for
propene to propan-1-ol and the only honest A-level route goes through propane
and free-radical bromination, which the output flags as low-yield.

## What is on the map

Alkanes, alkenes, halogenoalkanes (primary and secondary), alcohols (primary, secondary, tertiary),
aldehydes, ketones, carboxylic acids and their salts, esters, acyl chlorides,
acid anhydrides, amides, nitriles, hydroxynitriles, amines (primary to
quaternary), benzene, nitrobenzene, phenylamine, alkylbenzenes and aromatic
ketones.

Deliberately **not** on the map, so you get an honest "no route" instead of a
wrong answer:

- **Phenol** — its reactions do not connect to the aliphatic synthesis routes.
- **Side-chain oxidation** of arenes (so no benzoic acid) — not on the AQA spec.
- **Condensation polymers** (polyesters, polyamides) — these need two
  difunctional monomers, which a one-molecule-at-a-time map cannot represent.
- **Tertiary alcohols** are a dead end, because they are not oxidised.
- **KCN or NH3 on a secondary halogenoalkane** — these give branched products
  (2-methylpropanenitrile, propan-2-amine) that the straight-chain naming used
  here cannot describe, so those steps are left out.

Two steps are marked as *beyond the core specification* and flagged in the
output when a route uses them: making an acyl chloride with SOCl2, and making an
amino acid from a carboxylic acid.

## Tests

```bash
python3 -m unittest discover -s tests
```

The tests check the naming and formula generation, the name parser, the data in
the reaction map, and a set of known synthesis routes — including that each
step's product really is the next step's starting material.

## A warning

This is a revision aid, not a lab manual. It gives the reagents and conditions
that AQA expect in an answer; it does not know about safety, yield, purification
or anything that happens in a real fume cupboard. Always check your route
against your own textbook and mark scheme.
