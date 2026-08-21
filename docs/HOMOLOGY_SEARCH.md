# Why this application does not run a homology search

*A decision recorded like a result, in the form `NOT_PREREGISTERED_ROUND64.md`
uses: what was asked, what was measured, what was decided, and what would
reopen it.*

---

## The question

*Should the application build in BLAST searches, so a user can discover PIEZO
homologues and learn about protein structural relationships?*

Asked directly. It deserves a measured answer rather than a preference, because
both possible answers are defensible in the abstract and the project has a rule
about instruments that return plausible numbers.

## The answer

**No BLAST client, and no remote homology search of any kind. Yes to the
statistic BLAST reports, which was the part actually missing.**

Three findings support that, in order of weight.

---

## 1. The family is enumerable, so search is the wrong tool

A single UniProt query — `reviewed:true AND protein_name:piezo` — returns
**nine** entries, in under a second, and that is the whole reviewed family:

| key | accession | protein | organism | length | TM |
|---|---|---|---|---|---|
| `human` | Q92508 | PIEZO1 | *Homo sapiens* | 2521 | 38 |
| `mouse` | E2JF22 | PIEZO1 | *Mus musculus* | 2547 | 38 |
| `rat` | Q0KL00 | PIEZO1 | *Rattus norvegicus* | 2535 | 38 |
| `human_piezo2` | Q9H5I5 | PIEZO2 | *Homo sapiens* | 2752 | 38 |
| `mouse_piezo2` | Q8CD54 | PIEZO2 | *Mus musculus* | 2822 | 38 |
| `worm_piezo` | A0A061ACU2 | PEZO-1 | *C. elegans* | 2442 | 36 |
| `fly_piezo` | M9MSG8 | dPIEZO | *D. melanogaster* | 2551 | 40 |
| `plant_piezo` | F4IN58 | AtPIEZO | *A. thaliana* | 2462 | 35 |
| `dicty_piezo` | Q54S52 | pzoA | *D. discoideum* | 3080 | 35 |

BLAST is the right tool when the answer set is unknown. Here it is known, and
replacing a nine-row pinned table with a live search would trade something
reproducible for something that is not. Every number this project records
carries provenance and must survive a cold clone; a pinned accession does, and
a hit list that grows with the database does not.

All nine are now committed resources, distilled by the same script from the
same source, which is what makes a cross-species helix comparison a measurement
rather than a comparison of two annotation conventions.

## 2. Percent identity is nearly all noise at the distances that matter

This is the measurement that settles it, and it goes the opposite way to the
intuition that a search would help.

Each pair was aligned with BLOSUM62 and then re-aligned against **shuffles of
the same partner that preserve its amino-acid composition exactly** — what the
aligner would score with no homology at all. 20 replicates, registered
parameters, `piezo1.analysis.homology`.

| pair | identity | shuffled null | z | local score | null | z |
|---|---|---|---|---|---|---|
| human vs mouse | 0.831 | 0.233 | 63.9 | 10793 | 56 | 1277 |
| human vs PIEZO2 | 0.495 | 0.220 | 29.1 | 5707 | 50 | 1392 |
| human vs PEZO-1 | 0.317 | 0.221 | 10.2 | 2069 | 51 | 408 |
| human vs AtPIEZO | 0.263 | 0.222 | **6.0** | 627 | 52 | **85** |
| human vs pzoA | 0.238 | 0.208 | **3.5** | 647 | 47 | **102** |
| rat vs pzoA | 0.235 | 0.210 | **2.6** | 663 | 47 | **145** |
| PEZO-1 vs AtPIEZO | 0.238 | 0.225 | **1.5** | 391 | 51 | **64** |

Read the last row. PEZO-1 and Arabidopsis PIEZO are **23.8% identical**, and a
*scrambled* Arabidopsis sequence of the same composition is **22.5% identical**
to PEZO-1. The percentage carries essentially nothing. The local alignment
score on the identical pair is 64 standard deviations above its null.

**17 of the 45 pairs fall below Rost's 30% twilight line** (15 of 36 before
piezo3 joined the family in Round 93). Two of them have an identity
statistically indistinguishable from chance while their alignment score is
overwhelming.

The conclusion is not "these proteins are only distantly related". It is that
**percent identity is the wrong statistic below that line**, which is precisely
what the line means, and an application reporting one without the other would
be handing users a confident wrong reading — the exact failure this project
guards against everywhere else.

So the module refuses to report either statistic alone, and `Relationship.verdict`
says in words which of the two a given pair is entitled to.

## 3. The structural route is the one this project actually has

The generality question — *is the gating mechanism PIEZO1's, or the fold's?* —
is answered far better here by structure than by sequence, and this project
already owns that machinery. Generalising `analysis.paralogue` past its
PIEZO1-vs-PIEZO2 restriction (`analysis.homology_structure`) gives:

Run over **every** pair of deposited entries — three PIEZO1 entries against
each available partner — rather than one representative pair:

| partner | sequence identity | gating-mode overlap, 3 PIEZO1 entries | beats control | |
|---|---|---|---|---|
| PIEZO2 (6KG7, 9VEE) | 0.48 | **0.80 – 0.98** | 6 / 6 | stable |
| PEZO-1 (9UOY, 9ZIS) | 0.29 | **0.18 – 0.98** | 5 / 6 | not stable |
| dPIEZO (9W7X) | 0.30 | **0.19 – 0.98** | 2 / 3 | not stable |

**The range is the result, and getting to it required not stopping at the first
pair.** Comparing 7WLT with 9W7X gives 0.980, and it is tempting: at 30%
sequence identity, *Drosophila* PIEZO appearing to carry PIEZO1's gating
coordinate almost exactly is a striking sentence. Comparing **8YEZ** with the
same 9W7X gives **0.189**. Same two proteins, same method, same code. Which
deposited entry is used decides the number, so no single pair is a property of
the proteins and the first one run would have gone into this document as one.

PIEZO2 is the positive control that makes this a finding rather than a broken
instrument: six pairs, all between 0.80 and 0.98, every one clearing its
shuffled control. The method can say *stable*. It does not say it for the
invertebrates.

What survives is weaker than the cherry-picked row and still worth having: a
symmetric low-frequency mode of PIEZO1's network has a high-overlap counterpart
in *some* pairing with every catalogued homologue, and the paralogue result is
robust. And the coverage-matched dome radius of curvature is 9.25 nm for mouse
PIEZO1 against **9.24 nm** for *C. elegans* PEZO-1.

A sequence search cannot produce any of those rows — including the negative
ones. That is the sense in which this application should help a user "learn
about protein structural relationships": by measuring the structures, and by
saying when the measurement does not hold still.

---

## What was built instead

| Instead of | Built |
|---|---|
| A BLAST client | A pinned nine-member family, all committed, all regenerable |
| An E-value | A composition-matched shuffled null on **both** statistics, and a verdict naming which one a pair supports |
| A hit list | `analysis.homology_structure`, comparing any two catalogued homologues by dome, helix correspondence and mode overlap |
| Nothing | `analysis.homology_sites` — every curated functional residue looked up across all nine, gated by whether the alignment is in register there |

## What this does *not* claim

- **Not that BLAST is a bad tool.** BLAST reports a bit score and an E-value
  rather than a percentage for exactly the reason measured above. The criticism
  here is of the percentage, not of BLAST.
- **Not that no PIEZO homologue remains to be found.** The reviewed set is nine;
  TrEMBL holds many more, including the zebrafish `si:dkey-11f4.7` entries that
  a third vertebrate paralogue would be built from. Those are unreviewed, of
  varying length, and cataloguing them is a curation job with a provenance gate,
  not a search feature. **One of them has since been curated — see the
  Round 93 note below.**
- **Not that the family cannot grow.** `io.fetch.FAMILY_ACCESSIONS` is a pinned
  list precisely so that growth is something a person notices and curates.

## What would reopen this

1. **A reviewed PIEZO appearing that is not in the nine.** Re-running the query
   is a one-line check; adopting the result is a curation decision with the same
   provenance gate every other resource has.
2. **A search over something other than sequence.** Foldseek searches structure,
   which is the axis where this family actually has signal — the sister project
   `piezo_genes` has a working client. Worth revisiting if a *structural*
   neighbour outside the family is the question; it is a different question from
   the one asked here.
3. **A curated PIEZO3.** Dong et al.'s zebrafish paralogue would be the first
   new vertebrate PIEZO since PIEZO2 and would change what the family means. It
   needs the annotation discrepancy resolved first — see below. ✅ **This fired
   in Round 93; the record is the next section.**

## Round 93 — the third clause fired, and what it changed

Standing instruction 3 above was written expecting exactly this, so this is a
record rather than a revision.

The `piezo_genes` census settled the annotation discrepancy that the clause made
a precondition: **piezo3 is a real vertebrate paralogue**, as old as PIEZO1 and
PIEZO2 (both duplications on the jawed-vertebrate stem, ~460–560 Ma), transcribed
and correctly spliced, under purifying selection at its pore, and lost eleven
times since — including on the primate stem, which is why the human genome has
the pseudogene `PIEZO1P2` at the locus and no gene. `si:dkey-11f4.7`
(A0AB32U1Q1) is now the tenth committed reference.

**Three things this changes, and one it does not.**

- **The family is ten and the query returns nine.** The sentence carrying the
  argument above — one UniProt query returns exactly the family — was true when
  every member was reviewed. piezo3 is TrEMBL, so `reviewed:true AND
  family:piezo` still returns nine and now *misses a genuine vertebrate
  paralogue*. `analysis.homology.reviewed_family()` is the nine, kept separate
  from `family()` so the gap stays visible instead of being absorbed into a
  bumped count.
- **That does not weaken the enumerability argument; it sharpens it.** The
  problem the census names is not that the family is too large to enumerate — it
  is that the databases are *incomplete*, and a similarity search over the same
  databases inherits exactly that incompleteness. A BLAST run would not have
  found piezo3 either, for the same reason the reviewed query does not: the
  sequence is there and the *record* is the thing that is missing. What found it
  was a genome sweep with synteny, a phylogeny and a loss reconstruction — not a
  search feature.
- **An unreviewed member's annotation is not architecture.** piezo3's UniProt
  entry names **21** transmembrane helices where its two siblings have 38. That
  is an automatic annotation on an uncurated record, and `FamilyMember.reviewed`
  is carried so nothing reads it as a statement about the protein. A test pins
  it.
- **Unchanged: this project still ships no homology search.** All three grounds
  in the argument above hold, and the third — that the structural route answers
  the generality question better — is now stronger, because piezo3 arrived with
  a structure and `analysis.core_periphery` could measure it.

## A note owed to `piezo_genes`

The sister project at `../piezo_genes` was built around a lesson from Dong et
al.: databases disagree, and the disagreement is the signal. Its own run
flagged that AlphaFold reported a 709-residue maximum for PIEZO2 while every
sequence database said ~2,800.

That lesson paid off directly here. `io.fetch.fetch_alphafold` took
`entries[0]` from the AlphaFold API. For **Q9H5I5 (human PIEZO2)** the endpoint
returns two entries, neither of them canonical — isoform 3 (**709 aa**) and
isoform 2 (2,689 aa) — and `entries[0]` is the 709-residue one. It arrives as a
well-formed mmCIF of the right protein and nothing downstream could have told.
For **A0A061ACU2 (PEZO-1)** the endpoint returns **twelve** entries, one per
annotated isoform; the canonical happened to be first, so this was right by
luck rather than by construction.

The fetcher now selects by exact accession and **refuses** when there is no
canonical model, naming the isoforms it was offered. Human PIEZO2 consequently
has no AlphaFold model in this project, which is the correct answer.

---

*Recorded 2026-08-12. Numbers reproduce with `python -m piezo1.cli homology`
and `python -m piezo1.cli homology --structural 7WLT 9W7X` at the default
parameter registry.*
