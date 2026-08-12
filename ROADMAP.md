# ROADMAP — what is not done

Planned work, in ~20-minute rounds. Each round: implement, test, fix, update
the docs, commit. Items are `[ ]` planned, `[~]` in progress, `[x]` done.

**Status: no open items.** Block R is finished. The next block
has not been written; see `docs/CONCLUSION.md` before adding one.
Rounds **84b** through **84f** were added mid-block on request and completed
out of order; Round 84 itself is still open. Their records are in the archive.
Everything finished — 382 items across 83 rounds, each carrying the result it
measured — is in
[`docs/ROADMAP_COMPLETED.md`](docs/ROADMAP_COMPLETED.md).

Round 75 split the two apart. This file had grown to 2,702 lines of which 96%
was a record of work already done, duplicating `SESSION_LOG.md`, and the thing
it is supposed to answer — *what is left?* — took scrolling to find. Splitting
it also surfaced a defect that length had hidden: **two entries numbered Round
68**, adjacent, one ticked as superseded by Round 63 and one still open with the
same question. Round 67 had recorded the supersession by adding a heading rather
than ticking the original. Merged in the archive, with what happened recorded
rather than tidied away.

Nothing was lost in the split, and that is checked rather than asserted:
`tests/test_roadmap.py` compares both files against the pre-split original in
git, item by item.

**Before adding to this list, read [`docs/CONCLUSION.md`](docs/CONCLUSION.md).**
The variant-effect prediction this project was built for does not work — five
pre-registered tests, five nulls — and a sixth on the same variant set should
not be run whatever predictor goes into it. Items that would reopen that
question need the feasibility argument answered first, not a new predictor.

---

## The destination

Everything here builds toward one capability, which is the reason the project
exists:

> **Predict, from structure alone, whether a PIEZO1 variant is gain- or
> loss-of-function — and validate that prediction blind against the 68 curated
> variants whose phenotypes are already known.**

If elastic-network perturbation at a mutated residue systematically shifts the
predicted gating energetics in the direction of the measured phenotype, that is
a real, testable, publishable result. Every item below is either a component of
that pipeline or a tool for interrogating it.

The chain is: **structure → elastic network → gating coordinate → dome area
change → tension-dependent free energy → open probability → comparison with
measured P50 and inactivation kinetics.** Each round closes one link.

---

## Standing per-round checklist

1. Run the full test suite; fix anything red before adding features.
2. Run the scripted GUI smoke test — mechanical refactors of Qt code have
   broken it twice, silently.
3. Check no file exceeds 500 lines.
4. Update `INTERFACE.md` status, `SESSION_LOG.md` reasoning, `docs/SCIENCE.md`
   if a parameter or source changed, and this file's checkboxes.
5. Commit with a message explaining *why*, not just what.

---

## Deliberately not doing

- All-atom MD in the interactive loop. OpenMM is available for offline
  refinement, but interactive dynamics stay coarse-grained by design.
- Cooperativity between channels: measured P50 and open probability are
  invariant from 1 to 100 channels/µm² (Lewis & Grandl 2021), so modelling it
  would be inventing physics.
- Modelling the gain-of-function activation latency (344 ± 133 ms) inside the
  Markov scheme — the authors state explicitly that no Markov model reproduces
  it.
- Full C3 block-diagonalisation of the Hessian: benchmarked at only 1.76× on
  top of the sparse solve, not worth the complexity. Symmetry *labelling* gives
  the scientific payoff already.

## Block Q — what the last five rounds exposed (Rounds 76–80)

---

## Block R — what a comparable research tool has and this does not (Rounds 81–85)

*Added 2026-08-10, after reviewing the codebase against HOLE, CHAP, MOLEonline,
ProDy and the APBS/PDB2PQR electrostatics route. Three of the five items below
are gaps this project can state as measurements rather than opinions: an API
argument nothing supplies, a field on every structure nothing reads, and a
downloaded entry nothing compares.*

### Round 81 — The pore has no charge, and the API already says so
- [x] `solve_pnp` takes a `fixed_charge` argument, its documented equation
      carries a ρ_fixed term, and **no caller anywhere has ever supplied one**.
      Every permeation number this project has produced treats the pore as
      electrically neutral. Meanwhile `functional_residues.json` curates four
      sequence-verified acidic residues as *"acidic residues setting ion
      selectivity"* — human E2117, E2461, E2469, E2470, so twelve charges
      across the trimer — and `default_species()` offers a generic "cation" and
      "anion" that are not even named Na⁺, K⁺ and Ca²⁺. A cation channel
      modelled without charge cannot be selective, and this one is not.
      **Done.** `physics/pore_charge.py` builds the density; the species are
      now K⁺/Cl⁻/Ca²⁺ and a NaCl set exists for the published protocol.
- [x] Map the curated acidic and basic pore-lining residues onto the axial
      slices the pore profile already produces, and pass the result through.
      **Done, and it measured something on the way:** three of the four curated
      selectivity glutamates are not within side-chain reach of the lumen
      (E2117 at 12.9 Å past the wall), which is what Coste et al. 2015
      concluded about E2117 from function alone. Two routes are reported side
      by side — curated (6 charges, net −6) and geometric (46, net **+8**) —
      because they disagree in kind.
- [x] *Validate:* the charged pore must make the model **cation-selective**, as
      a permeability ratio compared with the published PIEZO1 value — and if it
      does not, that is the result and it gets reported. Zero fixed charge must
      reproduce today's numbers **exactly**, so the change is visible,
      reversible, and cannot silently move the existing conductance claim.
      **Yes, in direction; no, in value.** P_Cl/P_Na = 0.021 curated and 0.207
      geometric against a published 0.14, from an uncharged baseline of 0.904
      — the two routes bracket the measurement tenfold apart rather than
      reproducing it, and the curated route only gets there at an in-pore
      concentration of 13.9 M, which is flagged as outside the continuum
      model's validity. An explicitly zero charge reproduces the neutral pore
      bit for bit. Finding this required correcting an **inverted drift sign**
      in the Scharfetter-Gummel flux that no test could see while every
      current was computed between identical baths; the correction moved the
      recorded conductance by one part in 10¹⁴.

### Round 82 — The B-factors every structure carries, and nothing reads
- [x] `Structure` parses `b_factor` for every atom and no analysis has ever
      looked at one. Comparing predicted mean-square fluctuation against
      observed B-factor is the standard validation of an elastic network — the
      one ProDy users run first — and this project has never run it on the
      network its central mechanism claim rests on. There is no GNM here
      either; the ANM has no scalar-fluctuation counterpart.
      **Done, with one correction to the premise:** the *prediction* did
      exist — `ModeSet.msf` sums `|v|²/λ` and the feature table already
      consumes it. What was missing was the comparison, so this round is a
      missing validation rather than missing physics.
      `analysis/fluctuations.py`, reachable from the CLI and the GUI.
- [x] *Validate:* report the correlation **including if it is poor**, and
      calibrate before believing it: cryo-EM B-factors are sharpened, per-atom
      values in a backbone-only model may be uniform, and some entries carry
      pLDDT rather than B in that column. Establish on a case with a known
      answer that the comparison can distinguish a good network from a
      deliberately bad one, or record why the comparison cannot be made here.
      **All three, and the answer is split.** Calibrated on a planted
      fluctuation that the right network recovers exactly and one built on
      shuffled coordinates does not. All three bad-column cases are refused
      with the reason — and the pLDDT gate is shown to point the right way by
      measurement, not assertion: AlphaFold's own column anti-correlates at
      −0.57. Measured across the catalogue: **18 of 21 entries can answer**,
      median Spearman **0.74** against a contact-number control's **0.32**,
      winning 13 of 15 — but on Pearson **0.48 against 0.39**, winning only 9
      of 15. The network orders residues by mobility much better than burial
      does and predicts how far they move barely better; both are reported.
      Three entries have a *negative* control, meaning their B-factor rises
      with burial and is not a mobility at all.

### Round 83 — PIEZO2 is downloaded, classified, and only ever excluded
- [x] 6KG7 is fetched, entity-classified and then excluded from every ensemble
      as a paralogue. That exclusion is correct for a PIEZO1 ensemble and wrong
      as a final answer. PIEZO2 is the tactile paralogue with different
      inactivation kinetics, and it is the obvious control for the question
      this project never asks: **how much of this mechanism is PIEZO1, and how
      much is the fold?**
      **Done — and the registry was wrong about the entry.** Its note said
      6KG7 "resolves residues 8-823"; it resolves **8–2822 in 16 segments,
      1,817 C-alphas per protomer**, more than any PIEZO1 entry here and
      including all 38 TM helices. Corrected, and pinned against the file.
      `analysis/paralogue.py`, reachable from the CLI and the GUI.
- [x] *Validate:* the dome geometry and the gating-mode symmetry analysis must
      run on PIEZO2 and be reported beside PIEZO1. If the two are
      indistinguishable, say so plainly — that is a result about generality,
      not a failure. Cross-species numbering must go through the alignment, and
      PIEZO2 is 2,752 aa, so a constant offset is certain to be wrong.
      **They are indistinguishable, and saying so needed the coverage
      matching.** Measured naively PIEZO2's dome is 8.5 nm deep against
      PIEZO1's 4.9 — a **coverage artefact**, since 6KG7 resolves 38 helices
      against 22. Matched to the shared helices it gives 5.6 nm and R_c
      10.32 against 9.72, inside the PIEZO1 range on every quantity. The
      gating mode is the fold's: overlap **0.804** with one PIEZO2 symmetric
      mode, **0.925** of it inside PIEZO2's symmetric subspace, against a
      shuffled-correspondence control of 0.190. Numbering is measured rather
      than assumed — each entry matches one of four UniProt references at
      1.000 — and 6KG7 is **mouse** Piezo2 (Q8CD54, 2,822 aa), not the human
      2,752 the roadmap assumed. The TM-index pairing is confirmed by the
      alignment for 37 of 38 helices, and the protomer order is **(2, 0, 1)**,
      so chain labels would have been wrong.

### Round 84 — Nothing computed can leave the application
- [x] Conservation, mechanical coupling, PRS response, mode displacement and
      the wetting score are all per-residue scalars, and none of them could be
      exported in a form another viewer can colour by. **Done:**
      `core/export.py`, reachable as `piezo1.cli export` and from
      **File → Export coloured structure**, which writes whatever scalar is
      currently painted on the model.
- [x] *Validate:* the round trip checked **numerically**, atom by atom against
      the source array. **Done, and it caught two things.** Unscored atoms go
      out with **occupancy 0.00** rather than a zero score, and the GUI exports
      the raw residue map rather than `view.values`, whose unmeasured entries
      are filled to the map floor for display. And the B-factor field is six
      characters, so `-999.99` — the sentinel this first used — **overflows it
      and shifts every column after it**, producing a file that still parses
      and is wrong. The limit is now derived by formatting a value and
      measuring it rather than remembered.

### Rounds 84b–84f — requested mid-block, completed out of order

Five rounds added on request while Round 84 above was still open: replicating
Guo & MacKinnon 2017 panel by panel, fixing an ion animation that had never
drawn an ion, replicating Liu et al. 2025 and making the conduction pathway a
choice, a component viewer, and correcting how the two halves of the wetting
verdict are composed. **Their records are in
[`docs/ROADMAP_COMPLETED.md`](docs/ROADMAP_COMPLETED.md)**, with the reasoning
in `SESSION_LOG.md`.

### Round 85 — Review after Rounds 81–84
- [x] Five-round review. **Done, and it found a hole.** The discipline held for
      positives — every new instrument was calibrated first and each calibration
      caught something. It did **not** hold for negatives: Round 84d recorded
      "the lateral route does not separate open from closed" as an honest null
      and pinned it with a test, and it was an artefact of its own composition.
      Every guard in this repository is aimed at not over-claiming; nothing
      interrogates a "no". The rule added: **a null needs a positive control**,
      an input on which the instrument must return "yes", exactly as a checker
      needs an input on which it must say "no". Written up in
      `docs/ROADMAP_COMPLETED.md` and `docs/METHODS_NOTE.md`.

### Round 86 — Two deposited entries are not in the numbering we read them in
*Found by Round 83's identification instrument, which reports the shift that
would repair what it finds. Nothing applied it.*
- [x] **6LQI** is the Piezo1.1 isoform, deposited in the isoform's own
      continuous numbering across its 1382-1405 deletion. **Corrected:
      identity 0.447 -> 1.000 with a +24 shift over 765 residues.**
- [x] **8ZU3, 8YFC, 9VMX and 8YFG** carry residues **767-857 numbered 22 low**.
      **Corrected: 0.932 -> 1.000 on three of them, and 0.999 on 8YFG**, which
      is right — it carries the R2456H substitution, a genuine residue change
      that a numbering fix must not absorb.
- [x] *Validate:* full agreement on every corrected entry, **and the null**.
      8YEZ resolves the same 767-857 region correctly and must come back
      untouched; 7WLT, 3JAC and 6B3R likewise. A renumberer is a rewriter, and
      one that touched a correct file would corrupt it.
- [x] Say which published numbers move. **Done, and the answer is specific:
      the dome moves and the pore does not.** 8ZU3/8YFC/9VMX go 12.50 -> 11.91
      nm, 8YFG 11.12 -> 10.77, 6LQI 9.35 -> 11.03, and the transmembrane
      helices found go 25 -> 26 (23 -> 26 on 6LQI) because the correction
      recovers helices that were being missed. The pore bottleneck does not
      move by an Angstrom — it is geometric and reads no residue number. **No
      frozen claim uses any affected entry**, so nothing is superseded, and a
      test fails if one ever starts.

