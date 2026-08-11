# ROADMAP — what is not done

Planned work, in ~20-minute rounds. Each round: implement, test, fix, update
the docs, commit. Items are `[ ]` planned, `[~]` in progress, `[x]` done.

**Status: 6 open items in Block Q, plus a new Block R.**
Everything finished — 368 items across 77 rounds, each carrying the result it
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

### Round 78 — Retire `HALOTAG_CALCIUM_PLAN.md`
- [ ] It is marked 📋 and describes work Rounds 29–32 completed. A plan document
      that outlived its execution is the documentation equivalent of the four
      module rows Round 65 deleted.
- [ ] *Validate:* nothing may be lost — anything in it not carried by the
      implemented modules must move to `SCIENCE.md` first.

### Round 79 — `ARCHITECTURE.md`, or the row goes
- [ ] The last 📋 in INTERFACE. Either write why the code is shaped this way —
      the dependency arrow, the impostor rendering, structure-of-arrays — or
      delete the promise as Round 65 deleted four others.
- [ ] *Validate:* if written, it must not restate INTERFACE; it is the *why*,
      and a test should check it cites the constraints rather than the contents.

### Round 80 — What a reader should be told first, measured
- [ ] Five surfaces now state the record (README, CONCLUSION, tour, help,
      SCIENCE). Round 59 linked three by a test. Extend it to all five, and
      measure how many clicks or scrolls a new reader needs to reach the
      conclusion from each entry point.
- [ ] *Validate:* the answer must be one step from every entry point, or the
      surface is wrong rather than the reader.

---

## Block R — what a comparable research tool has and this does not (Rounds 81–85)

*Added 2026-08-10, after reviewing the codebase against HOLE, CHAP, MOLEonline,
ProDy and the APBS/PDB2PQR electrostatics route. Three of the five items below
are gaps this project can state as measurements rather than opinions: an API
argument nothing supplies, a field on every structure nothing reads, and a
downloaded entry nothing compares.*

### Round 81 — The pore has no charge, and the API already says so
- [ ] `solve_pnp` takes a `fixed_charge` argument, its documented equation
      carries a ρ_fixed term, and **no caller anywhere has ever supplied one**.
      Every permeation number this project has produced treats the pore as
      electrically neutral. Meanwhile `functional_residues.json` curates four
      sequence-verified acidic residues as *"acidic residues setting ion
      selectivity"* — human E2117, E2461, E2469, E2470, so twelve charges
      across the trimer — and `default_species()` offers a generic "cation" and
      "anion" that are not even named Na⁺, K⁺ and Ca²⁺. A cation channel
      modelled without charge cannot be selective, and this one is not.
- [ ] Map the curated acidic and basic pore-lining residues onto the axial
      slices the pore profile already produces, and pass the result through.
- [ ] *Validate:* the charged pore must make the model **cation-selective**, as
      a permeability ratio compared with the published PIEZO1 value — and if it
      does not, that is the result and it gets reported. Zero fixed charge must
      reproduce today's numbers **exactly**, so the change is visible,
      reversible, and cannot silently move the existing conductance claim.

### Round 82 — The B-factors every structure carries, and nothing reads
- [ ] `Structure` parses `b_factor` for every atom and no analysis has ever
      looked at one. Comparing predicted mean-square fluctuation against
      observed B-factor is the standard validation of an elastic network — the
      one ProDy users run first — and this project has never run it on the
      network its central mechanism claim rests on. There is no GNM here
      either; the ANM has no scalar-fluctuation counterpart.
- [ ] *Validate:* report the correlation **including if it is poor**, and
      calibrate before believing it: cryo-EM B-factors are sharpened, per-atom
      values in a backbone-only model may be uniform, and some entries carry
      pLDDT rather than B in that column. Establish on a case with a known
      answer that the comparison can distinguish a good network from a
      deliberately bad one, or record why the comparison cannot be made here.

### Round 83 — PIEZO2 is downloaded, classified, and only ever excluded
- [ ] 6KG7 is fetched, entity-classified and then excluded from every ensemble
      as a paralogue. That exclusion is correct for a PIEZO1 ensemble and wrong
      as a final answer. PIEZO2 is the tactile paralogue with different
      inactivation kinetics, and it is the obvious control for the question
      this project never asks: **how much of this mechanism is PIEZO1, and how
      much is the fold?**
- [ ] *Validate:* the dome geometry and the gating-mode symmetry analysis must
      run on PIEZO2 and be reported beside PIEZO1. If the two are
      indistinguishable, say so plainly — that is a result about generality,
      not a failure. Cross-species numbering must go through the alignment, and
      PIEZO2 is 2,752 aa, so a constant offset is certain to be wrong.

### Round 84 — Nothing computed can leave the application
- [ ] Conservation, mechanical coupling, PRS response, mode displacement and
      the wetting score are all per-residue scalars, and none of them can be
      exported in a form another viewer can colour by. `to_pdb` writes
      coordinates only. The standard interop route — write the scalar into the
      B-factor column and open it in PyMOL or ChimeraX — is about twenty lines
      and is missing, so every result this project computes is trapped inside
      it or inside a JSON blob.
- [ ] *Validate:* check the round trip **numerically**, not by eye — read the
      written file back and compare the B-factor column against the source
      array element by element. A file that merely opens proves nothing. State
      what is lost: unmeasured residues must be distinguishable from a genuine
      zero, which is exactly the trap `analysis_controller` already handles by
      using the map floor rather than zero.

### Round 85 — Review after Rounds 81–84
- [ ] Five-round review. The standing question for this block: the project has
      spent seventy rounds proving what it *cannot* establish, and these four
      items are the first in a while that could each return a positive result.
      Check that the same discipline held — pre-registration where a comparison
      is involved, a calibrated instrument before any cross-check is believed,
      and a null reported as a null.
