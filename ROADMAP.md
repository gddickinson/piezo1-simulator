# ROADMAP — what is not done

Planned work, in ~20-minute rounds. Each round: implement, test, fix, update
the docs, commit. Items are `[ ]` planned, `[~]` in progress, `[x]` done.

**Status: 10 open items across five rounds**, all of them Block Q.
Everything finished — 362 items across 75 rounds, each carrying the result it
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

### Round 76 — The hybrid model should be reachable
- [ ] `structure/hybrid.py` exists and nothing in the GUI or CLI can build one.
      It serves a stated project aim and is currently notebook-only, which is
      the exposure gap Round 58 found for the coupling score.
- [ ] *Validate:* the seam must be visibly rendered — a full-length model whose
      predicted 569 residues look like the experimental ones is precisely the
      confident-wrong-picture failure Round 50 audited for.

### Round 77 — A fetch that verifies what it downloaded
- [ ] Round 60 found a broken CDS endpoint; Round 65 found two 127-byte error
      pages stored as structures. `_download`'s size guard is necessary and not
      sufficient. Verify content type or parse-ability before writing.
- [ ] *Validate:* a planted error page must be rejected, and the test must show
      the guard rejecting something the size check would accept.

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
