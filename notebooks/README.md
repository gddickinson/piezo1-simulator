# Notebooks

Four worked examples of the headless engine. Read them in order the first time;
each one stands alone afterwards.

| | | Runs in |
|---|---|---|
| [`01_first_look`](01_first_look.ipynb) | What is in a deposited structure, how to frame it, and measuring the membrane dome | ~1 s |
| [`02_gating_motion`](02_gating_motion.ipynb) | The elastic network model, and the symmetry rule that says which motions can couple to tension | ~15 s |
| [`03_pore_to_current`](03_pore_to_current.ipynb) | Is the pore open, would water stay in it, and what current would flow | ~5 s |
| [`04_variants_and_the_null`](04_variants_and_the_null.ipynb) | The variant workflow, and the result that did not work | ~10 s |

## Running them

```bash
python -m piezo1.io.fetch          # the structures and the hydration grid
pip install jupyterlab             # or: pip install -e ".[notebooks]"
jupyter lab notebooks/
```

Everything is Qt-free and needs no display.

## Why they ship without outputs

A committed cell output is a number nobody checks. It goes stale silently,
because nothing recomputes it, and a stale number that looks authoritative is
the failure this project spends most of its machinery avoiding.

So the notebooks carry no stored outputs, and instead **assert** the numbers
they quote. Running one checks the science rather than only the syntax: if the
dome radius drifts, or a symmetry-forbidden mode stops scoring zero, the
notebook stops at that cell.

## How they are maintained

They are generated, not hand-edited:

```bash
python scripts/build_notebooks.py           # verify and write
python scripts/build_notebooks.py --check   # verify only
make notebooks
```

The cell content lives in `scripts/notebook_content.py` and
`scripts/notebook_content_analysis.py` as ordinary Python, so it can be
reviewed and diffed instead of being buried in JSON. **Every code cell is
executed, in order, in one namespace, before anything is written** — a notebook
that raises is not published.

Editing a `.ipynb` here directly will be overwritten. Edit the content module.

`tests/test_notebooks.py` executes the committed files, which is what a reader
actually runs, and checks they have not drifted from the content modules.

## Before you build on this

Notebook 04 walks the variant workflow to the point where it stops working.
The one-page account is [`../docs/CONCLUSION.md`](../docs/CONCLUSION.md): the
structural machinery reproduces the literature, the variant-effect prediction
does not work, and five pre-registered tests returned five nulls.

## See also

`docs/NOTEBOOK.md` is the full API reference, including a table of the things
that will bite you — non-constant residue numbering, unreliable chain labels,
and the difference between the two ways a pore can be shut.
