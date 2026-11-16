# Reproducing *A self-certifying simplex for polytope-hull contraction in interval branch and bound*

I. Araya, E. Pérez, V. Reyes.

DOI: <https://doi.org/10.5281/zenodo.24284785>

This archive is a snapshot of the branch `dfb-sc-paper` of
<https://github.com/rilianx/ibex-lib>, a fork of Ibex 2.8.9. It contains
everything needed to rebuild the code of the paper and to recompute its tables
and figure from the raw results.

## Where things are

| Path | Contents |
|---|---|
| `examples/dfb/` | the research code: the contractor `CtcDFBHull`, the simplex SC (`ibex_DFBSimplex*`, `ibex_DFBBasis*`), the drivers, and `build_paper.sh` |
| `examples/dfb/paper/` | the paper (`main.tex`, `main.pdf`) and the data of its figure |
| `results/escalera_2026-11-06/` | raw results of every experiment, the scripts that ran them, the scripts that turn them into the tables, and the banks of passes |
| `examples/sc_plugin/` | `CtcSC`, a standalone version of SC that builds against an installed Ibex (also at <https://github.com/rilianx/ibex-sc>) |
| everything else | Ibex 2.8.9 with the fork's changes, built in place by the experiments |

## Steps

1. **Build** Ibex in place and the binaries of the paper (GCC, CMake; GAOL
   4.2.0 and SoPlex 4.0.2 are bundled):

       cmake -DCMAKE_BUILD_TYPE=Release -DINTERVAL_LIB=gaol -DLP_LIB=soplex . && make
       examples/dfb/build_paper.sh

2. **Unpack the data** and **recompute the tables**: follow
   `results/escalera_2026-11-06/README.md`. It maps every table, figure and
   side experiment of the paper to the script that produces it, and explains
   the two sets of runs (`work/S` and `work/S100`, the latter used by the
   paper).

3. **Rerun the experiments** (optional, several days of CPU): the run scripts
   are in the same folder and are described in its README. Runs with the
   `compo` relaxation on systems of equations are not reproducible bit for
   bit, because SC then chooses its warm start by measured CPU time; the banks
   of passes are.

## License

Ibex and the code of this fork are distributed under the LGPL-3.0 (see
`LICENSE`). SoPlex and GAOL are bundled under their own licenses (see their
archives in `lp_lib_wrapper/` and `interval_lib_wrapper/`). The results and
scripts in `results/` may be reused under CC BY 4.0.
