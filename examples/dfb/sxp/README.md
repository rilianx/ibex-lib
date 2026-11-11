SoPlex 4.0.2 `spxsolve.cpp` with two experimental switches (2026-11-11), linked
only into `ibexreplay_sxp` (see the commands in `ESTADO.md`):

- `SOPLEX_SIN_REFACT_OPT=1`: no refactorization when the optimum is reached;
- `SOPLEX_SIN_PRECISION=1`: the precision check (`precisionReached`) is skipped.

With both switches off the code is the original. Build: compile `spxsolve.cpp`
with the flags of SoPlex, replace the member in a copy of `lib/libsoplex.a`
(`ar r libsoplex.a spxsolve.cpp.o`) and link with that copy.
