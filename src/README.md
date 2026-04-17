# Native Rebuild Skeleton

This directory is the starting point for the new native implementation.

Planned layout:

- `leaf/`
- `semantic/`
- `compile/`
- `eval/direct/`
- `eval/exact/`
- `runtime/`
- `bridge/`

The architecture this tree is intended to support is defined in:

- [dev/test_models/rebuild_design.md](/Users/nstevenson/Documents/2025/AccumulatR/dev/test_models/rebuild_design.md)
- [dev/test_models/rebuild_phase0_decisions.md](/Users/nstevenson/Documents/2025/AccumulatR/dev/test_models/rebuild_phase0_decisions.md)

Phase 1 now defines:

- `leaf/dist_kind.hpp`
- `leaf/channels.hpp`
- `semantic/model.hpp`
- `semantic/model.cpp`

These files only define the base native types and a minimal semantic-model
validator. There is still no evaluator logic here.
