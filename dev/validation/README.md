# Validation

This folder contains a hand-derived validation harness for the rebuild likelihood engine.

Run it from the repo root with:

```sh
Rscript dev/validation/run_validation.R
```

The default runner compares engine likelihoods against explicit manual references for 28 model shapes:

1. `independent_trigger_two_way`
2. `pool_vs_competitor`
3. `all_of_three_way`
4. `first_of_three_way`
5. `chained_onset_single_outcome`
6. `ranked_independent`
7. `ranked_chained_onset`
8. `shared_trigger_conditioning`
9. `stop_change_shared_trigger`
10. `stim_selective_stop`
11. `stim_selective_stop2`
12. `shared_gate_pair`
13. `guarded_positive_mass_tie`
14. `shared_gate_three_way_tie`
15. `nested_guard_pair`
16. `deep_guard_chain`
17. `pooled_shared_gate_tie`
18. `pooled_guarded_shared_gate_tie`
19. `density_lift_competitor_subset`
20. `overlapping_composite_competitors`
21. `guarded_overlapping_competitors`
22. `shared_gate_four_way_tie`
23. `none_of_conjunction`
24. `first_of_absence_choice`
25. `guarded_composite_vs_guarded_competitor`
26. `composite_blocker_guard`
27. `partial_overlap_composite_gates`
28. `nested_choice_guard_absence`

There is also a heavier adversarial runner:

```sh
Rscript dev/validation/run_validation.R --adversarial
```

To run one adversarial case:

```sh
Rscript dev/validation/run_validation.R --adversarial --case=oracle_deep_composite_blocker
```

That runner checks complex compositions against independent density,
order-statistic, and shared-gate formulas:

1. `oracle_repeated_shared_gate_six_way`
2. `oracle_deep_composite_blocker`
3. `oracle_pool_k2_shared_gate_guard`

It exits nonzero if any check fails.

Compiler structure and complexity budgets are checked separately with:

```sh
Rscript dev/validation/compiler_architecture_acceptance.R
```
