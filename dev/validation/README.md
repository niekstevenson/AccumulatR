# Validation

This folder contains a hand-derived validation harness for the rebuild likelihood engine.

Run it from the repo root with:

```sh
Rscript dev/validation/run_validation.R
```

The default runner compares engine likelihoods against explicit manual references for 24 model shapes:

1. `independent_trigger_two_way`
2. `pool_vs_competitor`
3. `all_of_three_way`
4. `first_of_three_way`
5. `chained_onset_single_outcome`
6. `ranked_independent`
7. `ranked_chained_onset`
8. `shared_trigger_conditioning`
9. `shared_gate_pair`
10. `guarded_positive_mass_tie`
11. `shared_gate_three_way_tie`
12. `nested_guard_pair`
13. `deep_guard_chain`
14. `pooled_shared_gate_tie`
15. `pooled_guarded_shared_gate_tie`
16. `overlapping_composite_competitors`
17. `guarded_overlapping_competitors`
18. `shared_gate_four_way_tie`
19. `none_of_conjunction`
20. `stop_change_shared_trigger`
21. `stim_selective_stop`
22. `stim_selective_stop2`
23. `density_lift_competitor_subset`
24. `partial_overlap_composite_gates`

There is also a heavier adversarial runner:

```sh
Rscript dev/validation/run_adversarial_validation.R
```

To run one adversarial case:

```sh
Rscript dev/validation/run_adversarial_validation.R --case=oracle_nested_choice_guard_absence
```

That runner includes intentionally expensive or currently suspicious compositions:

1. `oracle_repeated_shared_gate_six_way`
2. `oracle_nested_choice_guard_absence`
3. `oracle_deep_composite_blocker`
4. `oracle_pool_k2_shared_gate_guard`

The runner exits nonzero if any check fails. That is intentional. This folder is meant to expose correctness gaps, not hide them.
