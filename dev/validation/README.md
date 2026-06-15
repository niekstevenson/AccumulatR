# Validation

This folder contains a hand-derived validation harness for the rebuild likelihood engine.

Run it from the repo root with:

```sh
Rscript dev/validation/run_validation.R
```

The default runner compares engine likelihoods against explicit manual references for 30 model shapes:

1. `independent_trigger_two_way`
2. `censor_trunc_independent_two_way`
3. `truncated_latent_mixture_ratio`
4. `pool_vs_competitor`
5. `all_of_three_way`
6. `first_of_three_way`
7. `chained_onset_single_outcome`
8. `ranked_independent`
9. `ranked_chained_onset`
10. `shared_trigger_conditioning`
11. `stop_change_shared_trigger`
12. `stim_selective_stop`
13. `stim_selective_stop2`
14. `shared_gate_pair`
15. `guarded_positive_mass_tie`
16. `shared_gate_three_way_tie`
17. `nested_guard_pair`
18. `deep_guard_chain`
19. `pooled_shared_gate_tie`
20. `pooled_guarded_shared_gate_tie`
21. `density_lift_competitor_subset`
22. `overlapping_composite_competitors`
23. `guarded_overlapping_competitors`
24. `shared_gate_four_way_tie`
25. `none_of_conjunction`
26. `first_of_absence_choice`
27. `guarded_composite_vs_guarded_competitor`
28. `composite_blocker_guard`
29. `partial_overlap_composite_gates`
30. `nested_choice_guard_absence`

There is also a heavier adversarial runner:

```sh
Rscript dev/validation/run_adversarial_validation.R
```

To run one adversarial case:

```sh
Rscript dev/validation/run_adversarial_validation.R --case=oracle_deep_composite_blocker
```

That runner includes intentionally expensive or currently suspicious compositions:

1. `oracle_repeated_shared_gate_six_way`
2. `oracle_deep_composite_blocker`
3. `oracle_pool_k2_shared_gate_guard`

The runner exits nonzero if any check fails. That is intentional. This folder is meant to expose correctness gaps, not hide them.
