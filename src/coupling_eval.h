#pragma once

#include <string>

#include "forced_state.h"
#include "coupling_program.h"
#include "trial_params.h"

double evaluate_outcome_coupling_unified(
    const uuber::NativeContext &ctx, const uuber::VectorProgram &program,
    double upper, int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, double rel_tol, double abs_tol,
    int max_depth, bool include_na_donors, int outcome_idx_context,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid);

bool evaluate_outcome_coupling_unified_batch(
    const uuber::NativeContext &ctx, const uuber::VectorProgram &program,
    double upper, int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, double rel_tol, double abs_tol,
    int max_depth, bool include_na_donors, int outcome_idx_context,
    const uuber::BitsetState *forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    bool forced_survive_bits_valid,
    const std::vector<const uuber::TrialParamsSoA *> &trial_params_soa_batch,
    std::vector<double> &out);
