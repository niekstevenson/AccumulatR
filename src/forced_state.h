#pragma once

#include <Rcpp.h>

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "bitset_state.h"
#include "context.h"

std::vector<int> forced_vec_from_sexp(SEXP vec);
int forced_bit_capacity(const uuber::NativeContext &ctx);
void build_forced_bitset_strict(const uuber::NativeContext &ctx,
                                const std::vector<int> &forced_ids,
                                uuber::BitsetState &out_bits, bool &out_valid);
bool ensure_forced_bitset_capacity(const uuber::NativeContext &ctx,
                                   uuber::BitsetState &bits_out,
                                   bool &bits_valid);
void set_forced_id_bit_strict(const uuber::NativeContext &ctx, int id,
                              uuber::BitsetState &bits_out,
                              bool &bits_valid);
bool forced_bits_any(const uuber::BitsetState *bits, bool bits_valid);
bool forced_bits_contains_label_id_strict(const uuber::NativeContext &ctx,
                                          int label_id,
                                          const uuber::BitsetState &bits,
                                          bool bits_valid);

struct ForcedScopeFilter {
  const int *source_ids_data{nullptr};
  int source_ids_count{0};
  const std::uint64_t *source_mask_words{nullptr};
  int source_mask_count{0};
  const std::unordered_map<int, int> *label_id_to_bit_idx{nullptr};
  const ForcedScopeFilter *parent{nullptr};
};

bool scope_filter_allows_id(const ForcedScopeFilter *scope_filter, int id);
bool forced_bits_contains_scoped(
    int id, const ForcedScopeFilter *scope_filter,
    const uuber::BitsetState *forced_bits, bool forced_bits_valid,
    const std::unordered_map<int, int> *label_id_to_bit_idx);
bool forced_bits_intersects_scope(
    const ForcedScopeFilter *scope_filter,
    const uuber::BitsetState *forced_bits, bool forced_bits_valid,
    const std::unordered_map<int, int> *label_id_to_bit_idx);

struct ForcedStateView {
  const ForcedScopeFilter *scope_filter{nullptr};
  const uuber::BitsetState *forced_complete_bits{nullptr};
  const bool *forced_complete_bits_valid{nullptr};
  const uuber::BitsetState *forced_survive_bits{nullptr};
  const bool *forced_survive_bits_valid{nullptr};
  const std::unordered_map<int, int> *label_id_to_bit_idx{nullptr};
};

ForcedStateView make_forced_state_view(
    const ForcedScopeFilter *scope_filter,
    const uuber::BitsetState *forced_complete_bits,
    const bool *forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    const bool *forced_survive_bits_valid,
    const std::unordered_map<int, int> *label_id_to_bit_idx);
bool forced_state_complete_valid(const ForcedStateView &view);
bool forced_state_survive_valid(const ForcedStateView &view);
bool forced_state_contains_complete(const ForcedStateView &view, int id);
bool forced_state_contains_survive(const ForcedStateView &view, int id);
bool forced_state_intersects_scope_complete(const ForcedStateView &view);
bool forced_state_intersects_scope_survive(const ForcedStateView &view);
