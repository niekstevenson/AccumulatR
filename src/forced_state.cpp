#include "forced_state.h"

#include <Rcpp.h>

#include <algorithm>

std::vector<int> forced_vec_from_sexp(SEXP vec) {
  std::vector<int> out;
  if (Rf_isNull(vec)) {
    return out;
  }
  Rcpp::IntegerVector iv(vec);
  out.reserve(iv.size());
  for (int val : iv) {
    if (val == NA_INTEGER) {
      continue;
    }
    out.push_back(val);
  }
  return out;
}

int forced_bit_capacity(const uuber::NativeContext &ctx) {
  if (ctx.ir.source_mask_words > 0) {
    return ctx.ir.source_mask_words * 64;
  }
  return static_cast<int>(ctx.ir.label_id_to_bit_idx.size());
}

void build_forced_bitset_strict(const uuber::NativeContext &ctx,
                                const std::vector<int> &forced_ids,
                                uuber::BitsetState &out_bits,
                                bool &out_valid) {
  out_valid = false;
  if (forced_ids.empty()) {
    return;
  }
  const int bit_capacity = forced_bit_capacity(ctx);
  if (bit_capacity <= 0) {
    Rcpp::stop("IR forced-state bitset unavailable for non-empty forced ids");
  }
  out_bits.reset_size(bit_capacity);
  out_valid = true;
  for (int id : forced_ids) {
    if (id == NA_INTEGER || id < 0) {
      continue;
    }
    auto it = ctx.ir.label_id_to_bit_idx.find(id);
    if (it == ctx.ir.label_id_to_bit_idx.end()) {
      Rcpp::stop("IR forced-state id %d not found in bit index map", id);
    }
    out_bits.set(it->second);
  }
}

bool ensure_forced_bitset_capacity(const uuber::NativeContext &ctx,
                                   uuber::BitsetState &bits_out,
                                   bool &bits_valid) {
  if (bits_valid) {
    return true;
  }
  const int bit_capacity = forced_bit_capacity(ctx);
  if (bit_capacity <= 0) {
    return false;
  }
  bits_out.reset_size(bit_capacity);
  bits_valid = true;
  return true;
}

void set_forced_id_bit_strict(const uuber::NativeContext &ctx, int id,
                              uuber::BitsetState &bits_out, bool &bits_valid) {
  if (id == NA_INTEGER || id < 0) {
    return;
  }
  auto bit_it = ctx.ir.label_id_to_bit_idx.find(id);
  if (bit_it == ctx.ir.label_id_to_bit_idx.end()) {
    Rcpp::stop("IR forced-state id %d not found in bit index map", id);
  }
  if (!ensure_forced_bitset_capacity(ctx, bits_out, bits_valid)) {
    Rcpp::stop("IR forced-state bitset unavailable for id %d", id);
  }
  bits_out.set(bit_it->second);
}

bool forced_bits_any(const uuber::BitsetState *bits, bool bits_valid) {
  return bits_valid && bits && bits->any();
}

bool forced_bits_contains_label_id_strict(const uuber::NativeContext &ctx,
                                          int label_id,
                                          const uuber::BitsetState &bits,
                                          bool bits_valid) {
  if (!bits_valid || label_id == NA_INTEGER || label_id < 0) {
    return false;
  }
  auto bit_it = ctx.ir.label_id_to_bit_idx.find(label_id);
  if (bit_it == ctx.ir.label_id_to_bit_idx.end()) {
    Rcpp::stop("IR forced-state id %d not found in bit index map", label_id);
  }
  return bits.test(bit_it->second);
}

bool scope_filter_allows_id(const ForcedScopeFilter *scope_filter, int id) {
  if (id == NA_INTEGER) {
    return false;
  }
  for (const ForcedScopeFilter *it = scope_filter; it != nullptr;
       it = it->parent) {
    if (it->source_mask_words && it->source_mask_count > 0 &&
        it->label_id_to_bit_idx) {
      auto bit_it = it->label_id_to_bit_idx->find(id);
      if (bit_it == it->label_id_to_bit_idx->end()) {
        return false;
      }
      const int bit_idx = bit_it->second;
      const int word_idx = bit_idx / 64;
      const int bit = bit_idx % 64;
      if (word_idx < 0 || word_idx >= it->source_mask_count) {
        return false;
      }
      const std::uint64_t word =
          it->source_mask_words[static_cast<std::size_t>(word_idx)];
      if ((word & (1ULL << bit)) == 0ULL) {
        return false;
      }
      continue;
    }
    if (it->source_ids_data && it->source_ids_count > 0) {
      if (!std::binary_search(it->source_ids_data,
                              it->source_ids_data + it->source_ids_count,
                              id)) {
        return false;
      }
    }
  }
  return true;
}

bool forced_bits_contains_scoped(
    int id, const ForcedScopeFilter *scope_filter,
    const uuber::BitsetState *forced_bits, bool forced_bits_valid,
    const std::unordered_map<int, int> *label_id_to_bit_idx) {
  if (id == NA_INTEGER || !forced_bits || !forced_bits_valid ||
      !label_id_to_bit_idx) {
    return false;
  }
  if (!scope_filter_allows_id(scope_filter, id)) {
    return false;
  }
  auto bit_it = label_id_to_bit_idx->find(id);
  if (bit_it == label_id_to_bit_idx->end()) {
    return false;
  }
  return forced_bits->test(bit_it->second);
}

bool forced_bits_intersects_scope(
    const ForcedScopeFilter *scope_filter,
    const uuber::BitsetState *forced_bits, bool forced_bits_valid,
    const std::unordered_map<int, int> *label_id_to_bit_idx) {
  if (!forced_bits || !forced_bits_valid || !forced_bits->any()) {
    return false;
  }
  if (!scope_filter) {
    return true;
  }
  if (scope_filter->parent == nullptr && scope_filter->source_mask_words &&
      scope_filter->source_mask_count > 0 && scope_filter->label_id_to_bit_idx &&
      label_id_to_bit_idx &&
      scope_filter->label_id_to_bit_idx == label_id_to_bit_idx) {
    return forced_bits->intersects_words(scope_filter->source_mask_words,
                                         scope_filter->source_mask_count);
  }
  if (!label_id_to_bit_idx) {
    return false;
  }
  for (const ForcedScopeFilter *it = scope_filter; it != nullptr;
       it = it->parent) {
    if (!it->source_ids_data || it->source_ids_count <= 0) {
      continue;
    }
    for (int i = 0; i < it->source_ids_count; ++i) {
      const int id = it->source_ids_data[i];
      auto bit_it = label_id_to_bit_idx->find(id);
      if (bit_it != label_id_to_bit_idx->end() &&
          forced_bits->test(bit_it->second)) {
        return true;
      }
    }
  }
  return false;
}

ForcedStateView make_forced_state_view(
    const ForcedScopeFilter *scope_filter,
    const uuber::BitsetState *forced_complete_bits,
    const bool *forced_complete_bits_valid,
    const uuber::BitsetState *forced_survive_bits,
    const bool *forced_survive_bits_valid,
    const std::unordered_map<int, int> *label_id_to_bit_idx) {
  ForcedStateView view;
  view.scope_filter = scope_filter;
  view.forced_complete_bits = forced_complete_bits;
  view.forced_complete_bits_valid = forced_complete_bits_valid;
  view.forced_survive_bits = forced_survive_bits;
  view.forced_survive_bits_valid = forced_survive_bits_valid;
  view.label_id_to_bit_idx = label_id_to_bit_idx;
  return view;
}

bool forced_state_complete_valid(const ForcedStateView &view) {
  return view.forced_complete_bits != nullptr &&
         view.forced_complete_bits_valid != nullptr &&
         *view.forced_complete_bits_valid;
}

bool forced_state_survive_valid(const ForcedStateView &view) {
  return view.forced_survive_bits != nullptr &&
         view.forced_survive_bits_valid != nullptr &&
         *view.forced_survive_bits_valid;
}

bool forced_state_contains_complete(const ForcedStateView &view, int id) {
  return forced_bits_contains_scoped(id, view.scope_filter,
                                     view.forced_complete_bits,
                                     forced_state_complete_valid(view),
                                     view.label_id_to_bit_idx);
}

bool forced_state_contains_survive(const ForcedStateView &view, int id) {
  return forced_bits_contains_scoped(id, view.scope_filter,
                                     view.forced_survive_bits,
                                     forced_state_survive_valid(view),
                                     view.label_id_to_bit_idx);
}

bool forced_state_intersects_scope_complete(const ForcedStateView &view) {
  return forced_bits_intersects_scope(view.scope_filter, view.forced_complete_bits,
                                      forced_state_complete_valid(view),
                                      view.label_id_to_bit_idx);
}

bool forced_state_intersects_scope_survive(const ForcedStateView &view) {
  return forced_bits_intersects_scope(view.scope_filter, view.forced_survive_bits,
                                      forced_state_survive_valid(view),
                                      view.label_id_to_bit_idx);
}
