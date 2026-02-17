#include "context.h"

#include "prep_builder.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_set>

namespace {

using namespace uuber;

std::vector<std::string> extract_string_vector(SEXP obj) {
  if (Rf_isNull(obj))
    return {};
  Rcpp::CharacterVector chr(obj);
  std::vector<std::string> out;
  out.reserve(chr.size());
  for (R_xlen_t i = 0; i < chr.size(); ++i) {
    if (chr[i] == NA_STRING)
      continue;
    out.push_back(Rcpp::as<std::string>(chr[i]));
  }
  return out;
}

void populate_component_metadata(const Rcpp::List &prep, NativeContext &ctx) {
  Rcpp::List components = prep["components"];
  ctx.component_index.clear();
  ctx.component_info.clear();
  if (components.isNULL())
    return;
  Rcpp::CharacterVector comp_ids = components["ids"];
  Rcpp::List comp_attrs = components["attrs"];
  for (R_xlen_t i = 0; i < comp_ids.size(); ++i) {
    if (comp_ids[i] == NA_STRING)
      continue;
    std::string comp_label = Rcpp::as<std::string>(comp_ids[i]);
    ComponentContextInfo meta;
    Rcpp::RObject attrs_obj = comp_attrs[i];
    if (!attrs_obj.isNULL()) {
      Rcpp::List attrs(attrs_obj);
      if (attrs.containsElementNamed("guess")) {
        Rcpp::RObject guess_obj = attrs["guess"];
        if (!guess_obj.isNULL()) {
          Rcpp::List guess_list(guess_obj);
          if (guess_list.containsElementNamed("outcome")) {
            meta.guess.target = Rcpp::as<std::string>(guess_list["outcome"]);
          }
          if (meta.guess.target.empty()) {
            meta.guess.target = "__GUESS__";
          }
          if (!meta.guess.target.empty() && meta.guess.target != "__GUESS__") {
            auto it = ctx.outcome_index.find(meta.guess.target);
            if (it != ctx.outcome_index.end() && !it->second.empty()) {
              // For guess targets, if ambiguous, pick the first one?
              // Or should guesses follow component logic?
              // For now, default to first defined outcome.
              meta.guess.target_outcome_idx = it->second[0];
            }
          }
          if (guess_list.containsElementNamed("weights")) {
            Rcpp::NumericVector weights(guess_list["weights"]);
            Rcpp::CharacterVector weight_names = weights.names();
            if (!weight_names.isNULL()) {
              for (R_xlen_t j = 0; j < weights.size(); ++j) {
                if (weight_names[j] == NA_STRING)
                  continue;
                std::string source_label =
                    Rcpp::as<std::string>(weight_names[j]);
                double keep_prob = weights[j];
                if (source_label == "<NA>" || source_label == "NA") {
                  meta.guess.keep_weight_na = keep_prob;
                  meta.guess.has_keep_weight_na = true;
                  continue;
                }
                auto it = ctx.outcome_index.find(source_label);
                if (it == ctx.outcome_index.end() || it->second.empty())
                  continue;
                // Add weighted entry for potentially multiple outcomes of same
                // label? Probably yes. If "Left" is split, we guess "Left"
                // (either one).
                for (int idx : it->second) {
                  meta.guess.keep_weights.emplace_back(idx, keep_prob);
                }
              }
            }
          }
          meta.guess.valid =
              !meta.guess.keep_weights.empty() || meta.guess.has_keep_weight_na;
        }
      }
    }
    int idx = static_cast<int>(ctx.component_info.size());
    ctx.component_index[comp_label] = idx;
    ctx.component_info.push_back(std::move(meta));
  }
}

void populate_outcome_metadata(const Rcpp::List &prep, NativeContext &ctx) {
  Rcpp::List outcomes = prep["outcomes"];
  ctx.outcome_index.clear();
  ctx.outcome_labels.clear();
  ctx.outcome_label_ids.clear();
  ctx.outcome_info.clear();

  auto add_outcome_entry = [&](const std::string &label) -> int {
    int idx = static_cast<int>(ctx.outcome_info.size());
    ctx.outcome_index[label].push_back(idx);
    ctx.outcome_labels.push_back(label);
    auto it = ctx.label_to_id.find(label);
    ctx.outcome_label_ids.push_back((it == ctx.label_to_id.end()) ? NA_INTEGER
                                                                   : it->second);
    ctx.outcome_info.push_back(OutcomeContextInfo());
    return idx;
  };

  if (!outcomes.isNULL()) {
    Rcpp::CharacterVector outcome_names = outcomes.names();
    // Use loop to add all outcomes sequentially, allowing duplicates
    for (R_xlen_t i = 0; i < outcomes.size(); ++i) {
      if (outcome_names[i] == NA_STRING)
        continue;
      std::string outcome_label = Rcpp::as<std::string>(outcome_names[i]);

      // Always create a NEW outcome entry for each definition
      // We process definitions in order.
      Rcpp::RObject def_obj = outcomes[i];
      if (def_obj.isNULL())
        continue; // Skip NULLs but maybe shouldn't happen?

      int idx = add_outcome_entry(outcome_label);
      OutcomeContextInfo &info =
          ctx.outcome_info[static_cast<std::size_t>(idx)];
      info.node_id = -1;

      Rcpp::List def(def_obj);
      if (def.containsElementNamed("expr")) {
        Rcpp::RObject expr = def["expr"];
        if (!expr.isNULL()) {
          Rcpp::RObject lik_id = expr.attr(".lik_id");
          if (!lik_id.isNULL()) {
            info.node_id = Rcpp::as<int>(lik_id);
          }
        }
      }
      Rcpp::List options = def.containsElementNamed("options")
                               ? Rcpp::List(def["options"])
                               : Rcpp::List();
      if (!options.isNULL()) {
        if (options.containsElementNamed("component")) {
          Rcpp::RObject comp_obj = options["component"];
          if (!comp_obj.isNULL()) {
            info.allowed_components = extract_string_vector(comp_obj);
          }
        }
        if (options.containsElementNamed("map_outcome_to")) {
          Rcpp::RObject map_obj = options["map_outcome_to"];
          if (!map_obj.isNULL()) {
            Rcpp::CharacterVector map_cv(map_obj);
            if (map_cv.size() > 0) {
              Rcpp::String map_val = map_cv[0];
              if (map_val == NA_STRING) {
                info.maps_to_na = true;
              } else {
                std::string target = Rcpp::as<std::string>(
                    Rcpp::CharacterVector::create(map_val)[0]);

                auto it = ctx.outcome_index.find(target);
                if (it != ctx.outcome_index.end()) {
                  for (int t_idx : it->second) {
                    ctx.outcome_info[static_cast<std::size_t>(t_idx)]
                        .alias_sources.push_back(idx);
                  }
                }
              }
            }
          }
        }
        if (options.containsElementNamed("guess")) {
          Rcpp::RObject guess_obj = options["guess"];
          if (!guess_obj.isNULL()) {
            Rcpp::List guess_list(guess_obj);
            Rcpp::CharacterVector labels =
                guess_list.containsElementNamed("labels")
                    ? Rcpp::CharacterVector(guess_list["labels"])
                    : Rcpp::CharacterVector();
            Rcpp::NumericVector weights =
                guess_list.containsElementNamed("weights")
                    ? Rcpp::NumericVector(guess_list["weights"])
                    : Rcpp::NumericVector();
            std::string rt_policy = "keep";
            if (guess_list.containsElementNamed("rt_policy")) {
              Rcpp::RObject pol = guess_list["rt_policy"];
              if (!pol.isNULL()) {
                rt_policy = Rcpp::as<std::string>(pol);
              }
            }
            if (labels.size() == weights.size() && labels.size() > 0) {
              for (R_xlen_t j = 0; j < labels.size(); ++j) {
                if (labels[j] == NA_STRING)
                  continue;
                std::string target = Rcpp::as<std::string>(labels[j]);
                auto it = ctx.outcome_index.find(target);
                if (it != ctx.outcome_index.end()) {
                  for (int target_idx : it->second) {
                    OutcomeGuessDonor donor;
                    donor.outcome_idx = idx;
                    donor.weight = weights[j];
                    donor.rt_policy = rt_policy;
                    ctx.outcome_info[static_cast<std::size_t>(target_idx)]
                        .guess_donors.push_back(std::move(donor));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  Rcpp::List competitors = prep[".competitors"];
  if (!competitors.isNULL()) {
    bool by_index = false;
    Rcpp::RObject by_index_attr = competitors.attr("by_index");
    if (!by_index_attr.isNULL()) {
      by_index = Rcpp::as<bool>(by_index_attr);
    }
    if (by_index && competitors.size() == outcomes.size()) {
      for (R_xlen_t i = 0; i < competitors.size(); ++i) {
        Rcpp::RObject comp_obj = competitors[i];
        if (comp_obj.isNULL())
          continue;
        Rcpp::List comp_list(comp_obj);
        std::vector<int> comp_ids;
        comp_ids.reserve(comp_list.size());
        for (R_xlen_t j = 0; j < comp_list.size(); ++j) {
          Rcpp::RObject entry = comp_list[j];
          if (entry.isNULL())
            continue;
          Rcpp::RObject lik_id = entry.attr(".lik_id");
          if (lik_id.isNULL())
            continue;
          comp_ids.push_back(Rcpp::as<int>(lik_id));
        }
        if (i < static_cast<R_xlen_t>(ctx.outcome_info.size())) {
          ctx.outcome_info[static_cast<std::size_t>(i)].competitor_ids =
              comp_ids;
        }
      }
    } else {
      Rcpp::CharacterVector comp_names = competitors.names();
      for (R_xlen_t i = 0; i < competitors.size(); ++i) {
        if (comp_names[i] == NA_STRING)
          continue;
        Rcpp::RObject comp_obj = competitors[i];
        if (comp_obj.isNULL())
          continue;
        std::string outcome_label = Rcpp::as<std::string>(comp_names[i]);
        auto it = ctx.outcome_index.find(outcome_label);
        if (it == ctx.outcome_index.end())
          continue;

        Rcpp::List comp_list(comp_obj);
        std::vector<int> comp_ids;
        comp_ids.reserve(comp_list.size());
        for (R_xlen_t j = 0; j < comp_list.size(); ++j) {
          Rcpp::RObject entry = comp_list[j];
          if (entry.isNULL())
            continue;
          Rcpp::RObject lik_id = entry.attr(".lik_id");
          if (lik_id.isNULL())
            continue;
          comp_ids.push_back(Rcpp::as<int>(lik_id));
        }
        for (int idx : it->second) {
          ctx.outcome_info[static_cast<std::size_t>(idx)].competitor_ids =
              comp_ids;
        }
      }
    }
  }
}

void assign_component_indices(NativeContext &ctx) {
  for (auto &acc : ctx.accumulators) {
    acc.component_indices.clear();
    acc.component_indices.reserve(acc.components.size());
    for (const auto &comp : acc.components) {
      auto it = ctx.component_index.find(comp);
      if (it != ctx.component_index.end()) {
        acc.component_indices.push_back(it->second);
      }
    }
  }
}

void assign_label_outcome_indices(NativeContext &ctx) {
  for (auto &node : ctx.nodes) {
    if (node.source.empty())
      continue;
    auto it = ctx.outcome_index.find(node.source);
    if (it != ctx.outcome_index.end() && !it->second.empty()) {
      node.source_ref.outcome_idx = it->second[0];
    }
  }
  for (auto &pool : ctx.pools) {
    if (pool.members.empty() || pool.member_refs.empty())
      continue;
    std::size_t n = pool.members.size();
    if (pool.member_refs.size() != n) {
      pool.member_refs.resize(n);
    }
    for (std::size_t i = 0; i < n; ++i) {
      auto it = ctx.outcome_index.find(pool.members[i]);
      if (it != ctx.outcome_index.end() && !it->second.empty()) {
        pool.member_refs[i].outcome_idx = it->second[0];
      }
    }
  }
}

void populate_label_id_vectors(NativeContext &ctx) {
  ctx.accumulator_label_ids.assign(ctx.accumulators.size(), NA_INTEGER);
  for (std::size_t i = 0; i < ctx.accumulators.size(); ++i) {
    const std::string &label = ctx.accumulators[i].id;
    auto it = ctx.label_to_id.find(label);
    if (it != ctx.label_to_id.end()) {
      ctx.accumulator_label_ids[i] = it->second;
    }
  }

  ctx.pool_label_ids.assign(ctx.pools.size(), NA_INTEGER);
  for (std::size_t i = 0; i < ctx.pools.size(); ++i) {
    const std::string &label = ctx.pools[i].id;
    auto it = ctx.label_to_id.find(label);
    if (it != ctx.label_to_id.end()) {
      ctx.pool_label_ids[i] = it->second;
    }
  }

  if (ctx.outcome_label_ids.size() != ctx.outcome_labels.size()) {
    ctx.outcome_label_ids.assign(ctx.outcome_labels.size(), NA_INTEGER);
  }
  for (std::size_t i = 0; i < ctx.outcome_labels.size(); ++i) {
    if (ctx.outcome_label_ids[i] != NA_INTEGER) {
      continue;
    }
    auto it = ctx.label_to_id.find(ctx.outcome_labels[i]);
    if (it != ctx.label_to_id.end()) {
      ctx.outcome_label_ids[i] = it->second;
    }
  }
}

ComponentMap build_component_map(const Rcpp::List &prep,
                                 const NativeContext &ctx) {
  ComponentMap map;
  Rcpp::List components = prep["components"];
  if (!components.isNULL()) {
    Rcpp::CharacterVector comp_ids =
        components.containsElementNamed("ids")
            ? Rcpp::CharacterVector(components["ids"])
            : Rcpp::CharacterVector();
    Rcpp::NumericVector weights =
        components.containsElementNamed("weights")
            ? Rcpp::NumericVector(components["weights"])
            : Rcpp::NumericVector();
    map.ids.reserve(comp_ids.size());
    map.base_weights.reserve(comp_ids.size());
    for (R_xlen_t i = 0; i < comp_ids.size(); ++i) {
      if (comp_ids[i] == NA_STRING)
        continue;
      map.ids.push_back(Rcpp::as<std::string>(comp_ids[i]));
      double w = (i < weights.size()) ? static_cast<double>(weights[i]) : 1.0;
      if (!std::isfinite(w) || w < 0.0)
        w = 1.0;
      map.base_weights.push_back(w);
    }
  }
  if (map.ids.empty()) {
    map.ids.push_back("__default__");
    map.base_weights.push_back(1.0);
  }
  map.leader_idx.assign(map.ids.size(), -1);
  for (std::size_t c = 0; c < map.ids.size(); ++c) {
    const std::string &cid = map.ids[c];
    for (std::size_t a = 0; a < ctx.accumulators.size(); ++a) {
      const auto &comps = ctx.accumulators[a].components;
      if (std::find(comps.begin(), comps.end(), cid) != comps.end()) {
        map.leader_idx[c] = static_cast<int>(a);
        break;
      }
    }
    if (map.leader_idx[c] < 0 && !ctx.accumulators.empty()) {
      map.leader_idx[c] = 0;
    }
  }
  return map;
}

int resolve_na_cache_limit() {
  int default_limit = 128;
  Rcpp::Function getOption("getOption");
  Rcpp::RObject opt =
      getOption("uuber.cache.na.max_per_trial", Rcpp::wrap(default_limit));
  if (opt.isNULL())
    return default_limit;
  try {
    double val = Rcpp::as<double>(opt);
    if (!std::isfinite(val) || val < 0.0)
      return 0;
    return static_cast<int>(std::floor(val));
  } catch (...) {
    return default_limit;
  }
}

inline int rt_policy_code(const std::string &rt_policy) {
  if (rt_policy == "keep")
    return 0;
  if (rt_policy == "na")
    return 1;
  if (rt_policy == "drop")
    return 2;
  return 3;
}

inline int mask_word_count(int bit_count) {
  if (bit_count <= 0)
    return 0;
  return (bit_count + 63) / 64;
}

inline void set_mask_bit(std::vector<std::uint64_t> &mask, int bit_idx) {
  if (bit_idx < 0)
    return;
  int word = bit_idx / 64;
  int bit = bit_idx % 64;
  if (word < 0 || word >= static_cast<int>(mask.size()))
    return;
  mask[static_cast<std::size_t>(word)] |= (1ULL << bit);
}

inline int append_mask_words(std::vector<std::uint64_t> &flat,
                             const std::vector<std::uint64_t> &mask) {
  if (mask.empty())
    return -1;
  int begin = static_cast<int>(flat.size());
  flat.insert(flat.end(), mask.begin(), mask.end());
  return begin;
}

inline std::uint64_t mix_lookup_hash(std::uint64_t hash, int value) {
  std::uint32_t v = static_cast<std::uint32_t>(value);
  for (int i = 0; i < 4; ++i) {
    std::uint8_t b = static_cast<std::uint8_t>((v >> (8 * i)) & 0xFFU);
    hash ^= static_cast<std::uint64_t>(b);
    hash *= 1099511628211ULL;
  }
  return hash;
}

inline std::uint64_t shared_gate_lookup_key(int node_idx,
                                            const std::vector<int> &competitors) {
  std::uint64_t hash = 1469598103934665603ULL;
  hash = mix_lookup_hash(hash, node_idx);
  hash = mix_lookup_hash(hash, static_cast<int>(competitors.size()));
  for (int comp : competitors) {
    hash = mix_lookup_hash(hash, comp);
  }
  return hash;
}

bool mask_overlap(const IrContext &ir, const IrNode &a, const IrNode &b) {
  if (a.source_mask_count <= 0 || b.source_mask_count <= 0)
    return false;
  if (a.source_mask_count != b.source_mask_count)
    return false;
  for (int i = 0; i < a.source_mask_count; ++i) {
    std::uint64_t wa =
        ir.node_source_masks[static_cast<std::size_t>(a.source_mask_begin + i)];
    std::uint64_t wb =
        ir.node_source_masks[static_cast<std::size_t>(b.source_mask_begin + i)];
    if ((wa & wb) != 0ULL)
      return true;
  }
  return false;
}

bool and_event_pair(const IrContext &ir, int node_idx, int &event_a,
                    int &event_b) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ir.nodes.size()))
    return false;
  const IrNode &node = ir.nodes[static_cast<std::size_t>(node_idx)];
  if (node.op != IrNodeOp::And || node.child_count != 2)
    return false;
  int c0 = ir.node_children[static_cast<std::size_t>(node.child_begin)];
  int c1 = ir.node_children[static_cast<std::size_t>(node.child_begin + 1)];
  if (c0 < 0 || c1 < 0 || c0 >= static_cast<int>(ir.nodes.size()) ||
      c1 >= static_cast<int>(ir.nodes.size())) {
    return false;
  }
  const IrNode &n0 = ir.nodes[static_cast<std::size_t>(c0)];
  const IrNode &n1 = ir.nodes[static_cast<std::size_t>(c1)];
  if (n0.event_idx < 0 || n1.event_idx < 0)
    return false;
  event_a = n0.event_idx;
  event_b = n1.event_idx;
  return true;
}

bool detect_shared_gate_pair_ir(const IrContext &ir, int node_idx,
                                int competitor_idx,
                                IrSharedGateSpec &spec_out) {
  int ta = -1, tb = -1, ca = -1, cb = -1;
  if (!and_event_pair(ir, node_idx, ta, tb) ||
      !and_event_pair(ir, competitor_idx, ca, cb)) {
    return false;
  }
  if (ta < 0 || tb < 0 || ca < 0 || cb < 0)
    return false;
  if (ta >= static_cast<int>(ir.events.size()) ||
      tb >= static_cast<int>(ir.events.size()) ||
      ca >= static_cast<int>(ir.events.size()) ||
      cb >= static_cast<int>(ir.events.size())) {
    return false;
  }
  const int tl0 = ir.events[static_cast<std::size_t>(ta)].label_id;
  const int tl1 = ir.events[static_cast<std::size_t>(tb)].label_id;
  const int cl0 = ir.events[static_cast<std::size_t>(ca)].label_id;
  const int cl1 = ir.events[static_cast<std::size_t>(cb)].label_id;
  if (tl0 < 0 || tl1 < 0 || cl0 < 0 || cl1 < 0)
    return false;

  int pair_x = -1, pair_y = -1, pair_c = -1;
  if (tl0 == cl0) {
    pair_c = ta;
    pair_x = tb;
    pair_y = cb;
  } else if (tl0 == cl1) {
    pair_c = ta;
    pair_x = tb;
    pair_y = ca;
  } else if (tl1 == cl0) {
    pair_c = tb;
    pair_x = ta;
    pair_y = cb;
  } else if (tl1 == cl1) {
    pair_c = tb;
    pair_x = ta;
    pair_y = ca;
  } else {
    return false;
  }
  if (pair_x < 0 || pair_y < 0 || pair_c < 0)
    return false;
  const int x_label = ir.events[static_cast<std::size_t>(pair_x)].label_id;
  const int y_label = ir.events[static_cast<std::size_t>(pair_y)].label_id;
  const int c_label = ir.events[static_cast<std::size_t>(pair_c)].label_id;
  if (x_label < 0 || y_label < 0 || c_label < 0)
    return false;
  if (x_label == y_label || x_label == c_label || y_label == c_label)
    return false;

  spec_out = IrSharedGateSpec();
  spec_out.kind = IrSharedGateKind::Pair;
  spec_out.node_idx = node_idx;
  spec_out.pair_x_event_idx = pair_x;
  spec_out.pair_y_event_idx = pair_y;
  spec_out.pair_c_event_idx = pair_c;
  return true;
}

bool detect_shared_gate_nway_ir(const IrContext &ir, int node_idx,
                                const std::vector<int> &competitor_nodes,
                                std::vector<int> &nway_competitor_events,
                                IrSharedGateSpec &spec_out) {
  if (competitor_nodes.size() < 2)
    return false;
  std::vector<int> nodes;
  nodes.reserve(1 + competitor_nodes.size());
  nodes.push_back(node_idx);
  nodes.insert(nodes.end(), competitor_nodes.begin(), competitor_nodes.end());

  struct Pair {
    int e0{-1};
    int e1{-1};
    int l0{-1};
    int l1{-1};
  };
  std::vector<Pair> pairs(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    int e0 = -1, e1 = -1;
    if (!and_event_pair(ir, nodes[i], e0, e1))
      return false;
    if (e0 < 0 || e1 < 0 || e0 >= static_cast<int>(ir.events.size()) ||
        e1 >= static_cast<int>(ir.events.size())) {
      return false;
    }
    pairs[i] = Pair{e0, e1, ir.events[static_cast<std::size_t>(e0)].label_id,
                    ir.events[static_cast<std::size_t>(e1)].label_id};
    if (pairs[i].l0 < 0 || pairs[i].l1 < 0 || pairs[i].l0 == pairs[i].l1)
      return false;
  }

  const int cand1 = pairs[0].l0;
  const int cand2 = pairs[0].l1;
  bool cand1_all = true;
  bool cand2_all = true;
  for (std::size_t i = 1; i < pairs.size(); ++i) {
    const Pair &p = pairs[i];
    if (!(p.l0 == cand1 || p.l1 == cand1))
      cand1_all = false;
    if (!(p.l0 == cand2 || p.l1 == cand2))
      cand2_all = false;
  }
  if (cand1_all == cand2_all)
    return false;
  const int gate_label = cand1_all ? cand1 : cand2;

  int gate_event_idx = -1;
  int target_event_idx = -1;
  std::unordered_set<int> seen_x;
  seen_x.reserve(pairs.size());
  nway_competitor_events.clear();
  nway_competitor_events.reserve(competitor_nodes.size());

  for (std::size_t i = 0; i < pairs.size(); ++i) {
    const Pair &p = pairs[i];
    bool l0_is_gate = (p.l0 == gate_label);
    bool l1_is_gate = (p.l1 == gate_label);
    if (l0_is_gate == l1_is_gate)
      return false;
    int gate_e = l0_is_gate ? p.e0 : p.e1;
    int x_e = l0_is_gate ? p.e1 : p.e0;
    int x_label = l0_is_gate ? p.l1 : p.l0;
    if (x_label == gate_label || x_label < 0)
      return false;
    if (!seen_x.insert(x_label).second)
      return false;

    if (i == 0) {
      gate_event_idx = gate_e;
      target_event_idx = x_e;
    } else {
      nway_competitor_events.push_back(x_e);
    }
  }
  if (gate_event_idx < 0 || target_event_idx < 0 ||
      nway_competitor_events.empty()) {
    return false;
  }

  std::vector<int> event_node_indices;
  event_node_indices.reserve(2 + nway_competitor_events.size());
  event_node_indices.push_back(
      ir.events[static_cast<std::size_t>(gate_event_idx)].node_idx);
  event_node_indices.push_back(
      ir.events[static_cast<std::size_t>(target_event_idx)].node_idx);
  for (int ce : nway_competitor_events) {
    if (ce < 0 || ce >= static_cast<int>(ir.events.size()))
      return false;
    event_node_indices.push_back(ir.events[static_cast<std::size_t>(ce)].node_idx);
  }
  for (std::size_t i = 0; i < event_node_indices.size(); ++i) {
    int ni = event_node_indices[i];
    if (ni < 0 || ni >= static_cast<int>(ir.nodes.size()))
      continue;
    for (std::size_t j = i + 1; j < event_node_indices.size(); ++j) {
      int nj = event_node_indices[j];
      if (nj < 0 || nj >= static_cast<int>(ir.nodes.size()))
        continue;
      if (mask_overlap(ir, ir.nodes[static_cast<std::size_t>(ni)],
                       ir.nodes[static_cast<std::size_t>(nj)])) {
        return false;
      }
    }
  }

  spec_out = IrSharedGateSpec();
  spec_out.kind = IrSharedGateKind::NWay;
  spec_out.node_idx = node_idx;
  spec_out.nway_gate_event_idx = gate_event_idx;
  spec_out.nway_target_event_idx = target_event_idx;
  spec_out.nway_competitor_event_count =
      static_cast<int>(nway_competitor_events.size());
  return true;
}

} // namespace

namespace uuber {

void build_ir_context(NativeContext &ctx) {
  IrContext ir;
  ir.valid = false;
  ir.id_to_node_idx.clear();
  ir.id_to_node_idx.reserve(ctx.nodes.size());
  for (std::size_t i = 0; i < ctx.nodes.size(); ++i) {
    ir.id_to_node_idx[ctx.nodes[i].id] = static_cast<int>(i);
  }

  std::unordered_set<int> label_ids_set;
  label_ids_set.reserve(ctx.label_to_id.size() + ctx.nodes.size() * 2);
  for (const auto &kv : ctx.label_to_id) {
    if (kv.second >= 0)
      label_ids_set.insert(kv.second);
  }
  for (const auto &node : ctx.nodes) {
    if (node.source_ref.label_id >= 0) {
      label_ids_set.insert(node.source_ref.label_id);
    }
    for (int id : node.source_ids) {
      if (id >= 0)
        label_ids_set.insert(id);
    }
  }
  std::vector<int> label_ids(label_ids_set.begin(), label_ids_set.end());
  std::sort(label_ids.begin(), label_ids.end());
  ir.label_id_to_bit_idx.reserve(label_ids.size());
  for (std::size_t i = 0; i < label_ids.size(); ++i) {
    ir.label_id_to_bit_idx[label_ids[i]] = static_cast<int>(i);
  }
  ir.source_mask_words = mask_word_count(static_cast<int>(label_ids.size()));

  const int n_components = static_cast<int>(ctx.components.ids.size());
  ir.component_mask_words = mask_word_count(n_components);
  std::vector<std::uint64_t> all_component_mask(
      static_cast<std::size_t>(ir.component_mask_words), 0ULL);
  for (int i = 0; i < n_components; ++i) {
    set_mask_bit(all_component_mask, i);
  }

  ir.nodes.reserve(ctx.nodes.size());
  for (std::size_t i = 0; i < ctx.nodes.size(); ++i) {
    const NativeNode &node = ctx.nodes[i];
    IrNode ir_node;
    ir_node.node_id = node.id;
    ir_node.reference_idx = -1;
    ir_node.blocker_idx = -1;
    if (node.reference_id >= 0) {
      auto it = ir.id_to_node_idx.find(node.reference_id);
      if (it != ir.id_to_node_idx.end()) {
        ir_node.reference_idx = it->second;
      }
    }
    if (node.blocker_id >= 0) {
      auto it = ir.id_to_node_idx.find(node.blocker_id);
      if (it != ir.id_to_node_idx.end()) {
        ir_node.blocker_idx = it->second;
      }
    }
    if (node.needs_forced) {
      ir_node.flags |= IR_NODE_FLAG_NEEDS_FORCED;
    }
    if (node.scenario_sensitive) {
      ir_node.flags |= IR_NODE_FLAG_SCENARIO_SENSITIVE;
    }
    if (node.source == "__DEADLINE__") {
      ir_node.flags |= IR_NODE_FLAG_SPECIAL_DEADLINE;
    }
    if (node.source == "__GUESS__") {
      ir_node.flags |= IR_NODE_FLAG_SPECIAL_GUESS;
    }

    switch (node.kind_id) {
    case NODE_EVENT:
      if (node.source_ref.pool_idx >= 0) {
        ir_node.op = IrNodeOp::EventPool;
      } else {
        ir_node.op = IrNodeOp::EventAcc;
      }
      break;
    case NODE_AND:
      ir_node.op = IrNodeOp::And;
      break;
    case NODE_OR:
      ir_node.op = IrNodeOp::Or;
      break;
    case NODE_NOT:
      ir_node.op = IrNodeOp::Not;
      break;
    case NODE_GUARD:
      ir_node.op = IrNodeOp::Guard;
      break;
    default:
      ir_node.op = IrNodeOp::And;
      break;
    }

    std::vector<int> dense_children;
    if (node.kind_id == NODE_AND || node.kind_id == NODE_OR) {
      dense_children.reserve(node.args.size());
      for (int child_id : node.args) {
        auto it = ir.id_to_node_idx.find(child_id);
        if (it != ir.id_to_node_idx.end()) {
          dense_children.push_back(it->second);
        }
      }
    } else if (node.kind_id == NODE_NOT && node.arg_id >= 0) {
      auto it = ir.id_to_node_idx.find(node.arg_id);
      if (it != ir.id_to_node_idx.end()) {
        dense_children.push_back(it->second);
      }
    }
    ir_node.child_begin =
        dense_children.empty() ? -1 : static_cast<int>(ir.node_children.size());
    ir_node.child_count = static_cast<int>(dense_children.size());
    ir.node_children.insert(ir.node_children.end(), dense_children.begin(),
                            dense_children.end());

    std::vector<std::uint64_t> source_mask(
        static_cast<std::size_t>(ir.source_mask_words), 0ULL);
    if (ir.source_mask_words > 0) {
      if (node.source_ref.label_id >= 0) {
        auto it = ir.label_id_to_bit_idx.find(node.source_ref.label_id);
        if (it != ir.label_id_to_bit_idx.end()) {
          set_mask_bit(source_mask, it->second);
        }
      }
      for (int source_id : node.source_ids) {
        auto it = ir.label_id_to_bit_idx.find(source_id);
        if (it != ir.label_id_to_bit_idx.end()) {
          set_mask_bit(source_mask, it->second);
        }
      }
    }
    ir_node.source_mask_begin =
        append_mask_words(ir.node_source_masks, source_mask);
    ir_node.source_mask_count = ir.source_mask_words;

    if (node.kind_id == NODE_EVENT) {
      IrEvent event;
      event.node_idx = static_cast<int>(i);
      event.acc_idx = node.source_ref.acc_idx;
      event.pool_idx = node.source_ref.pool_idx;
      event.label_id = node.source_ref.label_id;
      event.outcome_idx = node.source_ref.outcome_idx;
      std::vector<std::uint64_t> component_mask(
          static_cast<std::size_t>(ir.component_mask_words), 0ULL);
      if (ir.component_mask_words > 0) {
        if (event.acc_idx >= 0 &&
            event.acc_idx < static_cast<int>(ctx.accumulators.size()) &&
            !ctx.accumulators[static_cast<std::size_t>(event.acc_idx)]
                 .component_indices.empty()) {
          const auto &idxs =
              ctx.accumulators[static_cast<std::size_t>(event.acc_idx)]
                  .component_indices;
          for (int comp_idx : idxs) {
            set_mask_bit(component_mask, comp_idx);
          }
        } else {
          component_mask = all_component_mask;
        }
      }
      event.component_mask_offset =
          append_mask_words(ir.component_masks, component_mask);
      ir_node.event_idx = static_cast<int>(ir.events.size());
      ir.events.push_back(std::move(event));
    }

    ir.nodes.push_back(std::move(ir_node));
  }

  ir.outcomes.reserve(ctx.outcome_info.size());
  for (std::size_t i = 0; i < ctx.outcome_info.size(); ++i) {
    const OutcomeContextInfo &info = ctx.outcome_info[i];
    IrOutcome ir_out;
    ir_out.label_id =
        (i < ctx.outcome_label_ids.size()) ? ctx.outcome_label_ids[i] : -1;
    if (ir_out.label_id == NA_INTEGER) {
      ir_out.label_id = -1;
    }
    auto node_it = ir.id_to_node_idx.find(info.node_id);
    ir_out.node_idx = (node_it == ir.id_to_node_idx.end()) ? -1 : node_it->second;

    ir_out.competitor_begin =
        info.competitor_ids.empty()
            ? -1
            : static_cast<int>(ir.outcome_competitors.size());
    for (int comp_id : info.competitor_ids) {
      auto it = ir.id_to_node_idx.find(comp_id);
      if (it != ir.id_to_node_idx.end()) {
        ir.outcome_competitors.push_back(it->second);
      }
    }
    ir_out.competitor_count =
        (ir_out.competitor_begin < 0)
            ? 0
            : static_cast<int>(ir.outcome_competitors.size()) -
                  ir_out.competitor_begin;

    std::vector<std::uint64_t> allowed_mask(
        static_cast<std::size_t>(ir.component_mask_words), 0ULL);
    if (ir.component_mask_words > 0) {
      if (info.allowed_components.empty()) {
        allowed_mask = all_component_mask;
      } else {
        for (const auto &comp : info.allowed_components) {
          auto it = ctx.component_index.find(comp);
          if (it != ctx.component_index.end()) {
            set_mask_bit(allowed_mask, it->second);
          }
        }
      }
    }
    ir_out.allowed_component_mask_offset =
        append_mask_words(ir.component_masks, allowed_mask);
    ir_out.maps_to_na = info.maps_to_na;

    ir_out.alias_begin = info.alias_sources.empty()
                             ? -1
                             : static_cast<int>(ir.outcome_alias_sources.size());
    ir.outcome_alias_sources.insert(ir.outcome_alias_sources.end(),
                                    info.alias_sources.begin(),
                                    info.alias_sources.end());
    ir_out.alias_count =
        (ir_out.alias_begin < 0)
            ? 0
            : static_cast<int>(ir.outcome_alias_sources.size()) -
                  ir_out.alias_begin;

    ir_out.guess_begin = info.guess_donors.empty()
                             ? -1
                             : static_cast<int>(ir.outcome_guess_donors.size());
    for (const auto &donor : info.guess_donors) {
      IrOutcomeGuessDonor ir_donor;
      ir_donor.outcome_idx = donor.outcome_idx;
      ir_donor.weight = donor.weight;
      ir_donor.rt_policy_code = rt_policy_code(donor.rt_policy);
      ir.outcome_guess_donors.push_back(std::move(ir_donor));
    }
    ir_out.guess_count =
        (ir_out.guess_begin < 0)
            ? 0
            : static_cast<int>(ir.outcome_guess_donors.size()) -
                  ir_out.guess_begin;
    ir.outcomes.push_back(std::move(ir_out));
    if (ir.outcomes.back().label_id >= 0) {
      ir.label_id_to_outcomes[ir.outcomes.back().label_id].push_back(
          static_cast<int>(i));
    }
    if (ir.outcomes.back().node_idx >= 0) {
      ir.node_idx_to_outcomes[ir.outcomes.back().node_idx].push_back(
          static_cast<int>(i));
    }
  }

  for (std::size_t oi = 0; oi < ir.outcomes.size(); ++oi) {
    IrOutcome &out = ir.outcomes[oi];
    if (out.node_idx < 0 || out.competitor_count <= 0)
      continue;
    std::vector<int> competitor_nodes;
    competitor_nodes.reserve(static_cast<std::size_t>(out.competitor_count));
    for (int k = 0; k < out.competitor_count; ++k) {
      competitor_nodes.push_back(
          ir.outcome_competitors[static_cast<std::size_t>(out.competitor_begin +
                                                          k)]);
    }
    IrSharedGateSpec spec;
    bool ok = false;
    if (competitor_nodes.size() == 1) {
      ok = detect_shared_gate_pair_ir(ir, out.node_idx, competitor_nodes[0], spec);
    } else if (competitor_nodes.size() >= 2) {
      std::vector<int> nway_events;
      ok = detect_shared_gate_nway_ir(ir, out.node_idx, competitor_nodes,
                                      nway_events, spec);
      if (ok && !nway_events.empty()) {
        spec.nway_competitor_event_begin =
            static_cast<int>(ir.nway_competitor_event_indices.size());
        ir.nway_competitor_event_indices.insert(
            ir.nway_competitor_event_indices.end(), nway_events.begin(),
            nway_events.end());
      }
    }
    if (ok) {
      spec.competitor_begin = out.competitor_begin;
      spec.competitor_count = out.competitor_count;
      int spec_idx = static_cast<int>(ir.shared_gate_specs.size());
      ir.shared_gate_specs.push_back(std::move(spec));
      out.shared_gate_spec_idx = spec_idx;
      std::vector<int> competitors;
      competitors.reserve(static_cast<std::size_t>(out.competitor_count));
      for (int k = 0; k < out.competitor_count; ++k) {
        competitors.push_back(
            ir.outcome_competitors[static_cast<std::size_t>(out.competitor_begin + k)]);
      }
      ir.shared_gate_lookup.emplace(
          shared_gate_lookup_key(out.node_idx, competitors), spec_idx);
    }
  }

  ir.valid = true;
  ctx.ir = std::move(ir);
}

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep) {
  NativePrepProto proto = build_prep_proto(prep);
  std::unique_ptr<NativeContext> ctx = build_context_from_proto(proto);
  populate_outcome_metadata(prep, *ctx);
  populate_component_metadata(prep, *ctx);
  ctx->components = build_component_map(prep, *ctx);
  assign_component_indices(*ctx);
  assign_label_outcome_indices(*ctx);
  populate_label_id_vectors(*ctx);
  build_ir_context(*ctx);
  ctx->na_cache_limit = resolve_na_cache_limit();
  return Rcpp::XPtr<NativeContext>(ctx.release(), true);
}

} // namespace uuber
