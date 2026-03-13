#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

struct SourceTimeConstraint {
  bool has_exact{false};
  double exact_time{std::numeric_limits<double>::quiet_NaN()};
  bool has_lower{false};
  double lower{0.0}; // exclusive lower bound
  bool has_upper{false};
  double upper{std::numeric_limits<double>::infinity()}; // inclusive upper bound
};

using TimeConstraintMap = std::map<int, SourceTimeConstraint>;

inline bool time_constraint_same_time(double lhs, double rhs) {
  if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
    return false;
  }
  const double scale =
      std::max(1.0, std::max(std::fabs(lhs), std::fabs(rhs)));
  return std::fabs(lhs - rhs) <= 1e-12 * scale;
}

inline bool time_constraint_has_bounds(const SourceTimeConstraint &constraint) {
  return constraint.has_lower || constraint.has_upper;
}

inline bool time_constraint_empty(const SourceTimeConstraint &constraint) {
  return !constraint.has_exact && !constraint.has_lower && !constraint.has_upper;
}

inline bool time_constraint_valid(const SourceTimeConstraint &constraint) {
  if (constraint.has_exact) {
    if (!std::isfinite(constraint.exact_time)) {
      return false;
    }
    if (constraint.has_lower &&
        (!(constraint.exact_time > constraint.lower) &&
         !time_constraint_same_time(constraint.exact_time, constraint.lower))) {
      return false;
    }
    if (constraint.has_upper && constraint.exact_time > constraint.upper &&
        !time_constraint_same_time(constraint.exact_time, constraint.upper)) {
      return false;
    }
    return true;
  }
  if (constraint.has_lower && constraint.has_upper &&
      (constraint.upper < constraint.lower ||
       time_constraint_same_time(constraint.upper, constraint.lower))) {
    return false;
  }
  return true;
}

inline const SourceTimeConstraint *
time_constraints_find(const TimeConstraintMap *time_constraints, int source_id) {
  if (time_constraints == nullptr) {
    return nullptr;
  }
  auto it = time_constraints->find(source_id);
  if (it == time_constraints->end()) {
    return nullptr;
  }
  return &it->second;
}

inline bool time_constraints_any(const TimeConstraintMap *time_constraints) {
  return time_constraints != nullptr && !time_constraints->empty();
}

inline SourceTimeConstraint &
time_constraints_ensure(int source_id, TimeConstraintMap &time_constraints) {
  auto it = time_constraints.find(source_id);
  if (it == time_constraints.end()) {
    it = time_constraints.emplace(source_id, SourceTimeConstraint{}).first;
  }
  return it->second;
}

inline void time_constraints_cleanup(int source_id,
                                     TimeConstraintMap &time_constraints) {
  auto it = time_constraints.find(source_id);
  if (it == time_constraints.end()) {
    return;
  }
  if (time_constraint_empty(it->second)) {
    time_constraints.erase(it);
  }
}

inline void time_constraints_remove(int source_id,
                                    TimeConstraintMap &time_constraints) {
  auto it = time_constraints.find(source_id);
  if (it != time_constraints.end()) {
    time_constraints.erase(it);
  }
}

inline bool time_constraints_bind_exact(int source_id, double t,
                                        TimeConstraintMap &time_constraints) {
  if (source_id < 0) {
    return true;
  }
  SourceTimeConstraint &constraint =
      time_constraints_ensure(source_id, time_constraints);
  if (constraint.has_exact) {
    return time_constraint_same_time(constraint.exact_time, t);
  }
  if (constraint.has_lower &&
      (!(t > constraint.lower) &&
       !time_constraint_same_time(t, constraint.lower))) {
    return false;
  }
  if (constraint.has_upper && t > constraint.upper &&
      !time_constraint_same_time(t, constraint.upper)) {
    return false;
  }
  constraint.has_exact = true;
  constraint.exact_time = t;
  constraint.has_lower = false;
  constraint.lower = 0.0;
  constraint.has_upper = false;
  constraint.upper = std::numeric_limits<double>::infinity();
  return true;
}

inline bool time_constraints_add_upper(int source_id, double upper,
                                       TimeConstraintMap &time_constraints) {
  if (source_id < 0) {
    return true;
  }
  SourceTimeConstraint &constraint =
      time_constraints_ensure(source_id, time_constraints);
  if (constraint.has_exact) {
    return constraint.exact_time <= upper ||
           time_constraint_same_time(constraint.exact_time, upper);
  }
  if (constraint.has_upper) {
    constraint.upper = std::min(constraint.upper, upper);
  } else {
    constraint.has_upper = true;
    constraint.upper = upper;
  }
  if (!time_constraint_valid(constraint)) {
    return false;
  }
  time_constraints_cleanup(source_id, time_constraints);
  return true;
}

inline bool time_constraints_add_lower(int source_id, double lower,
                                       TimeConstraintMap &time_constraints) {
  if (source_id < 0) {
    return true;
  }
  SourceTimeConstraint &constraint =
      time_constraints_ensure(source_id, time_constraints);
  if (constraint.has_exact) {
    return constraint.exact_time > lower &&
           !time_constraint_same_time(constraint.exact_time, lower);
  }
  if (constraint.has_lower) {
    constraint.lower = std::max(constraint.lower, lower);
  } else {
    constraint.has_lower = true;
    constraint.lower = lower;
  }
  if (!time_constraint_valid(constraint)) {
    return false;
  }
  time_constraints_cleanup(source_id, time_constraints);
  return true;
}

inline bool time_constraints_mark_complete(int source_id, double t,
                                           bool bind_exact_current_time,
                                           TimeConstraintMap &time_constraints) {
  if (bind_exact_current_time) {
    return time_constraints_bind_exact(source_id, t, time_constraints);
  }
  return time_constraints_add_upper(source_id, t, time_constraints);
}

inline bool time_constraints_mark_survive(int source_id, double t,
                                          TimeConstraintMap &time_constraints) {
  return time_constraints_add_lower(source_id, t, time_constraints);
}

inline bool time_constraint_guarantees_complete_at(
    const SourceTimeConstraint &constraint, double t) {
  if (constraint.has_exact) {
    return constraint.exact_time <= t ||
           time_constraint_same_time(constraint.exact_time, t);
  }
  if (constraint.has_upper) {
    return constraint.upper <= t || time_constraint_same_time(constraint.upper, t);
  }
  return false;
}

inline bool time_constraint_guarantees_survive_at(
    const SourceTimeConstraint &constraint, double t) {
  if (constraint.has_exact) {
    return constraint.exact_time > t &&
           !time_constraint_same_time(constraint.exact_time, t);
  }
  if (constraint.has_lower) {
    return constraint.lower > t || time_constraint_same_time(constraint.lower, t);
  }
  return false;
}

inline bool time_constraints_contains_complete_at(
    const TimeConstraintMap *time_constraints, int source_id, double t) {
  const SourceTimeConstraint *constraint =
      time_constraints_find(time_constraints, source_id);
  return constraint != nullptr &&
         time_constraint_guarantees_complete_at(*constraint, t);
}

inline bool time_constraints_contains_survive_at(
    const TimeConstraintMap *time_constraints, int source_id, double t) {
  const SourceTimeConstraint *constraint =
      time_constraints_find(time_constraints, source_id);
  return constraint != nullptr &&
         time_constraint_guarantees_survive_at(*constraint, t);
}
