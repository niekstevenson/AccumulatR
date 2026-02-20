#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace uuber {

class BitsetState {
public:
  BitsetState() = default;
  explicit BitsetState(int bit_count) { reset_size(bit_count); }

  void reset_size(int bit_count) {
    bit_count_ = bit_count < 0 ? 0 : bit_count;
    const int words = word_count(bit_count_);
    if (words <= 1) {
      words_.clear();
      small_word_ = 0ULL;
    } else {
      words_.assign(static_cast<std::size_t>(words), 0ULL);
      small_word_ = 0ULL;
    }
  }

  int bit_count() const { return bit_count_; }

  void clear_all() {
    small_word_ = 0ULL;
    std::fill(words_.begin(), words_.end(), 0ULL);
  }

  void set(int bit) {
    if (bit < 0 || bit >= bit_count_) {
      return;
    }
    const int word = bit / 64;
    const int off = bit % 64;
    if (words_.empty()) {
      small_word_ |= (1ULL << off);
    } else {
      words_[static_cast<std::size_t>(word)] |= (1ULL << off);
    }
  }

  void clear(int bit) {
    if (bit < 0 || bit >= bit_count_) {
      return;
    }
    const int word = bit / 64;
    const int off = bit % 64;
    if (words_.empty()) {
      small_word_ &= ~(1ULL << off);
    } else {
      words_[static_cast<std::size_t>(word)] &= ~(1ULL << off);
    }
  }

  bool test(int bit) const {
    if (bit < 0 || bit >= bit_count_) {
      return false;
    }
    const int word = bit / 64;
    const int off = bit % 64;
    if (words_.empty()) {
      return (small_word_ & (1ULL << off)) != 0ULL;
    }
    return (words_[static_cast<std::size_t>(word)] & (1ULL << off)) != 0ULL;
  }

  bool intersects_words(const std::uint64_t *mask_words, int mask_word_count) const {
    if (!mask_words || mask_word_count <= 0 || bit_count_ <= 0) {
      return false;
    }
    if (words_.empty()) {
      return (small_word_ & mask_words[0]) != 0ULL;
    }
    const int n = std::min<int>(mask_word_count, static_cast<int>(words_.size()));
    for (int i = 0; i < n; ++i) {
      if ((words_[static_cast<std::size_t>(i)] & mask_words[i]) != 0ULL) {
        return true;
      }
    }
    return false;
  }

  void and_not_words(const std::uint64_t *mask_words, int mask_word_count) {
    if (!mask_words || mask_word_count <= 0 || bit_count_ <= 0) {
      return;
    }
    if (words_.empty()) {
      small_word_ &= ~mask_words[0];
      return;
    }
    const int n = std::min<int>(mask_word_count, static_cast<int>(words_.size()));
    for (int i = 0; i < n; ++i) {
      words_[static_cast<std::size_t>(i)] &= ~mask_words[i];
    }
  }

  const std::vector<std::uint64_t> &words() const { return words_; }
  std::uint64_t small_word() const { return small_word_; }

private:
  static int word_count(int bit_count) { return (bit_count + 63) / 64; }

  int bit_count_{0};
  std::uint64_t small_word_{0ULL};
  std::vector<std::uint64_t> words_;
};

} // namespace uuber

