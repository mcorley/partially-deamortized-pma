// pma.cc
// Packed-Memory Array
// Written by Michael Corley
// Stony Brook University
// Cse638 Spring 2012

#include <stdint.h>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <vector>
#include "pma.h"

using namespace std;

pma::pma() 
  : _implicit_tree_height(log2(INITIAL_CAPACITY)),    
    _segment_size(INITIAL_CAPACITY / _implicit_tree_height),
    _size(0)
{
  _free_index_bitmap.resize(INITIAL_CAPACITY);
  _storage.resize(INITIAL_CAPACITY);  
}

pma::~pma() {
}

int& pma::operator[] (uint32_t& n) {
  return _storage[n];
}

const int& pma::operator[] (uint32_t& n) const {
  return _storage[n];
}

void pma::insert(const int& x)
{
  const uint32_t segment = segment_to_insert(x);
  const uint32_t pos = position_to_insert(segment, x);
  if(index_is_free(pos)) {
    _storage[pos] = x;
    _free_index_bitmap[pos].flip();
    _size++;
    return;
  }

  // TODO:
  // Rearrange the elements within the leaf (segment) 
  // to make room for x.
  // 1a.) Find first free index y after pos
  // 1a.) Find first free index y before pos
  // 2.) Compare distances and shift elements right or left
  // 3.) insert x into pos

  // If segment density exceeds its upper density threshold from
  // inserting x, start the rebalance algorithm.
  const double density = window_size(segment, _segment_size) / _segment_size;
  if (upper_density_threshold(0) <= density)
    rebalance(segment);
}

uint32_t pma::segment_to_insert(const int& x) const
{
  // O(logn) steps to narrow down a segment to scan.
  for (uint32_t seg = 0; seg < capacity(); seg += _segment_size) {
    for (uint32_t i = seg; seg < _segment_size; ++i) {
      if (index_is_free(i)) 
        continue;
      if (x  > _storage[i]) 
        break;
      if (x == _storage[i]) 
        return seg;
      else // (x  < _storage[i]) 
        return seg - _segment_size;
    }
  }
}

uint32_t pma::position_to_insert(const uint32_t& segment, const int& x) const
{
  // O(logn) steps to locate the index position to insert x.
  uint32_t pos = 0;
  for (uint32_t i = segment; i < _segment_size - 1; ++i) {
    if (index_is_free(i)) 
      continue;
    if (x  > _storage[i])
      pos = i + 1;
    if (x == _storage[i])
      pos = i + 1;
    else // (x  < _storage[i]) 
      return pos;
  }
  return pos;
}

void pma::rebalance(const uint32_t& segment)
{
  uint32_t window = (segment % number_of_segments() == 0) ? 
    segment : segment - _segment_size;
  uint32_t length = _segment_size * 2;

  // Find closest ancenstor whose density is within the permitted 
  // threshold.
  for (int height = 1; height < _implicit_tree_height; ++height) {
    const double udt = upper_density_threshold(height);
    const uint32_t cap = window_capacity(height);
    const uint32_t sz = window_size(window, length);

    if (sz/cap < udt) // This ancestor is within permitted threshold
      break;

    // This ancestor is also out of balance!
    if (window == 0 && window + length == capacity() - 1) {
      resize();
      return;
    }

    if (window == 0) {
      length += cap;
    } else if (window + length == capacity() - 1) {
      window -= cap;
    } else {
      window -= cap / 2;
      length += cap / 2;
    }
  }
  // Perform the naive two phase rebalance algorithm.
  naive_rebalance(window, length);

  // TODO:
  // one_phase_rebalance(window, length);
}

void pma::clear_window(const uint32_t& window, const uint32_t& length)
{
  for (int i = window; i < window + length; ++i) {
    _storage[i] = 0;
    _free_index_bitmap[i] = false;
  }
}

void pma::naive_rebalance(const uint32_t& window, const uint32_t& length)
{
  // We rebalance a node as follows:
  //    1.) Compress all the elements to the left part of the node without 
  //        adding empty spaces.
  //    2.) Evenly space out those elements, proceeding from left to right.
  // This rebalance algorithm requires two phases, and each phase needs to scan
  // the whole node.
  const int size = window_size(window, length);
  int next_index = 0;
  for (int i = window; i < window + length; ++i) {
    if (index_is_free(i))
      continue;
    if (next_index == i)
      next_index++;
    else { 
      _storage[next_index] = _storage[i];
      _free_index_bitmap[next_index].flip();
      _storage[i] = 0;
      _free_index_bitmap[i].flip();
      next_index++;
    }
  }
  next_index = size - 1;
  const int gap = length / size;
  for (int i = window + length - 1; i >= window; i -= gap) {
    _storage[i] = _storage[next_index];
    _free_index_bitmap[i].flip();
    _storage[next_index] = 0;
    _free_index_bitmap[next_index].flip();
    next_index--;
  }
}

void pma::resize()
{
  const int new_capacity = capacity() << 1;
  _implicit_tree_height = log2(new_capacity);
  _segment_size = new_capacity / _implicit_tree_height;
  _free_index_bitmap.resize(new_capacity);
  _storage.resize(new_capacity);
  naive_rebalance(0, new_capacity);
}

uint32_t pma::capacity() const {
  return _storage.size();
}

uint32_t pma::size() const {
  return _size;
}

uint32_t pma::segment_size() const {
  return _segment_size;
}

int pma::tree_height() const {
  return _implicit_tree_height;
}

int pma::number_of_segments() const {
  return _storage.size() / _segment_size;
}

bool pma::index_is_free(uint32_t indexno) const {
  return _free_index_bitmap[indexno] == false;
}

double pma::upper_density_threshold(int height) const
{
  return ROOT_UPPER_DENSITY + (LEAF_UPPER_DENSITY - ROOT_UPPER_DENSITY) * 
        (_implicit_tree_height - height) / _implicit_tree_height;
}

double pma::lower_density_threshold(int height) const
{
  return ROOT_LOWER_DENSITY - (ROOT_LOWER_DENSITY - LEAF_LOWER_DENSITY) * 
        (_implicit_tree_height - height) / _implicit_tree_height;
}

uint32_t pma::window_capacity(int height) const {
  return _segment_size << height;
}

uint32_t pma::window_size(const uint32_t& window, const uint32_t length) const
{
  uint32_t sz = 0;
  for (int i = window; i < window + length; ++i) {
    if (!index_is_free(i)) 
      sz++;
  }
  return sz;
}

