# Bloom Filter

[![Build Status](https://travis-ci.com/sagalasan/bloom-filter.svg?branch=master)](https://travis-ci.com/sagalasan/bloom-filter)

A Bloom Filter implementation in Rust.

## Installation

```toml
[dependencies]
bloom-filter = "0.1"
```

## Usage

```rust,no_run
extern crate bloom_filter;

use std::vec::Vec;
use bloom_filter::{BloomFilter, Murmur3};

let words = vec!["Hello", "I", "am", "some", "words"];

let mut bloom_filter = BloomFilter::optimal(Murmur3, words.len() as u64, 0.01);

bloom_filter.insert_all(&words);

for word in words.iter() {
    assert!(bloom_filter.contains(&word));
}
```