use murmur3::murmur3_32;
use bit_vec::BitVec;

use std::io::{Cursor, Read};
use std::f64::consts::{LN_2, E};

/// A trait for hashing an arbitrary stream of bytes into a bloom filter.
pub trait BloomHasher {
    /// Returns the hashed value of the bytes given some seed.
    #[inline]
    fn hash(&self, seed: u32, bytes: &[u8]) -> u32;
}

/// A unit struct for the murmur3 hash function.
pub struct Murmur3;

impl BloomHasher for Murmur3 {
    fn hash(&self, seed: u32, bytes: &[u8]) -> u32 {
        let mut cursor = Cursor::new(bytes);
        murmur3_32(cursor.by_ref(), seed)
    }
}

/// BloomFilter
///
/// An implementation of a bloom filter
pub struct BloomFilter<T> {
    hasher: T,
    k: u32,
    bit_vec: BitVec,
    insert_count: u64,
}

impl<T: BloomHasher> BloomFilter<T> {
    /// Create a new `BloomFilter` given a `hasher`,
    /// the number of hash functions to use,
    /// and the size of the underlying bit array.
    ///
    /// Typically, this function should not be called directly unless,
    /// the optimal number of hash functions and optimal array size are
    /// already known.
    pub fn new(hasher: T, k: u32, array_size: u64) -> Self {
        Self {
            hasher,
            k,
            bit_vec: BitVec::from_elem(array_size as usize, false),
            insert_count: 0,
        }
    }

    /// Create a `BloomFilter` by computing its optimal parameters.
    ///
    /// This function computes the optimal array size using
    /// ```text
    /// -(n * ln(p)) / ln(2) ^ 2
    /// ```
    /// and computes the optimal number of hash functions using
    /// ```text
    /// m / n * ln(2)
    /// ```
    pub fn optimal(hasher: T, max_elements: u64, error_rate: f64) -> Self {
        // Check error rate is valid
        if error_rate <= 0_f64 || error_rate >= 1_f64 {
            panic!("Error rate must be 0 <= error_rate < 1");
        }

        // Calculate the length of the bit vector
        let m = optimal_vec_size(max_elements as u64, error_rate);

        // Calculate the number of hash functions to use
        let k = optimal_hash_functions(m, max_elements as u64);

        // Create the bloom filter
        Self::new(hasher, k, m)
    }

    /// Insert a slice of bytes into the `BloomFilter`.
    pub fn insert(&mut self, bytes: &[u8]) {
        for seed in 0..self.k {
            let hash = self.hasher.hash(seed, bytes) as usize % self.bit_vec.len();
            self.bit_vec.set(hash, true);
        }
        self.insert_count += 1;
    }

    /// Insert a slice of slices of bytes into the `BloomFilter`.
    pub fn insert_all<B: AsRef<[u8]>>(&mut self, slice: &[B]) {
        for item in slice {
            self.insert(item.as_ref());
        }
    }

    /// Check whether a slice of bytes exists in the `BloomFilter`.
    ///
    /// This is a probabilistic function that may return a false positive but will
    /// never return a false negative.
    ///
    /// # Examples
    ///
    /// ```
    /// extern crate bloom_filter;
    ///
    /// use std::vec::Vec;
    /// use bloom_filter::{BloomFilter, Murmur3};
    ///
    /// let words = vec!["Hello", "I", "am", "some", "words"];
    ///
    /// let mut bloom_filter = BloomFilter::optimal(Murmur3, words.len() as u64, 0.01);
    ///
    /// bloom_filter.insert_all(&words);
    ///
    /// for word in words.iter() {
    ///     assert!(bloom_filter.contains(&word));
    /// }
    /// ```
    pub fn contains<B: AsRef<[u8]>>(&self, bytes: B) -> bool {
        for seed in 0..self.k {
            let hash = self.hasher.hash(seed, bytes.as_ref()) as usize % self.bit_vec.len();
            if !self.bit_vec[hash] {
                return false;
            }
        }

        true
    }

    /// Calculate the expected false positive rate given the current state of
    /// the `BloomFilter`.
    pub fn false_positive_rate(&self) -> f64 {
        false_positive_rate(self.insert_count, self.bit_vec.len() as u64, self.k)
    }
}

/// This function computes the false positive rate given n, m, and k.
#[inline]
fn false_positive_rate(n: u64, m: u64, k: u32) -> f64 {
    (1_f64 - E.powf(-(n as f64) * m as f64 / k as f64)).powf(k as f64)
}

#[inline]
fn optimal_hash_functions(m: u64, n: u64) -> u32 {
    1_f64.max(m as f64 / n as f64 * LN_2).ceil() as u32
}

#[inline]
fn optimal_vec_size(n: u64, p: f64) -> u64 {
    (-(n as f64 * p.ln()) / LN_2.powi(2)).ceil() as u64
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufReader, BufRead};
    use super::*;

    #[test]
    fn test_optimal_hash_functions() {
        assert_eq!(1, optimal_hash_functions(1, 10));
        assert_eq!(7, optimal_hash_functions(95851, 10000));
        assert_eq!(7, optimal_hash_functions(9586, 1000));
        assert_eq!(7, optimal_hash_functions(90000, 10000));
    }

    #[test]
    fn test_optimal_vec_size() {
        assert_eq!(95851, optimal_vec_size(10000, 0.01));
        assert_eq!(9586, optimal_vec_size(1000, 0.01));
    }

    #[test]
    #[should_panic]
    fn test_error_rate_too_low() {
        BloomFilter::optimal(Murmur3, 10000, 0_f64);
    }

    #[test]
    #[should_panic]
    fn test_error_rate_too_high() {
        BloomFilter::optimal(Murmur3, 10000, 1_f64);
    }

    #[test]
    fn test_no_false_negatives() {
        let words: Vec<String> = BufReader::new(File::open("./resources/1000.txt").unwrap())
            .lines()
            .map(|s| s.unwrap())
            .collect();

        let mut bloom_filter = BloomFilter::optimal(Murmur3, words.len() as u64, 0.01);

        bloom_filter.insert_all(&words);

        for word in words.iter() {
            assert!(bloom_filter.contains(&word));
        }
//
//        let bloom_filter = BloomFilter::from_iter(words.iter(), Murmur3, 0.01);
//
//        for word in words.iter() {
//            assert!(bloom_filter.contains(&word));
//        }
    }
}