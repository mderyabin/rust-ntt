# rust-ntt

![Codecov](https://img.shields.io/codecov/c/github/mderyabin/rust-ntt)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/mderyabin/rust-ntt/ci.yml)

A high-performance Number Theoretic Transform (NTT) library for Rust,
implementing fast polynomial multiplication in finite fields.

The Number Theoretic Transform (NTT) is a discrete analogue of the Fast Fourier
Transform (FFT) that operates over finite fields instead of complex numbers.
It enables fast polynomial multiplication in `O(N log N)` time, significantly
improving over the naive `O(N^2)` method.

NTT is particularly useful in:

**Cryptography**: Post-quantum schemes like lattice-based cryptography
**Fully Homomorphic Encryption (FHE)**: Efficient polynomial operations
**Computer Algebra**: Fast polynomial arithmetic over finite fields.

This library implements negacyclic convolution in the ring `Z_q[x]/(x^N + 1)`,
where polynomials "wrap around" with a sign flip: `x^N = -1`.

[A Complete Beginner Guide to the Number Theoretic
Transform (NTT)](https://eprint.iacr.org/2024/585.pdf)

```rust
use rand::{SeedableRng, rngs::StdRng};
use rust_ntt::*;
use std::sync::Arc;

// Set up NTT context for degree-4 polynomials
const N: usize = 4;
let q = find_first_prime_up(10, N); // Find NTT-friendly prime
let ctx = NttContext::<N>::new(q);

// Create random polynomials
let mut rng = StdRng::seed_from_u64(42);
let a = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
let b = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

// Fast polynomial multiplication via NTT
let c = &a * &b;  // Uses negacyclic convolution

// Alternative explicit methods:
let c_ntt = a.negacyclic_convolution(&b);           // Barrett reduction
let c_shoup = a.negacyclic_convolution_shoup(&b);   // Shoup optimization

// All methods produce identical results
assert_eq!(c.coeffs(), c_ntt.coeffs());
assert_eq!(c.coeffs(), c_shoup.coeffs());
```

## Mathematical Background

This implementation is based on several key optimizations:


### Barrett Reduction

- [Barrett Reduction](https://en.wikipedia.org/wiki/Barrett_reduction)
- [Montgomery modular multiplication](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication)
  Precompute `μ = ⌊2^k/q⌋`` to replace division with multiplication and bit shifts

### Shoup Multiplication

[Shoup Multiplication](https://www.shoup.net/ntb/) - Further optimization for repeated multiplications with the same operand.

### Negacyclic NTT

[Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography](https://eprint.iacr.org/2016/504.pdf)
`x^N ≡ -1 (mod q)` creates the negacyclic structure needed for lattice cryptography.

- Prime modulus: `q` must be prime and satisfy `q ≡ 1 (mod 2N)``
- Polynomial degree: `N` must be a power of 2

Use the provided utilities to find suitable primes:

```
let q = find_first_prime_up(target_bits, N);  // Search upward
let q = find_first_prime_down(target_bits, N); // Search downward
```

---

## ✨ Features


- Multiple algorithms: Barrett reduction and Shoup optimization variants
- Const generics: Compile-time polynomial degrees for optimal performance
- Comprehensive testing: Property-based tests with proptest
- Educational examples: Clear demonstrations of NTT concepts

---

## 📦 Usage

### ✅ Run tests

```bash
cargo test
```

### ▶️ Run examples

```
cargo run --example convolution
cargo run --example primes
```
### ⚡Run benchmarks

```bash
cargo bench
```
