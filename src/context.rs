use crate::congruence::CongruenceClass;
use crate::math::find_generator;
use std::sync::Arc;

/// Shared NTT context containing precomputed values for a specific degree and modulus.
///
/// This context contains all the expensive-to-compute values needed for NTT operations:
/// - Modular arithmetic context (Barrett reduction parameters)
/// - Forward and inverse twiddle factors
/// - Shoup precomputed values for faster multiplication
/// - Normalization factor for inverse NTT
///
/// # Modulus Requirements
///
/// The modulus `q` must satisfy:
/// - `q` is prime
/// - `q ≡ 1 (mod 2*DEGREE)` (ensures primitive 2n-th root of unity exists)
/// - `q < 2^63` (fits in signed 64-bit for safe arithmetic)
///
/// Use utility functions like `find_first_prime_up(logq, DEGREE)` to find suitable moduli.
#[derive(Debug, Clone)]
pub struct NttContext<const DEGREE: usize> {
    /// Modular arithmetic context with Barrett reduction parameters
    pub(crate) class: CongruenceClass,
    /// Inverse of DEGREE modulo q, for NTT normalization
    pub(crate) inv_n: u64,
    /// Shoup precomputed value for inv_n
    pub(crate) inv_n_shoup: u64,
    /// Forward twiddle factors for NTT (bit-reversed order)
    pub(crate) tf: [u64; DEGREE],
    /// Shoup precomputed values for forward twiddle factors
    pub(crate) tf_shoup: [u64; DEGREE],
    /// Inverse twiddle factors for INTT (bit-reversed order)
    pub(crate) itf: [u64; DEGREE],
    /// Shoup precomputed values for inverse twiddle factors
    pub(crate) itf_shoup: [u64; DEGREE],
}

impl<const DEGREE: usize> NttContext<DEGREE> {
    pub fn tf(&self) -> &[u64; DEGREE] {
        &self.tf
    }

    pub fn itf(&self) -> &[u64; DEGREE] {
        &self.itf
    }

    pub fn class(&self) -> &CongruenceClass {
        &self.class
    }
}

impl<const DEGREE: usize> NttContext<DEGREE> {
    /// Create a new NTT context for the given modulus.
    ///
    /// # Arguments
    /// * `q` - Prime modulus satisfying q ≡ 1 (mod 2*DEGREE)
    ///
    /// # Panics
    /// * If DEGREE is not a power of 2
    /// * If q doesn't satisfy the modulus requirements
    ///
    /// # Examples
    /// ```
    /// use rust_ntt::{NttContext, find_first_prime_up};
    ///
    /// const N: usize = 1024;
    /// let q = find_first_prime_up(17, N); // Find suitable 17-bit prime
    /// let ctx = NttContext::<N>::new(q);
    /// ```
    pub fn new(q: u64) -> Arc<Self> {
        // Validate that DEGREE is a power of 2
        assert!(
            DEGREE.is_power_of_two() && DEGREE > 0,
            "DEGREE must be a power of 2, got {DEGREE}"
        );

        // Validate modulus requirements
        assert!(q >= 3, "Modulus must be at least 3, got {q}");
        assert!(q < (1u64 << 63), "Modulus must be < 2^63, got {q}");
        assert_eq!(
            (q - 1) % (2 * DEGREE as u64),
            0,
            "Modulus {q} must satisfy q ≡ 1 (mod 2*DEGREE={})",
            2 * DEGREE
        );

        let class = CongruenceClass::new(q);

        // Find generator (primitive 2n-th root of unity)
        let g = find_generator(q, DEGREE);

        // Compute twiddle factors
        let tf = compute_twiddle_factors::<DEGREE>(&class, g, false);
        let itf = compute_twiddle_factors::<DEGREE>(&class, g, true);

        // Precompute Shoup values for twiddle factors
        let mut tf_shoup = [0u64; DEGREE];
        for (i, &twiddle) in tf.iter().enumerate() {
            tf_shoup[i] = class.precompute_shoup(twiddle);
        }

        let mut itf_shoup = [0u64; DEGREE];
        for (i, &twiddle) in itf.iter().enumerate() {
            itf_shoup[i] = class.precompute_shoup(twiddle);
        }

        // Compute normalization factor (inverse of DEGREE)
        let inv_n = class.modinv(DEGREE as u64);
        let inv_n_shoup = class.precompute_shoup(inv_n);

        Arc::new(Self {
            class,
            inv_n,
            inv_n_shoup,
            tf,
            tf_shoup,
            itf,
            itf_shoup,
        })
    }

    /// Get the modulus for this context
    pub fn modulus(&self) -> u64 {
        self.class.q()
    }

    /// Get the polynomial degree for this context
    pub fn degree(&self) -> usize {
        DEGREE
    }

    /// Get the generator used for this context (for debugging/verification)
    pub fn generator(&self) -> u64 {
        // Reconstruct generator from first non-trivial twiddle factor
        // tf[1] = g^(bit_reverse(1)) = g^(2^(log2(DEGREE)-1)) for DEGREE > 1
        if DEGREE > 1 {
            // We'd need to compute discrete log to get g, but for debugging
            // we can verify the generator property instead
            find_generator(self.class.q(), DEGREE)
        } else {
            1 // Trivial case
        }
    }
}

/// Compute twiddle factors for NTT/INTT in bit-reversed order.
///
/// For forward NTT: twiddle factors are powers of g (primitive 2n-th root of unity)
/// For inverse NTT: twiddle factors are powers of g^(-1)
///
/// The factors are stored in bit-reversed order to match the NTT algorithm's
/// memory access pattern (Cooley-Tukey decimation-in-time).
///
/// # Arguments
/// * `class` - Modular arithmetic context
/// * `g` - Primitive 2n-th root of unity modulo q
/// * `is_inverse` - If true, compute factors for inverse NTT
fn compute_twiddle_factors<const DEGREE: usize>(
    class: &CongruenceClass,
    g: u64,
    is_inverse: bool,
) -> [u64; DEGREE] {
    let mut tf = [0u64; DEGREE];
    let mut tf_direct = [0u64; DEGREE];

    // Calculate bit-reversal length (log2 of DEGREE)
    let log_degree = (DEGREE.trailing_zeros()) as usize;

    // Use g or g^(-1) depending on direction
    let base = if is_inverse { class.modinv(g) } else { g };

    // Compute powers of base: base^0, base^1, base^2, ...
    tf_direct[0] = 1;
    for i in 1..DEGREE {
        tf_direct[i] = class.modmul(tf_direct[i - 1], base);
    }

    // Reorder in bit-reversed order for NTT algorithm
    for (i, tf_elem) in tf.iter_mut().enumerate() {
        let bit_rev_i = bit_reverse(i, log_degree);
        *tf_elem = tf_direct[bit_rev_i];
    }

    tf
}

/// Compute bit-reversal of a number within specified bit length.
///
/// Used to reorder twiddle factors for efficient NTT memory access.
///
/// # Arguments
/// * `number` - Input number to bit-reverse
/// * `bit_length` - Number of bits to consider (log2 of DEGREE)
///
/// # Examples
/// ```
/// # use rust_ntt::context::bit_reverse;
/// assert_eq!(bit_reverse(0b101, 3), 0b101); // 5 -> 5 (palindromic)
/// assert_eq!(bit_reverse(0b001, 3), 0b100); // 1 -> 4
/// assert_eq!(bit_reverse(0b010, 3), 0b010); // 2 -> 2 (palindromic)
/// ```
pub fn bit_reverse(number: usize, bit_length: usize) -> usize {
    let mut reversed = 0;
    for i in 0..bit_length {
        if (number >> i) & 1 != 0 {
            reversed |= 1 << (bit_length - 1 - i);
        }
    }
    reversed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{find_first_prime_down, find_first_prime_up};

    #[test]
    fn test_context_creation() {
        const N: usize = 4;
        let q = find_first_prime_up(5, N);
        let ctx = NttContext::<N>::new(q);

        assert_eq!(ctx.modulus(), q);
        assert_eq!(ctx.degree(), N);
        assert!(ctx.tf.len() == N);
        assert!(ctx.itf.len() == N);
    }

    #[test]
    fn test_context_larger_degree() {
        const N: usize = 1024;
        let q = find_first_prime_down(58, N);
        let ctx = NttContext::<N>::new(q);

        assert_eq!(ctx.modulus(), q);
        assert_eq!(ctx.degree(), N);
    }

    #[test]
    #[should_panic(expected = "DEGREE must be a power of 2")]
    fn test_invalid_degree() {
        let q = 17; // Just any prime
        let _ctx = NttContext::<6>::new(q); // 6 is not power of 2
    }

    #[test]
    #[should_panic(expected = "must satisfy q ≡ 1")]
    fn test_invalid_modulus() {
        let q = 19; // 19 ≡ 3 (mod 8), doesn't satisfy q ≡ 1 (mod 8) for DEGREE=4
        let _ctx = NttContext::<4>::new(q);
    }

    #[test]
    fn test_twiddle_factors_properties() {
        const N: usize = 8;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        // First twiddle factor should always be 1
        assert_eq!(ctx.tf[0], 1);
        assert_eq!(ctx.itf[0], 1);

        // Verify generator property: g^(2n) ≡ 1 (mod q)
        let g = ctx.generator();
        let g_to_2n = ctx.class.modexp(g, (2 * N) as u64);
        assert_eq!(g_to_2n, 1);

        // Verify g^n ≡ -1 (mod q) (primitive 2n-th root property)
        let g_to_n = ctx.class.modexp(g, N as u64);
        assert_eq!(g_to_n, q - 1); // -1 ≡ q-1 (mod q)
    }

    #[test]
    fn test_bit_reverse() {
        assert_eq!(bit_reverse(0, 3), 0);
        assert_eq!(bit_reverse(1, 3), 4); // 001 -> 100
        assert_eq!(bit_reverse(2, 3), 2); // 010 -> 010
        assert_eq!(bit_reverse(3, 3), 6); // 011 -> 110
        assert_eq!(bit_reverse(4, 3), 1); // 100 -> 001
        assert_eq!(bit_reverse(5, 3), 5); // 101 -> 101
        assert_eq!(bit_reverse(6, 3), 3); // 110 -> 011
        assert_eq!(bit_reverse(7, 3), 7); // 111 -> 111
    }

    #[test]
    fn test_inverse_relationship() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        // For each forward twiddle factor, verify inverse relationship
        for i in 0..N {
            let forward_power = ctx.tf[i];
            let inverse_power = ctx.itf[i];

            // forward * inverse should equal 1 (mod q)
            let product = ctx.class.modmul(forward_power, inverse_power);
            assert_eq!(product, 1, "tf[{}] * itf[{}] should equal 1", i, i);
        }
    }
}
