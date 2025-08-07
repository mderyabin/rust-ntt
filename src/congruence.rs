/// Precomputed modular arithmetic context using Barrett reduction.
///
/// This struct encapsulates all precomputed values needed for efficient modular
/// arithmetic operations modulo a fixed prime `q`. Barrett reduction avoids
/// expensive division operations by precomputing reduction parameters.
///
/// # Mathematical Background
///
/// Barrett reduction replaces `a * b mod q` with:
/// 1. Compute `μ = ⌊2^(2*logq) / q⌋` (precomputed)
/// 2. Estimate quotient: `⌊(a*b) >> (logq-2) * μ >> (logq+2)⌋`
/// 3. Compute remainder: `a*b - estimated_quotient * q`
/// 4. Apply final correction if needed
///
/// # Examples
///
/// ```rust
/// use rust_ntt::CongruenceClass;
///
/// let q = 97; // Prime modulus
/// let class = CongruenceClass::new(q);
///
/// let result = class.modmul(15, 23); // 15 * 23 mod 97
/// assert_eq!(result, (15 * 23) % 97);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct CongruenceClass {
    mu: u64,   // Barrett parameter μ = ⌊2^(2*logq) / q⌋
    q: u64,    // Prime modulus
    logq: u64, // Bit length of q (⌈log₂(q)⌉)
}

// Here are getters
impl CongruenceClass {
    #[inline]
    pub fn q(&self) -> u64 {
        self.q
    }
}

impl CongruenceClass {
    /// Creates a new modular arithmetic context for the given modulus.
    ///
    /// Precomputes Barrett reduction parameters for efficient modular operations.
    ///
    /// # Arguments
    /// * `q` - Prime modulus, must satisfy 2 ≤ q < 2^63
    ///
    /// # Panics
    /// * If q < 2 or q ≥ 2^63
    pub fn new(q: u64) -> Self {
        assert!(q >= 2, "modulus must be ≥ 2");
        assert!(q < (1u64 << 63), "modulus must be < 2^63");

        // let mu: u128 = (1u128 << (2 * 63)) / (q as u128);

        let logq: u64 = 64 - (q.leading_zeros() as u64);
        let mu: u64 = ((1u128 << (2 * logq)) / (q as u128)) as u64;

        Self { q, mu, logq }
    }
    // mu = (2^126 / q)

    /// Precomputes Shoup parameter for repeated multiplications with `b`.
    ///
    /// When multiplying many values by the same `b`, precomputing this parameter
    /// and using `modmul_shoup()` is faster than repeated `modmul()` calls.
    ///
    /// # Returns
    /// Shoup parameter `⌊b * 2^64 / q⌋` for use with `modmul_shoup()`
    ///
    /// # Example
    /// ```rust
    /// let prec = class.precompute_shoup(42);
    /// let result = class.modmul_shoup(15, 42, prec); // Faster than modmul(15, 42)
    /// ```
    pub fn precompute_shoup(&self, b: u64) -> u64 {
        (((b as u128) << 64) / (self.q as u128)) as u64
    }

    /// Alternative Shoup multiplication using 64-bit arithmetic.
    ///
    /// **Experimental**: This variant attempts to use only 64-bit operations
    /// instead of 128-bit arithmetic, but currently performs worse in practice
    /// due to additional overflow checks. Kept for research purposes.
    ///
    /// # Note
    /// Despite being faster in microbenchmarks, this performs poorly in actual
    /// NTT operations due to wrapping_mul overhead and correctness concerns.
    #[inline]
    pub fn modmul_shoup_as64(&self, a: u64, b: u64, b_prec: u64) -> u64 {
        // эта функция как альтернатива modmul_shoup
        // по идее mul и tmp (умножение на q) можно делать в u64 так как верхняя часть обрезается.
        // это должно быть намного быстрее.
        // бенчмарк этой функции действительно быстрее,
        // однако когда делаю бенчмарк ntt на ее основе, он значительно медленнее.
        // скорее всего из-за проверок и wrapping_mul
        // однако если делать без wrapping_mul то тесты падают и функция возвращает неверное значение
        let mul = a.wrapping_mul(b);
        let tmp =
            ((((a as u128) * (b_prec as u128)) >> 64) as u64).wrapping_mul(self.q);

        let r = mul.wrapping_sub(tmp);

        if r < self.q {
            r
        } else {
            r.wrapping_sub(self.q)
        }
    }

    /// Optimized modular multiplication using precomputed Shoup parameter.
    ///
    /// Faster than `modmul()` when multiplying many values by the same `b`.
    /// Must use `b_prec` from `precompute_shoup(b)`.
    ///
    /// # Arguments
    /// * `a` - First operand
    /// * `b` - Second operand (same as used in `precompute_shoup()`)
    /// * `b_prec` - Precomputed Shoup parameter from `precompute_shoup(b)`
    #[inline]
    pub fn modmul_shoup(&self, a: u64, b: u64, b_prec: u64) -> u64 {
        // let mul = a * b;
        // let tmp = ((((a as u128) * (b_prec as u128)) >> 64) as u64) * (self.q as u64);

        let mul = (a as u128) * (b as u128);
        let tmp = (((a as u128) * (b_prec as u128)) >> 64) * (self.q as u128);

        // let r = mul.wrapping_sub(tmp);
        let r = (mul - tmp) as u64;

        if r < self.q {
            r
        } else {
            r.wrapping_sub(self.q)
        }
    }

    #[inline]
    pub fn modmul_shoup_eq(&self, a: &mut u64, b: u64, b_prec: u64) {
        // let mul = (*a as u64).wrapping_mul(b as u64);
        // let tmp = ((((*a as u128) * (b_prec as u128)) >> 64) as u64).wrapping_mul(self.q);
        let mul = (*a as u128) * (b as u128);
        let tmp = (((*a as u128) * (b_prec as u128)) >> 64) * (self.q as u128);

        // let r = mul.wrapping_sub(tmp);
        let r = (mul - tmp) as u64;

        *a = if r < self.q {
            r
        } else {
            r.wrapping_sub(self.q)
        };
    }

    /// Fast modular multiplication: (a * b) mod q.
    ///
    /// Uses Barrett reduction to avoid expensive division operations.
    /// This is the primary method for modular multiplication.
    ///
    /// # Arguments
    /// * `a`, `b` - Operands, should be in range [0, q)
    ///
    /// # Returns
    /// `(a * b) mod q`
    #[inline]
    pub fn modmul(&self, a: u64, b: u64) -> u64 {
        let mul = (a as u128) * (b as u128);

        let tmp1 = mul >> (self.logq - 2); // (ab / 2^62)
        let tmp2 = (tmp1 * (self.mu as u128)) >> (self.logq + 2);
        // (ab / 2^62) * (2^126 / q) / 2^64 = (ab 2^64 / q) / 2^64 = floor(ab/q)

        let r = (mul.wrapping_sub(tmp2 * (self.q as u128))) as u64;
        // ab - floor(ab/q) * q = ab mod q

        // return if r < self.q { r } else { r.wrapping_sub(self.q) };
        if r < self.q {
            r
        } else {
            r.wrapping_sub(self.q)
        }
    }

    #[inline]
    pub fn modsquare(&self, a: u64) -> u64 {
        let mul = (a as u128) * (a as u128);

        let tmp1 = mul >> (self.logq - 2); // (ab / 2^62)
        let tmp2 = (tmp1 * (self.mu as u128)) >> (self.logq + 2);
        // (ab / 2^62) * (2^126 / q) / 2^64 = (ab 2^64 / q) / 2^64 = floor(ab/q)

        let r = (mul.wrapping_sub(tmp2 * (self.q as u128))) as u64;
        // ab - floor(ab/q) * q = ab mod q

        // return if r < self.q { r } else { r.wrapping_sub(self.q) };
        if r < self.q {
            r
        } else {
            r.wrapping_sub(self.q)
        }
    }

    #[inline]
    pub fn modmul_eq(&self, a: &mut u64, b: u64) {
        let mul = (*a as u128) * (b as u128);

        let tmp1 = mul >> (self.logq - 2); // (ab / 2^62)
        let tmp2 = (tmp1 * (self.mu as u128)) >> (self.logq + 2);

        let r = (mul.wrapping_sub(tmp2 * (self.q as u128))) as u64;

        *a = if r < self.q {
            r
        } else {
            r.wrapping_sub(self.q)
        };
    }

    #[inline]
    pub fn modsquare_eq(&self, a: &mut u64) {
        let mul = (*a as u128) * (*a as u128);

        let tmp1 = mul >> (self.logq - 2); // (ab / 2^62)
        let tmp2 = (tmp1 * (self.mu as u128)) >> (self.logq + 2);

        let r = (mul.wrapping_sub(tmp2 * (self.q as u128))) as u64;

        *a = if r < self.q {
            r
        } else {
            r.wrapping_sub(self.q)
        };
    }

    #[inline]
    pub fn modadd(&self, a: u64, b: u64) -> u64 {
        let t = a + b;
        if t <= self.q {
            t
        } else {
            t.wrapping_sub(self.q)
        }
    }

    #[inline]
    pub fn modadd_eq(&self, a: &mut u64, b: u64) {
        let t = *a + b;
        *a = if t <= self.q {
            t
        } else {
            t.wrapping_sub(self.q)
        };
    }

    #[inline]
    pub fn modsub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a.wrapping_sub(b)
        } else {
            (self.q + a).wrapping_sub(b)
        }
    }

    #[inline]
    pub fn modsub_eq(&self, a: &mut u64, b: u64) {
        *a = if *a >= b {
            (*a).wrapping_sub(b)
        } else {
            (self.q + *a).wrapping_sub(b)
        };
    }

    #[inline]
    pub fn modneg(&self, a: u64) -> u64 {
        // process edge case when a = 0
        if a == 0 { 0 } else { self.q.wrapping_sub(a) }
    }

    #[inline]
    pub fn modneg_eq(&self, a: &mut u64) {
        (*a) = self.q.wrapping_sub(*a);
    }

    #[inline]
    pub fn modexp(&self, a: u64, e: u64) -> u64 {
        let mut base = a;
        let mut exp = e;

        let mut result = 1u64;

        while exp > 0 {
            if exp % 2 == 1 {
                result = self.modmul(result, base);
            }
            base = self.modsquare(base);
            exp >>= 1;
        }

        result
    }

    pub fn modexp_eq(&self, a: &mut u64, e: u64) {
        let mut base = *a;
        let mut exp = e;

        *a = 1u64;

        while exp > 0 {
            if exp % 2 == 1 {
                self.modmul_eq(&mut *a, base);
            }
            self.modsquare_eq(&mut base);
            exp >>= 1;
        }
    }

    pub fn modinv(&self, a: u64) -> u64 {
        self.modexp(a, self.q - 2)
    }

    pub fn modinv_eq(&self, a: &mut u64) {
        self.modexp_eq(&mut *a, self.q - 2);
    }
}
