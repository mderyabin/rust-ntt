// Number Theoretic functions necessary for implementation
// down to a scalar number operations

/** naive version to test performance **/
#[inline]
pub fn modadd_naive(a: u64, b: u64, q: u64) -> u64 {
    (a + b) % q
}

#[inline]
pub fn modmul_naive(a: u64, b: u64, q: u64) -> u64 {
    let prod: u128 = ((a as u128) * (b as u128)) % (q as u128);
    prod as u64
}

/** efficient section **/
#[inline]
pub fn modnegate(a: u64, q: u64) -> u64 {
    q.wrapping_sub(a)// q - a
}

#[inline]
pub fn modadd(a: u64, b: u64, q: u64) -> u64 {
    let t = a + b;
    if t < q { t } else { t.wrapping_sub(q) }
}

#[inline]
pub fn modsub(a: u64, b: u64, q: u64) -> u64 {
    if a >= b {
        a.wrapping_sub(b)
    } else {
        (q + a).wrapping_sub(b)
    }
}

// возможна 2ая версия когда mu принимает значение в 128 бит диапазоне, а logq фиксируется в 63
pub fn barrett_precompute_old(q: u64) -> (u64, u64) {
    let logq: u64 = if q == 0 {
        0
    } else {
        64 - (q.leading_zeros() as u64)
    };
    let mu: u64 = ((1u128 << (2 * logq)) / (q as u128)) as u64;
    (mu, logq)
}

// возможна 2ая версия когда mu принимает значение в 128 бит диапазоне, а logq фиксируется в 63
#[inline]
pub fn modmul_barrett_old(a: u64, b: u64, q: u64, mu: u64, logq: u64) -> u64 {
    let mul = (a as u128) * (b as u128);

    let tmp1 = mul >> (logq - 1);
    let tmp2 = (tmp1 * (mu as u128)) >> (logq + 1);

    let r = (mul.wrapping_sub(tmp2 * (q as u128))) as u64;

    if r < q { r } else { r.wrapping_sub(q) }
}

// возможна 2ая версия когда mu принимает значение в 128 бит диапазоне, а logq фиксируется в 63
#[inline]
pub fn modmul_barrett_old_eq(a: &mut u64, b: u64, q: u64, mu: u64, logq: u64) {
    let aa = *a as u128;
    let mul = aa * (b as u128);

    let tmp1 = mul >> (logq - 1);
    let tmp2 = (tmp1 * (mu as u128)) >> (logq + 1);

    let r = (mul - tmp2 * (q as u128)) as u64;

    *a = if r < q { r } else { r.wrapping_sub(q) };
}

// возможна 2ая версия когда mu принимает значение в 128 бит диапазоне, а logq фиксируется в 63
pub fn barrett_precompute(q: u64) -> u128 {
    let mu: u128 = (1u128 << (2 * 63)) / (q as u128);
    mu
}

// возможна 2ая версия когда mu принимает значение в 128 бит диапазоне, а logq фиксируется в 63
#[inline]
pub fn modmul_barrett(a: u64, b: u64, q: u64, mu: u128) -> u64 {
    let mul = (a as u128) * (b as u128);

    let tmp1 = mul >> 62;
    let tmp2 = (tmp1 * mu) >> 64;

    let r = (mul.wrapping_sub(tmp2 * (q as u128))) as u64;

    if r < q { r } else { r.wrapping_sub(q) }
}

// возможна 2ая версия когда mu принимает значение в 128 бит диапазоне, а logq фиксируется в 63
#[inline]
pub fn modmul_barrett_eq(a: &mut u64, b: u64, q: u64, mu: u128) {
    let aa = *a as u128;
    let mul = aa * (b as u128);

    let tmp1 = mul >> 62;
    let tmp2 = (tmp1 * mu) >> 64;

    let r = (mul - tmp2 * (q as u128)) as u64;

    *a = if r < q { r } else { r.wrapping_sub(q) };
}

#[derive(Debug, Clone, Copy)]
#[repr(C, align(128))]
pub struct CongruenceClass {
    // mu: u128,
    mu: u64,
    pub q: u64,
    logq: u64, // можно меньше, но с выравниванием может быть лучше 64 бита
}

impl CongruenceClass {
    pub fn new(q: u64) -> Self {
        assert!(q >= 2, "modulus must be ≥ 2");
        assert!(q < (1u64 << 63), "modulus must be < 2^63");

        // let mu: u128 = (1u128 << (2 * 63)) / (q as u128);

        let logq: u64 =  64 - (q.leading_zeros() as u64);
        let mu: u64 = ((1u128 << (2 * logq)) / (q as u128)) as u64;


        Self { q, mu, logq }
    }
    // mu = (2^126 / q)

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

    pub fn modadd(&self, a: u64, b: u64) -> u64 {
        let t = a + b;
        if t <= self.q {
            t
        } else {
            t.wrapping_sub(self.q)
        }
    }

    pub fn modadd_eq(&self, a: &mut u64, b: u64) {
        let t = *a + b;
        *a = if t <= self.q {
            t
        } else {
            t.wrapping_sub(self.q)
        };
    }

    pub fn modsub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a.wrapping_sub(b)
        } else {
            (self.q + a).wrapping_sub(b)
        }
    }

    pub fn modsub_eq(&self, a: &mut u64, b: u64) {
        *a = if *a >= b {
            (*a).wrapping_sub(b)
        } else {
            (self.q + *a).wrapping_sub(b)
        };
    }

    pub fn modneg(&self, a: u64) -> u64 {
        self.q.wrapping_sub(a)
    }

    pub fn modneg_eq(&self, a: &mut u64) {
        (*a) = self.q.wrapping_sub(*a);
    }
}
