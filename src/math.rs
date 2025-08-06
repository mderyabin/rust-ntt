// Number Theoretic functions necessary for implementation
// down to a scalar number operations

use crate::congruence::CongruenceClass;
use primal::{Sieve, is_prime};

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
    q.wrapping_sub(a) // q - a
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
#[inline]
pub fn barrett_precompute(q: u64) -> u128 {
    (1u128 << (2 * 63)) / (q as u128)
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

pub fn find_first_prime_up(logq: usize, n: usize) -> u64 {
    let mut q: u64 = (1u64 << logq) + 1;
    let m: u64 = (n as u64) << 1;

    while !is_prime(q) {
        q += m;
    }

    q
}

pub fn find_next_prime_up(prev_q: u64, n: usize) -> u64 {
    let m: u64 = (n as u64) << 1;
    let mut q = prev_q + m;

    while !is_prime(q) {
        q += m;
    }

    q
}

pub fn find_first_prime_down(logq: usize, n: usize) -> u64 {
    let m: u64 = (n as u64) << 1;
    let mut q: u64 = (1u64 << logq) + 1 - m;

    while !is_prime(q) {
        q -= m;
    }

    q
}

pub fn find_next_prime_down(prev_q: u64, n: usize) -> u64 {
    let m: u64 = (n as u64) << 1;
    let mut q = prev_q - m;

    while !is_prime(q) {
        q -= m;
    }

    q
}

pub fn find_primitive_root(q: u64) -> u64 {
    assert!(is_prime(q), "primitive root search: modulus must prime");

    let phi = q - 1;
    let logq = 64 - q.leading_zeros();

    let sieve = Sieve::new(1usize << (1 + logq / 2));
    let class = CongruenceClass::new(q);

    let phi_factorized = sieve.factor(phi as usize).unwrap();

    let mut gen_found = false;
    let mut r = 1;
    while !gen_found {
        r += 1;

        gen_found = phi_factorized
            .iter()
            .all(|(prime, _)| class.modexp(r, phi / (*prime as u64)) != 1);
    }

    r
}

pub fn find_generator(q: u64, n: usize) -> u64 {
    let class = CongruenceClass::new(q);

    let m = (n << 1) as u64;

    let g0 = find_primitive_root(q);
    let g = class.modexp(g0, (q - 1) / m);

    g
}
