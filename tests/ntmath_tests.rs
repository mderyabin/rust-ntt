use rand::{Rng, rng};
use rust_ntt::math::{
    find_first_prime_up, find_generator, find_next_prime_up, modnegate,
};
use rust_ntt::*;

// const Q : u64 = 741507920154517877;
const Q: u64 = 1u64 << 62 - 1;

#[test]
fn test_modadd_simple() {
    let mut generator = rng();

    let a: u64 = generator.random_range(1..Q);
    let b: u64 = generator.random_range(1..Q);

    let expected = modadd_naive(a, b, Q);

    assert_eq!(modadd(a, b, Q), expected);
}

#[test]
fn test_modadd_zero() {
    let mut generator = rng();
    let a: u64 = generator.random_range(1..Q);

    assert_eq!(modadd(a, 0, Q), a);
    assert_eq!(modadd(0, a, Q), a);
}

#[test]
fn test_modadd_one() {
    let mut generator = rng();
    let a: u64 = generator.random_range(1..Q);

    assert_eq!(modadd(a, Q - 1, Q), a - 1);
    assert_eq!(modadd(Q - 1, a, Q), a - 1);
}

#[test]
fn test_modsub_simple() {
    let mut generator = rng();
    let a: u64 = generator.random_range(1..Q);
    let b: u64 = generator.random_range(1..Q);

    let expected = (Q + a - b) % Q;

    assert_eq!(modsub(a, b, Q), expected);
}

#[test]
fn test_modsub_zero() {
    let mut generator = rng();
    let a: u64 = generator.random_range(1..Q);

    assert_eq!(modsub(a, 0, Q), a);
    assert_eq!(modsub(0, a, Q), Q - a);
}

#[test]
fn test_modsub_one() {
    let mut generator = rng();
    let a: u64 = generator.random_range(1..Q);

    assert_eq!(modsub(a, Q - 1, Q), (a + 1) % Q);
    assert_eq!(modsub(Q - 1, a, Q), Q - (a + 1));
}

#[test]
fn test_modmul_naive() {
    let mut generator = rng();

    let a: u64 = generator.random_range(1..Q);
    // let b : u64 = generator.random_range(1..Q);

    assert_eq!(modmul_naive(a, 2, Q), (a * 2) % Q);
    assert_eq!(modmul_naive(a, 1, Q), a);
    assert_eq!(modmul_naive(2, a, Q), (a * 2) % Q);
    assert_eq!(modmul_naive(1, a, Q), a);
    assert_eq!(modmul_naive(a, 0, Q), 0);
    assert_eq!(modmul_naive(0, a, Q), 0);
    assert_eq!(modmul_naive(a, Q - 1, Q), modnegate(a, Q));
    assert_eq!(modmul_naive(Q - 1, a, Q), modnegate(a, Q));
}

#[test]
fn test_modmul_barrett_old() {
    let mut generator = rng();

    let a: u64 = generator.random_range(1..Q);
    let b: u64 = generator.random_range(1..Q);

    let (mu, logq) = barrett_precompute_old(Q);

    let expected = modmul_naive(a, b, Q);

    assert_eq!(modmul_barrett_old(a, b, Q, mu, logq), expected);
    assert_eq!(modmul_barrett_old(b, a, Q, mu, logq), expected);
    assert_eq!(modmul_barrett_old(a, 2, Q, mu, logq), (a * 2) % Q);
    assert_eq!(modmul_barrett_old(a, 1, Q, mu, logq), a);
    assert_eq!(modmul_barrett_old(2, a, Q, mu, logq), (a * 2) % Q);
    assert_eq!(modmul_barrett_old(1, a, Q, mu, logq), a);
    assert_eq!(modmul_barrett_old(a, 0, Q, mu, logq), 0);
    assert_eq!(modmul_barrett_old(0, a, Q, mu, logq), 0);
    assert_eq!(modmul_barrett_old(a, Q - 1, Q, mu, logq), modnegate(a, Q));
    assert_eq!(modmul_barrett_old(Q - 1, a, Q, mu, logq), modnegate(a, Q));
}

#[test]
fn test_modmul_barrett_old_eq() {
    let mut generator = rng();

    let a_in: u64 = generator.random_range(1..Q);
    let b: u64 = generator.random_range(1..Q);
    let mut a;

    let (mu, logq) = barrett_precompute_old(Q);

    let expected = modmul_naive(a_in, b, Q);

    a = a_in;
    modmul_barrett_old_eq(&mut a, b, Q, mu, logq);
    assert_eq!(a, expected);

    a = b;
    modmul_barrett_old_eq(&mut a, a_in, Q, mu, logq);
    assert_eq!(a, expected);

    a = a_in;
    modmul_barrett_old_eq(&mut a, 2, Q, mu, logq);
    assert_eq!(a, (a_in * 2) % Q);

    a = a_in;
    modmul_barrett_old_eq(&mut a, 1, Q, mu, logq);
    assert_eq!(a, a_in);

    a = 2;
    modmul_barrett_old_eq(&mut a, a_in, Q, mu, logq);
    assert_eq!(a, (a_in * 2) % Q);

    a = 1;
    modmul_barrett_old_eq(&mut a, a_in, Q, mu, logq);
    assert_eq!(a, a_in);

    a = a_in;
    modmul_barrett_old_eq(&mut a, 0, Q, mu, logq);
    assert_eq!(a, 0);

    a = 0;
    modmul_barrett_old_eq(&mut a, a_in, Q, mu, logq);
    assert_eq!(a, 0);

    a = a_in;
    modmul_barrett_old_eq(&mut a, Q - 1, Q, mu, logq);
    assert_eq!(a, modnegate(a_in, Q));

    a = Q - 1;
    modmul_barrett_old_eq(&mut a, a_in, Q, mu, logq);
    assert_eq!(a, modnegate(a_in, Q));
}

#[test]
fn test_modmul_barrett() {
    let mut generator = rng();

    let a: u64 = generator.random_range(1..Q);
    let b: u64 = generator.random_range(1..Q);

    let mu = barrett_precompute(Q);

    let expected = modmul_naive(a, b, Q);

    assert_eq!(modmul_barrett(a, b, Q, mu), expected);
    assert_eq!(modmul_barrett(b, a, Q, mu), expected);
    assert_eq!(modmul_barrett(a, 2, Q, mu), (a * 2) % Q);
    assert_eq!(modmul_barrett(a, 1, Q, mu), a);
    assert_eq!(modmul_barrett(2, a, Q, mu), (a * 2) % Q);
    assert_eq!(modmul_barrett(1, a, Q, mu), a);
    assert_eq!(modmul_barrett(a, 0, Q, mu), 0);
    assert_eq!(modmul_barrett(0, a, Q, mu), 0);
    assert_eq!(modmul_barrett(a, Q - 1, Q, mu), modnegate(a, Q));
    assert_eq!(modmul_barrett(Q - 1, a, Q, mu), modnegate(a, Q));
}

#[test]
fn test_modmul_barrett_eq() {
    let mut generator = rng();

    let a_in: u64 = generator.random_range(1..Q);
    let b: u64 = generator.random_range(1..Q);
    let mut a;

    let mu = barrett_precompute(Q);

    let expected = modmul_naive(a_in, b, Q);

    a = a_in;
    modmul_barrett_eq(&mut a, b, Q, mu);
    assert_eq!(a, expected);

    a = b;
    modmul_barrett_eq(&mut a, a_in, Q, mu);
    assert_eq!(a, expected);

    a = a_in;
    modmul_barrett_eq(&mut a, 2, Q, mu);
    assert_eq!(a, (a_in * 2) % Q);

    a = a_in;
    modmul_barrett_eq(&mut a, 1, Q, mu);
    assert_eq!(a, a_in);

    a = 2;
    modmul_barrett_eq(&mut a, a_in, Q, mu);
    assert_eq!(a, (a_in * 2) % Q);

    a = 1;
    modmul_barrett_eq(&mut a, a_in, Q, mu);
    assert_eq!(a, a_in);

    a = a_in;
    modmul_barrett_eq(&mut a, 0, Q, mu);
    assert_eq!(a, 0);

    a = 0;
    modmul_barrett_eq(&mut a, a_in, Q, mu);
    assert_eq!(a, 0);

    a = a_in;
    modmul_barrett_eq(&mut a, Q - 1, Q, mu);
    assert_eq!(a, modnegate(a_in, Q));

    a = Q - 1;
    modmul_barrett_eq(&mut a, a_in, Q, mu);
    assert_eq!(a, modnegate(a_in, Q));
}

#[test]
fn test_modmul_barrett_struct() {
    let mut generator = rng();

    let a: u64 = generator.random_range(1..Q);
    let b: u64 = generator.random_range(1..Q);

    let barrett64 = CongruenceClass::new(Q);

    let expected = modmul_naive(a, b, Q);

    assert_eq!(barrett64.modmul(a, b), expected);
    assert_eq!(barrett64.modmul(b, a), expected);
    assert_eq!(barrett64.modmul(a, 2), (a * 2) % Q);
    assert_eq!(barrett64.modmul(a, 1), a);
    assert_eq!(barrett64.modmul(2, a), (a * 2) % Q);
    assert_eq!(barrett64.modmul(1, a), a);
    assert_eq!(barrett64.modmul(a, 0), 0);
    assert_eq!(barrett64.modmul(0, a), 0);
    assert_eq!(barrett64.modmul(a, Q - 1), modnegate(a, Q));
    assert_eq!(barrett64.modmul(Q - 1, a), modnegate(a, Q));
}

#[test]
fn test_modmul_struct_vs_naive() {
    let mut rng = rng();
    let q = 741507920154517877;
    let class = CongruenceClass::new(q);

    for _ in 0..1000 {
        let a: u64 = rng.random_range(1..q);
        let b: u64 = rng.random_range(1..q);

        let fast = class.modmul(a, b); // Barrett
        let slow = modmul_naive(a, b, q); // Reference

        assert_eq!(
            fast, slow,
            "a = {}, b = {}, got {}, expected {}",
            a, b, fast, slow
        );
    }
}

#[test]
fn test_modmul_eq_struct_vs_naive() {
    let mut rng = rng();
    let q = 741507920154517877;
    let class = CongruenceClass::new(q);

    for _ in 0..1000 {
        let mut a: u64 = rng.random_range(1..q);
        let b: u64 = rng.random_range(1..q);

        let a_copy = a;

        let slow = modmul_naive(a, b, q); // Reference
        class.modmul_eq(&mut a, b); // Barrett

        assert_eq!(
            a, slow,
            "a = {}, b = {}, got {}, expected {}",
            a_copy, b, a, slow
        );
    }
}

#[test]
fn test_generator() {
    // warning: slow!
    let n = 1usize << 10;
    let mut q = find_first_prime_up(52, n);

    for _ in 1..10 {
        q = find_next_prime_up(q, n);

        let class = CongruenceClass::new(q);

        let g = find_generator(q, n);

        let check1 = class.modexp(g, n as u64);
        let check2 = class.modexp(g, (n << 1) as u64);

        assert_eq!(
            check1,
            q - 1,
            "power n g = {}, got {}, expected {}, q = {}",
            g,
            check1,
            q - 1,
            q
        );
        assert_eq!(
            check2, 1,
            "power 2n g = {}, got {}, expected {}, q = {}",
            g, check2, 1, q
        );
    }
}

#[test]
fn test_modmul_shoup_struct() {
    let mut generator = rng();

    let class = CongruenceClass::new(Q);

    for _ in 0..100 {
        let a: u64 = generator.random_range(1..Q);
        let b: u64 = generator.random_range(1..Q);

        let expected = modmul_naive(a, b, Q);

        let prec = class.precompute_shoup(b);

        assert_eq!(class.modmul_shoup(a, b, prec), expected);
    }
}

#[test]
fn test_modmul_shoup_eq_struct() {
    let mut generator = rng();

    let class = CongruenceClass::new(Q);

    for _ in 0..100 {
        let mut a: u64 = generator.random_range(1..Q);
        let b: u64 = generator.random_range(1..Q);

        let expected = modmul_naive(a, b, Q);

        let prec = class.precompute_shoup(b);

        class.modmul_shoup_eq(&mut a, b, prec);

        assert_eq!(a, expected);
    }
}

#[test]
fn test_modmul_shoup_as64_struct() {
    let mut generator = rng();

    let class = CongruenceClass::new(Q);

    for _ in 0..100 {
        let a: u64 = generator.random_range(1..Q);
        let b: u64 = generator.random_range(1..Q);

        let expected = modmul_naive(a, b, Q);

        let prec = class.precompute_shoup(b);

        assert_eq!(class.modmul_shoup_as64(a, b, prec), expected);
    }
}
