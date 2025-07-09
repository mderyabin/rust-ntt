use rand::{Rng, rng};
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

    for _ in 0..10 {
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
