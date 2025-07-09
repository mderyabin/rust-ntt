use rust_ntt::*;

const Q: u64 = 741507920154517877;
// const Q : u64 = 1u64<<62-1;
const N: usize = 1usize << 2;

#[test]
fn test_add() {
    let ring = PolyRing::new(Q, N);

    let ax: Vec<u64> = ring.sample_random();
    let bx: Vec<u64> = ring.sample_random();

    let expected: Vec<u64> = ax
        .iter()
        .zip(bx.clone())
        .map(|(&a, b)| modadd_naive(a, b, Q))
        .collect();

    assert_eq!(ring.add(&ax, &bx), expected);
}

#[test]
fn test_add_eq() {
    let ring = PolyRing::new(Q, N);

    let mut ax: Vec<u64> = ring.sample_random();
    let bx: Vec<u64> = ring.sample_random();

    let expected: Vec<u64> = ax
        .iter()
        .zip(bx.clone())
        .map(|(&a, b)| modadd_naive(a, b, Q))
        .collect();

    ring.add_eq(&mut ax, &bx);

    assert_eq!(ax, expected);
}

#[test]
fn test_mul() {
    let ring = PolyRing::new(Q, N);

    let ax: Vec<u64> = ring.sample_random();
    let bx: Vec<u64> = ring.sample_random();

    let expected: Vec<u64> = ax
        .iter()
        .zip(bx.clone())
        .map(|(&a, b)| modmul_naive(a, b, Q))
        .collect();

    assert_eq!(ring.mul(&ax, &bx), expected);
}

#[test]
fn test_mul_eq() {
    let ring = PolyRing::new(Q, N);

    let mut ax: Vec<u64> = ring.sample_random();
    let bx: Vec<u64> = ring.sample_random();

    let ax_orig = ax.clone();

    let expected: Vec<u64> = ax
        .iter()
        .zip(&bx)
        .map(|(&a, &b)| modmul_naive(a, b, Q))
        .collect();

    ring.mul_eq(&mut ax, &bx);

    assert_eq!(ax_orig, expected);
}
