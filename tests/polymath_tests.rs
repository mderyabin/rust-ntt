use rust_ntt::*;

const Q: u64 = 741507920154517877;
// const Q : u64 = 1u64<<62-1;
const N: usize = 1usize << 2;

#[test]
fn test_add() {
    let ring = PolyRing::<N>::new(Q);

    let ax = ring.sample_random();
    let bx = ring.sample_random();

    let mut expected = [0u64; N];
    for i in 0..N {
        expected[i] = modadd_naive(ax[i], bx[i], Q)
    }

    assert_eq!(ring.add(&ax, &bx), expected);
}

#[test]
fn test_add_eq() {
    let ring = PolyRing::<N>::new(Q);

    let mut ax = ring.sample_random();
    let bx = ring.sample_random();

    let mut expected = [0u64; N];
    for i in 0..N {
        expected[i] = modadd_naive(ax[i], bx[i], Q)
    }

    ring.add_eq(&mut ax, &bx);

    assert_eq!(ax, expected);
}

#[test]
fn test_mul() {
    let ring = PolyRing::<N>::new(Q);

    let ax = ring.sample_random();
    let bx = ring.sample_random();

    let mut expected = [0u64; N];
    for i in 0..N {
        expected[i] = modmul_naive(ax[i], bx[i], Q)
    }

    assert_eq!(ring.mul(&ax, &bx), expected);
}

#[test]
fn test_mul_eq() {
    let ring = PolyRing::<N>::new(Q);

    let mut ax = ring.sample_random();
    let bx = ring.sample_random();

    // let ax_orig = ax.clone();

    let mut expected = [0u64; N];
    for i in 0..N {
        expected[i] = modmul_naive(ax[i], bx[i], Q)
    }

    ring.mul_eq(&mut ax, &bx);

    assert_eq!(ax, expected);
}

#[test]
fn test_ntt_intt() {
    let q = find_first_prime_down(58, N);
    let ring = PolyRing::<N>::new(q);

    for _ in 1..10 {
        let ax = ring.sample_random();

        let mut ax_ntt = ax.clone();

        ring.ntt_forward(&mut ax_ntt);
        ring.ntt_inverse(&mut ax_ntt);

        assert_eq!(ax, ax_ntt);
    }
}

#[test]
fn test_convolution() {
    let q = find_first_prime_down(58, N);
    let ring = PolyRing::<N>::new(q);

    for _ in 1..10 {
        let ax = ring.sample_random();
        let bx = ring.sample_random();

        let cx = ring.ntt_negacyclic_convolution(&ax, &bx);

        let expected = ring.naive_negacyclic_convolution(&ax, &bx);

        assert_eq!(cx, expected);
    }
}

#[test]
fn test_ntt_intt_shoup() {
    let q = find_first_prime_down(58, N);
    let ring = PolyRing::<N>::new(q);

    for _ in 1..10 {
        let ax = ring.sample_random();

        let mut ax_ntt = ax.clone();

        ring.ntt_forward_shoup(&mut ax_ntt);
        ring.ntt_inverse_shoup(&mut ax_ntt);

        assert_eq!(ax, ax_ntt);
    }
}


#[test]
fn test_convolution_shoup() {
    let q = find_first_prime_down(58, N);
    let ring = PolyRing::<N>::new(q);

    for _ in 1..10 {
        let ax = ring.sample_random();
        let bx = ring.sample_random();

        let cx = ring.ntt_negacyclic_convolution_shoup(&ax, &bx);

        let expected = ring.naive_negacyclic_convolution(&ax, &bx);

        assert_eq!(cx, expected);
    }
}