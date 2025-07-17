use rust_ntt::*;

fn main() {
    const N: usize = 4;
    let q = find_first_prime_up(5, N);

    let ring = PolyRing::<N>::new(q);

    let ax = ring.sample_random();
    let bx = ring.sample_random();

    let cx = ring.naive_negacyclic_convolution(&ax, &bx);

    dbg!(&cx);

    let cx_ntt = ring.ntt_negacyclic_convolution(&ax, &bx);

    dbg!(&cx_ntt);
}
