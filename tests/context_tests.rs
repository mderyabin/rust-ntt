use proptest::prelude::*;
use rust_ntt::NttContext;
use rust_ntt::math::find_first_prime_up;

proptest! {
    #[test]
    fn context_creation_is_deterministic_degree_8(
        log_bits in 10u32..16u32,
    ) {
        const DEGREE: usize = 8;
        let q = find_first_prime_up(log_bits as usize, DEGREE);

        let ctx1 = NttContext::<DEGREE>::new(q);
        let ctx2 = NttContext::<DEGREE>::new(q);

        prop_assert_eq!(ctx1.modulus(), ctx2.modulus());
        prop_assert_eq!(ctx1.degree(), ctx2.degree());
        prop_assert_eq!(ctx1.tf(), ctx2.tf());
        prop_assert_eq!(ctx1.itf(), ctx2.itf());
    }
}

proptest! {
    #[test]
    fn generator_satisfies_primitive_root_property_degree_8(
        log_bits in 10u32..16u32,
    ) {
        const DEGREE: usize = 8;
        let q = find_first_prime_up(log_bits as usize, DEGREE);
        let ctx = NttContext::<DEGREE>::new(q);

        let g = ctx.generator();

        // g^(2n) ≡ 1 (mod q)
        let g_to_2n = ctx.class().modexp(g, (2 * DEGREE) as u64);
        prop_assert_eq!(g_to_2n, 1);

        // g^n ≡ -1 ≡ q-1 (mod q)
        let g_to_n = ctx.class().modexp(g, DEGREE as u64);
        prop_assert_eq!(g_to_n, q - 1);
    }
}

// Property: Inverse twiddle factors are multiplicative inverses
proptest! {
    #[test]
    fn inverse_twiddle_factors_are_correct_degree_8(
        log_bits in 10u32..15u32,
        index in 0usize..4usize
    ) {
        const DEGREE: usize = 8;
        let q = find_first_prime_up(log_bits as usize, DEGREE);
        let ctx = NttContext::<DEGREE>::new(q);

        // tf[i] * itf[i] ≡ 1 (mod q)
        let product = ctx.class().modmul(ctx.tf()[index], ctx.itf()[index]);
        prop_assert_eq!(product, 1);
    }
}
