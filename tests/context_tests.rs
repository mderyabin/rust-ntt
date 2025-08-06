/* use proptest::prelude::*;
use rust_ntt::NttContext;
use rust_ntt::math::find_first_prime_up;

// Property: Context creation should be deterministic
proptest! {
    #[test]
    fn context_creation_is_deterministic(
        log_bits in 10u32..20u32,
        degree_log in 2u32..12u32
    ) {
        let degree = 1usize << degree_log;
        let q = find_first_prime_up(log_bits as usize, degree);

        let ctx1 = NttContext::<degree>::new(q);
        let ctx2 = NttContext::new(q);

        prop_assert_eq!(ctx1.modulus(), ctx2.modulus());
        prop_assert_eq!(ctx1.degree(), ctx2.degree());

        // Twiddle factors should be identical
        prop_assert_eq!(ctx1.tf(), ctx2.tf());
        prop_assert_eq!(ctx1.itf(), ctx2.itf());
    }
}

// Property: Generator properties must hold
proptest! {
    #[test]
    fn generator_satisfies_primitive_root_property(
        log_bits in 10u32..16u32,
        degree_log in 2u32..8u32
    ) {
        let degree = 1usize << degree_log;
        let q = find_first_prime_up(log_bits as usize, degree);
        let ctx = NttContext::new(q);

        let g = ctx.generator();
        let class = &ctx.class;

        // g^(2n) ≡ 1 (mod q)
        let g_to_2n = class.modexp(g, (2 * degree) as u64);
        prop_assert_eq!(g_to_2n, 1);

        // g^n ≡ -1 ≡ q-1 (mod q)
        let g_to_n = class.modexp(g, degree as u64);
        prop_assert_eq!(g_to_n, q - 1);
    }
}

// Property: Inverse twiddle factors are multiplicative inverses
proptest! {
    #[test]
    fn inverse_twiddle_factors_are_correct(
        log_bits in 10u32..15u32,
        degree_log in 2u32..6u32,
        index in 0usize..64usize
    ) {
        let degree = 1usize << degree_log;
        let q = find_first_prime_up(log_bits as usize, degree);
        let ctx = NttContext::new(q);

        let idx = index % degree;

        // tf[i] * itf[i] ≡ 1 (mod q)
        let product = ctx.class.modmul(ctx.tf[idx], ctx.itf[idx]);
        prop_assert_eq!(product, 1);
    }
} */
