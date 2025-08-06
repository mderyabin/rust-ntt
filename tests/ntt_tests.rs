use proptest::prelude::*;
use rust_ntt::math::find_first_prime_up;
use rust_ntt::{NttContext, NttPolynomial};
use std::sync::Arc;

// Custom strategy for generating valid contexts
fn valid_context_strategy() -> impl Strategy<Value = Arc<NttContext<4>>> {
    (10u32..16u32).prop_map(|log_bits| {
        let q = find_first_prime_up(log_bits as usize, 4);
        NttContext::<4>::new(q)
    })
}

// Property: NTT followed by INTT is identity
proptest! {
    #[test]
    fn ntt_intt_is_identity(
        coeffs in prop::array::uniform4(1u64..1000u64),
        ctx in valid_context_strategy()
    ) {
        // Ensure coefficients are in valid range
        let valid_coeffs = coeffs.map(|c| c % ctx.modulus());

        let original = NttPolynomial::from_coeffs(valid_coeffs, Arc::clone(&ctx));
        let mut test_poly = original.clone();

        test_poly.ntt_forward();
        test_poly.ntt_inverse();

        prop_assert_eq!(test_poly.coeffs(), original.coeffs());
    }
}

// Property: Shoup NTT/INTT is equivalent to regular version
proptest! {
    #[test]
    fn shoup_ntt_equivalent_to_regular(
        coeffs in prop::array::uniform4(1u64..1000u64),
        ctx in valid_context_strategy()
    ) {
        let valid_coeffs = coeffs.map(|c| c % ctx.modulus());

        let mut regular = NttPolynomial::from_coeffs(valid_coeffs, Arc::clone(&ctx));
        let mut shoup = NttPolynomial::from_coeffs(valid_coeffs, Arc::clone(&ctx));

        // Forward transforms should be equivalent
        regular.ntt_forward();
        shoup.ntt_forward_shoup();
        prop_assert_eq!(regular.coeffs(), shoup.coeffs());

        // Inverse transforms should be equivalent
        regular.ntt_inverse();
        shoup.ntt_inverse_shoup();
        prop_assert_eq!(regular.coeffs(), shoup.coeffs());
    }
}

// Property: Addition is commutative and associative
proptest! {
    #[test]
    fn addition_properties(
        coeffs_a in prop::array::uniform4(0u64..1000u64),
        coeffs_b in prop::array::uniform4(0u64..1000u64),
        coeffs_c in prop::array::uniform4(0u64..1000u64),
        ctx in valid_context_strategy()
    ) {
        let a = NttPolynomial::from_coeffs(
            coeffs_a.map(|c| c % ctx.modulus()),
            Arc::clone(&ctx)
        );
        let b = NttPolynomial::from_coeffs(
            coeffs_b.map(|c| c % ctx.modulus()),
            Arc::clone(&ctx)
        );
        let c = NttPolynomial::from_coeffs(
            coeffs_c.map(|c| c % ctx.modulus()),
            Arc::clone(&ctx)
        );

        // Commutativity: a + b = b + a
        let ab = &a + &b;
        let ba = &b + &a;
        prop_assert_eq!(ab.coeffs(), ba.coeffs());

        // Associativity: (a + b) + c = a + (b + c)
        let abc_left = &(&a + &b) + &c;
        let abc_right = &a + &(&b + &c);
        prop_assert_eq!(abc_left.coeffs(), abc_right.coeffs());
    }
}

// Property: Zero is additive identity
proptest! {
    #[test]
    fn zero_is_additive_identity(
        coeffs in prop::array::uniform4(0u64..1000u64),
        ctx in valid_context_strategy()
    ) {
        let a = NttPolynomial::from_coeffs(
            coeffs.map(|c| c % ctx.modulus()),
            Arc::clone(&ctx)
        );
        let zero = NttPolynomial::zero(Arc::clone(&ctx));

        // a + 0 = a
        let a_plus_zero = &a + &zero;
        prop_assert_eq!(a_plus_zero.coeffs(), a.coeffs());

        // 0 + a = a
        let zero_plus_a = &zero + &a;
        prop_assert_eq!(zero_plus_a.coeffs(), a.coeffs());
    }
}

// Property: Negacyclic convolution properties
proptest! {
    #[test]
    fn convolution_properties(
        coeffs_a in prop::array::uniform4(1u64..100u64),
        coeffs_b in prop::array::uniform4(1u64..100u64),
        ctx in valid_context_strategy()
    ) {
        let a = NttPolynomial::from_coeffs(
            coeffs_a.map(|c| c % ctx.modulus()),
            Arc::clone(&ctx)
        );
        let b = NttPolynomial::from_coeffs(
            coeffs_b.map(|c| c % ctx.modulus()),
            Arc::clone(&ctx)
        );

        // Both convolution methods should give same result
        let result_regular = a.negacyclic_convolution(&b);
        let result_shoup = a.negacyclic_convolution_shoup(&b);

        prop_assert_eq!(result_regular.coeffs(), result_shoup.coeffs());

        // Convolution via multiplication operator should match
        let result_mul = &a * &b;
        prop_assert_eq!(result_regular.coeffs(), result_mul.coeffs());
    }
}

// Property: Negacyclic property x^n = -1
proptest! {
    #[test]
    fn negacyclic_property_holds(ctx in valid_context_strategy()) {
        // Create x^(n-1) polynomial: [0, 0, 0, 1]
        let x_to_n_minus_1 = NttPolynomial::from_coeffs([0, 0, 0, 1], Arc::clone(&ctx));
        // Create x polynomial: [0, 1, 0, 0]
        let x = NttPolynomial::from_coeffs([0, 1, 0, 0], Arc::clone(&ctx));

        // x^(n-1) * x should equal x^n = -1 ≡ q-1 (mod q)
        let result = &x_to_n_minus_1 * &x;
        let expected = [ctx.modulus() - 1, 0, 0, 0]; // -1 in first coefficient

        prop_assert_eq!(result.coeffs(), &expected);
    }
}

// Property: Assignment operators work correctly
proptest! {
    #[test]
    fn assignment_operators_consistent(
        coeffs_a in prop::array::uniform4(1u64..100u64),
        coeffs_b in prop::array::uniform4(1u64..100u64),
        ctx in valid_context_strategy()
    ) {
        let a = NttPolynomial::from_coeffs(
            coeffs_a.map(|c| c % ctx.modulus()),
            Arc::clone(&ctx)
        );
        let b = NttPolynomial::from_coeffs(
            coeffs_b.map(|c| c % ctx.modulus()),
            Arc::clone(&ctx)
        );

        // Test AddAssign
        let expected_add = &a + &b;
        let mut a_copy = a.clone();
        a_copy += &b;
        prop_assert_eq!(a_copy.coeffs(), expected_add.coeffs());

        // Test SubAssign
        let expected_sub = &a - &b;
        let mut a_copy = a.clone();
        a_copy -= &b;
        prop_assert_eq!(a_copy.coeffs(), expected_sub.coeffs());

        // Test MulAssign
        let expected_mul = &a * &b;
        let mut a_copy = a.clone();
        a_copy *= &b;
        prop_assert_eq!(a_copy.coeffs(), expected_mul.coeffs());
    }
}
