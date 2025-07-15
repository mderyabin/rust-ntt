// Number Theoretic operations over polynomials and vectors
use crate::ntmath::*;
use rand::{Rng, rng};

#[derive(Debug, Clone, Copy)]
pub struct PolyRing<const N: usize>  {
    pub class: CongruenceClass,

    //  tf: [u64; N],
    // itf: [u64; N],
}

impl<const N: usize>  PolyRing<N> {
    pub fn new(q: u64) -> Self {
        let class = CongruenceClass::new(q);
        Self { class }
    }

    pub fn sample_random(&self) -> [u64; N] {
        let mut generator = rng();
        let mut res = [0u64; N];
        for i in 0..N {
            // .map(|_| generator.random_range(1..self.class.q))
            // .collect()
            res[i] = generator.random_range(1..self.class.q);
        }
        res
    }

    pub fn add(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        for i in 0..N {
            cx[i] = self.class.modadd(ax[i], bx[i])
        }

        cx
    }

    pub fn add_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modadd_eq(&mut *a, b);
        }
    }

    pub fn sub(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        for i in 0..N {
            cx[i] = self.class.modsub(ax[i], bx[i])
        }

        cx
    }

    pub fn sub_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modsub_eq(&mut *a, b);
        }
    }

    pub fn neg(&self, ax: &[u64]) -> Vec<u64> {
        assert_eq!(ax.len(), N, "input must have len = {}", N);

        ax.iter().map(|&a| self.class.modneg(a)).collect()
    }

    pub fn neg_eq(&self, ax: &mut [u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);

        for a in ax.iter_mut() {
            self.class.modneg_eq(&mut *a);
        }
    }

    pub fn mul(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        for i in 0..N {
            cx[i] = self.class.modmul(ax[i], bx[i])
        }

        cx
    }

    pub fn mul_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modmul_eq(&mut *a, b);
        }
    }

    pub fn naive_negacyclic_convolution(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        for i in 0..N {
            for j in 0..=i {
                self.class
                    .modadd_eq(&mut cx[i], self.class.modmul(ax[j], bx[i - j]));
            }
            for j in (i + 1)..N {
                self.class
                    .modsub_eq(&mut cx[i], self.class.modmul(ax[j], bx[N + i - j]));
            }
        }

        cx
    }
}
