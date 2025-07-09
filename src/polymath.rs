// Number Theoretic operations over polynomials and vectors
use crate::ntmath::*;
use rand::{Rng, rng};

#[derive(Debug, Clone, Copy)]
pub struct PolyRing {
    pub class: CongruenceClass,
    pub n: usize,
}

impl PolyRing {
    pub fn new(q: u64, n: usize) -> Self {
        let class = CongruenceClass::new(q);
        Self { class, n }
    }

    pub fn sample_random(&self) -> Vec<u64> {
        let mut generator = rng();
        (0..self.n)
            .map(|_| generator.random_range(1..self.class.q))
            .collect()
    }

    pub fn add(&self, ax: &[u64], bx: &[u64]) -> Vec<u64> {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);
        assert_eq!(bx.len(), self.n, "input must have len = {}", self.n);

        ax.iter()
            .zip(bx)
            .map(|(&a, &b)| self.class.modadd(a, b))
            .collect()
    }

    pub fn add_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);
        assert_eq!(bx.len(), self.n, "input must have len = {}", self.n);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modadd_eq(&mut *a, b);
        }
    }

    pub fn sub(&self, ax: &[u64], bx: &[u64]) -> Vec<u64> {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);
        assert_eq!(bx.len(), self.n, "input must have len = {}", self.n);

        ax.iter()
            .zip(bx)
            .map(|(&a, &b)| self.class.modsub(a, b))
            .collect()
    }

    pub fn sub_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);
        assert_eq!(bx.len(), self.n, "input must have len = {}", self.n);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modsub_eq(&mut *a, b);
        }
    }

    pub fn neg(&self, ax: &[u64]) -> Vec<u64> {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);

        ax.iter().map(|&a| self.class.modneg(a)).collect()
    }

    pub fn neg_eq(&self, ax: &mut [u64]) {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);

        for a in ax.iter_mut() {
            self.class.modneg_eq(&mut *a);
        }
    }

    pub fn mul(&self, ax: &[u64], bx: &[u64]) -> Vec<u64> {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);
        assert_eq!(bx.len(), self.n, "input must have len = {}", self.n);

        ax.iter()
            .zip(bx)
            .map(|(&a, &b)| self.class.modmul(a, b))
            .collect()
    }

    pub fn mul_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);
        assert_eq!(bx.len(), self.n, "input must have len = {}", self.n);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modmul_eq(&mut *a, b);
        }
    }

    pub fn naive_negacyclic_convolution(&self, ax: &[u64], bx: &[u64]) -> Vec<u64> {
        assert_eq!(ax.len(), self.n, "input must have len = {}", self.n);
        assert_eq!(bx.len(), self.n, "input must have len = {}", self.n);

        let mut cx = Vec::with_capacity(self.n + 1);

        for i in 0..self.n {
            cx.push(0);

            for j in 0..=i {
                self.class
                    .modadd_eq(&mut cx[i], self.class.modmul(ax[j], bx[i - j]));
            }
            for j in (i + 1)..self.n {
                self.class
                    .modsub_eq(&mut cx[i], self.class.modmul(ax[j], bx[self.n + i - j]));
            }
        }

        cx
    }
}
