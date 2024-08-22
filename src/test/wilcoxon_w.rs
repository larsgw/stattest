use core::fmt::Debug;
use std::ops::Sub;

use crate::distribution::SignedRank;
use crate::statistics::*;
use crate::traits::abs::Abs;
use crate::traits::zero::Zero;
use statrs::distribution::ContinuousCDF;

use super::StatisticalTest;

/// Implements the [Wilcoxon signed rank test](https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test).
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct WilcoxonWTest {
    estimate: (f64, f64),
    effect_size: f64,
    p_value: f64,
}

impl WilcoxonWTest {
    /// Run Wilcoxon signed rank test on samples `x` and `y`.
    pub fn paired<I, J>(x: I, y: J) -> statrs::Result<WilcoxonWTest>
    where
        I: IntoIterator,
        J: IntoIterator<Item = I::Item>,
        I::IntoIter: ExactSizeIterator,
        J::IntoIter: ExactSizeIterator,
        I::Item: Copy + Debug + Sub<I::Item>,
        <I::Item as Sub<I::Item>>::Output:
            Abs<Output = <I::Item as Sub<I::Item>>::Output> + PartialOrd + Zero + Copy + Debug,
    {
        let mut x_iter = x.into_iter();
        let mut y_iter = y.into_iter();
        let x_len: usize = x_iter.len();
        let y_len = y_iter.len();

        assert_eq!(x_len, y_len, "Samples must have the same length");

        let mut deltas: Vec<<I::Item as Sub<I::Item>>::Output> =
            x_iter.zip(y_iter).map(|(x, y)| x - y).collect();
        deltas.sort_unstable_by(|a, b| {
            a.abs()
                .partial_cmp(&b.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut tie_solver = ResolveTies::new(
            deltas.iter().copied(),
            <<I::Item as Sub<I::Item>>::Output as Abs>::abs,
        );

        let mut estimate = (0.0, 0.0);
        let mut zeroes = 0;

        for (rank, delta) in &mut tie_solver {
            if delta < <<I::Item as Sub<I::Item>>::Output as Zero>::ZERO {
                estimate.0 += rank;
            } else if delta > <<I::Item as Sub<I::Item>>::Output as Zero>::ZERO {
                estimate.1 += rank;
            } else {
                zeroes += 1;
            }
        }

        let estimate_small = if estimate.0 < estimate.1 {
            estimate.0
        } else {
            estimate.1
        };
        let distribution = SignedRank::new(x_len, zeroes, tie_solver.tie_correction())?;
        let p_value = distribution.cdf(estimate_small);

        let n = x_len as f64;
        let rank_sum = n * (n + 1.0) / 2.0;
        let effect_size = estimate_small / rank_sum;

        Ok(WilcoxonWTest {
            effect_size,
            estimate,
            p_value,
        })
    }
}

impl StatisticalTest for WilcoxonWTest {
    type Estimate = (f64, f64);

    fn estimate(&self) -> (f64, f64) {
        self.estimate
    }

    fn p_value(&self) -> f64 {
        self.p_value
    }

    fn effect_size(&self) -> f64 {
        self.effect_size
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paired() {
        let x: Vec<f64> = vec![8.0, 6.0, 5.5, 11.0, 8.5, 5.0, 6.0, 6.0];
        let y: Vec<f64> = vec![8.5, 9.0, 6.5, 10.5, 9.0, 7.0, 6.5, 7.0];
        let test = WilcoxonWTest::paired(&x, &y).unwrap();
        assert_eq!(test.estimate(), (33.5, 2.5));
        assert_eq!(test.p_value(), 0.027785782704095215);
        assert_eq!(test.effect_size(), 0.06944444444444445);
    }

    #[test]
    fn paired_2() {
        let x: Vec<f64> = vec![209.0, 200.0, 177.0, 169.0, 159.0, 169.0, 187.0, 198.0];
        let y: Vec<f64> = vec![151.0, 168.0, 147.0, 164.0, 166.0, 163.0, 176.0, 188.0];
        let test = WilcoxonWTest::paired(&x, &y).unwrap();
        assert_eq!(test.estimate(), (3.0, 33.0));
        assert_eq!(test.p_value(), 0.0390625);
        assert_eq!(test.effect_size(), 0.08333333333333333);
    }
}
