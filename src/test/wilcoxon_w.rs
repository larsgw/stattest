use crate::distribution::SignedRank;
use crate::statistics::*;
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
    pub fn paired(x: &[f64], y: &[f64]) -> statrs::Result<WilcoxonWTest> {
        let d: Vec<_> = x.iter().zip(y).map(|(x, y)| (x - y).abs()).collect();
        let (ranks, tie_correction) = (&d).ranks();
        let mut estimate = (0.0, 0.0);
        let mut zeroes = 0;

        for ((x, y), rank) in x.iter().zip(y).zip(ranks) {
            if x < y {
                estimate.0 += rank;
            } else if x > y {
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
        let distribution = SignedRank::new(d.len(), zeroes, tie_correction)?;
        let p_value = distribution.cdf(estimate_small);

        let n = (&d).n();
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
        let x = vec![8.0, 6.0, 5.5, 11.0, 8.5, 5.0, 6.0, 6.0];
        let y = vec![8.5, 9.0, 6.5, 10.5, 9.0, 7.0, 6.5, 7.0];
        let test = WilcoxonWTest::paired(&x, &y).unwrap();
        assert_eq!(test.estimate(), (33.5, 2.5));
        assert_eq!(test.p_value(), 0.027785782704095215);
        assert_eq!(test.effect_size(), 0.06944444444444445);
    }

    #[test]
    fn paired_2() {
        let x = vec![209.0, 200.0, 177.0, 169.0, 159.0, 169.0, 187.0, 198.0];
        let y = vec![151.0, 168.0, 147.0, 164.0, 166.0, 163.0, 176.0, 188.0];
        let test = WilcoxonWTest::paired(&x, &y).unwrap();
        assert_eq!(test.estimate(), (3.0, 33.0));
        assert_eq!(test.p_value(), 0.0390625);
        assert_eq!(test.effect_size(), 0.08333333333333333);
    }
}
