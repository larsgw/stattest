use crate::statistics::*;
use statrs::distribution::{ContinuousCDF, Normal};

/// Implements the [Mann-Whitney U test](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test),
/// also known as the Wilcoxon rank-sum test.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct MannWhitneyUTest {
    estimate: (f64, f64),
    effect_size: f64,
    p_value: f64,
}

impl MannWhitneyUTest {
    /// Run Mann-Whitney U test/Wilcoxon rank-sum test on samples `x` and `y`.
    pub fn independent(x: &[f64], y: &[f64]) -> statrs::Result<MannWhitneyUTest> {
        let (ranks, tie_correction) = x.iter().chain(y).ranks();
        let n_x = x.n();
        let n_y = y.n();
        let n_xy = n_x * n_y;

        let estimate = (n_x * (n_x + 1.0)) / 2.0 - ranks[0..x.len()].iter().sum::<f64>();
        let estimate_x = n_xy + estimate;
        let estimate_y = -estimate;
        let estimate_small = if estimate < 0.0 {
            estimate_x
        } else {
            estimate_y
        };

        let n = n_x + n_y;
        let distribution_mean = n_xy / 2.0;
        let distribution_var = (n_xy * (n + 1.0 - tie_correction as f64 / (n * (n - 1.0)))) / 12.0;

        let normal = Normal::new(distribution_mean, distribution_var.sqrt())?;
        let p_value = 2.0 * normal.cdf(estimate_small);
        let effect_size = 1.0 - (2.0 * estimate_small) / n_xy;

        Ok(MannWhitneyUTest {
            effect_size,
            estimate: (estimate_x, estimate_y),
            p_value,
        })
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn mann_whitney_u() {
        let x = vec![
            134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0,
        ];
        let y = vec![70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0];
        let test = super::MannWhitneyUTest::independent(&x, &y).unwrap();
        assert_eq!(test.estimate, (21.5, 62.5));
        assert_eq!(test.effect_size, 0.48809523809523814);
        assert_eq!(test.p_value, 0.08303763193135497);
    }

    #[test]
    fn mann_whitney_u_2() {
        let x = vec![68.0, 68.0, 59.0, 72.0, 64.0, 67.0, 70.0, 74.0];
        let y = vec![60.0, 67.0, 61.0, 62.0, 67.0, 63.0, 56.0, 58.0];
        let test = super::MannWhitneyUTest::independent(&x, &y).unwrap();
        assert_eq!(test.estimate, (9.0, 55.0));
        assert_eq!(test.effect_size, 0.71875);
        assert_eq!(test.p_value, 0.01533316211294691);
    }
}
