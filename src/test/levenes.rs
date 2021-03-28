use crate::statistics::StatisticsExt;
use statrs::distribution::{ContinuousCDF, FisherSnedecor};
use statrs::statistics::Statistics;

/// Implements [Levene's test](https://en.wikipedia.org/wiki/Levene%27s_test) (Brown & Forsythe, 1974).
///
/// # References
///
/// Brown, M. B., & Forsythe, A. B. (1974). Robust Tests for the Equality of Variances.
///     Journal of the American Statistical Association, 69(346), 364–367. <https://doi.org/10.2307/2285659>
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct LevenesTest {
    df: f64,
    estimate: f64,
    p_value: f64,
}

impl LevenesTest {
    /// Run Levene's test on the samples `x` and `y`.
    pub fn new(x: &[f64], y: &[f64]) -> statrs::Result<LevenesTest> {
        let n_x = x.n();
        let n_y = y.n();
        let diff_x = x.iter().map(|xi| (xi - x.mean()).abs());
        let diff_y = y.iter().map(|yi| (yi - y.mean()).abs());

        let mean_diff_x = diff_x.clone().mean();
        let mean_diff_y = diff_y.clone().mean();
        let mean_diff = Iterator::chain(diff_x.clone(), diff_y.clone()).mean();

        let a: f64 =
            n_x * (mean_diff_x - mean_diff).powi(2) + n_y * (mean_diff_y - mean_diff).powi(2);
        let b: f64 = Iterator::chain(
            diff_x.map(|diff| (diff - mean_diff_x).powi(2)),
            diff_y.map(|diff| (diff - mean_diff_y).powi(2)),
        )
        .sum();

        let df = n_x + n_y - 2.0;
        let estimate = df * a / b;
        let distribution = FisherSnedecor::new(1.0, df)?;
        let p_value = 1.0 - distribution.cdf(estimate);

        Ok(LevenesTest {
            df,
            estimate,
            p_value,
        })
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn levenes_test() {
        let x = vec![
            134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0,
        ];
        let y = vec![70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0];
        let result = super::LevenesTest::new(&x, &y).unwrap();
        assert_eq!(result.df, 17.0);
        assert_eq!(result.estimate, 0.014721055064513417);
        assert_eq!(result.p_value, 0.9048519802923365);
    }
}
