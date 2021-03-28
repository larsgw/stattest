use crate::statistics::StatisticsExt;
use statrs::distribution::{ContinuousCDF, FisherSnedecor};

/// Implements the [F-test of equality of variances](https://en.wikipedia.org/wiki/F-test_of_equality_of_variances).
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct FTest {
    df: (f64, f64),
    estimate: f64,
    p_value: f64,
}

impl FTest {
    /// Carry out the F-test of equality of variances on the samples `x` and `y`.
    pub fn new(x: &[f64], y: &[f64]) -> statrs::Result<FTest> {
        let f = x.variance_ratio(y);
        let df = (x.df(), y.df());

        let distribution = FisherSnedecor::new(df.0, df.1)?;
        let probability = distribution.cdf(f);
        let p_value = if f.gt(&1.0) {
            1.0 - probability
        } else {
            probability
        };

        Ok(FTest {
            df,
            estimate: f,
            p_value,
        })
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn f_test() {
        let x = vec![
            134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0,
        ];
        let y = vec![70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0];
        let result = super::FTest::new(&x, &y).unwrap();
        assert_eq!(result.df, (11.0, 6.0));
        assert_eq!(result.estimate, 1.0755200911940725);
        assert_eq!(result.p_value, 0.4893961256182331);
    }

    #[test]
    fn f_test_reverse() {
        let x = vec![
            134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0,
        ];
        let y = vec![70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0];
        let result = super::FTest::new(&y, &x).unwrap();
        assert_eq!(result.df, (6.0, 11.0));
        assert_eq!(result.estimate, 0.9297827239003709);
        assert_eq!(result.p_value, 0.48939612561823265);
    }
}
