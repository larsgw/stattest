use statrs::distribution::{Normal, ContinuousCDF};
use crate::statistics::*;

/// Implements the [Wilcoxon rank-sum test](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test).
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct WilcoxonUTest {
    estimate: (f64, f64),
    p_value: f64
}

impl WilcoxonUTest {
    /// Run Wilcoxon rank-sum test on samples `x` and `y`.
    pub fn new (x: &Vec<f64>, y: &Vec<f64>) -> statrs::Result<WilcoxonUTest> {
        let (ranks, tie_correction) = x.iter().chain(y).ranks();
        let n_x = x.n();
        let n_y = y.n();
        let n_xy = n_x * n_y;

        let u = (n_x * (n_x + 1.0)) / 2.0 - ranks[0..x.len()].iter().sum::<f64>();
        let u_x = n_xy + u;
        let u_y = -u;

        let m = n_xy / 2.0;
        let s = ((n_xy * (n_x + n_y + 1.0 - tie_correction)) / 12.0).sqrt();

        let normal = Normal::new(m, s)?;
        let p_value = 2.0 * normal.cdf(if u < 0.0 { u_x } else { u_y });

        Ok(WilcoxonUTest {
            estimate: (u_x, u_y),
            p_value: p_value
        })
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn wilcoxon_u() {
        let x = vec!(134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0);
        let y = vec!(70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0);
        let test = super::WilcoxonUTest::new(&x, &y).unwrap();
        assert_eq!(test.p_value, 0.08303763193135497);
    }
}
