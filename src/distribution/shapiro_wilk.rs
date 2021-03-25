use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use statrs::function::evaluate::polynomial;
use statrs::statistics::*;
use statrs::Result;
use rand::Rng;

// Polynomials used to calculate the mean and standard deviation
// of a normal distribution that approximates the distribution of
// the W statistic of the Shapiro-Wilk test (Royston, 1992).
static C3: [f64; 4] = [ 0.5440, -0.39978 ,  0.025054 , -0.0006714];
static C4: [f64; 4] = [ 1.3822, -0.77857 ,  0.062767 , -0.0020322];
static C5: [f64; 4] = [-1.5861, -0.31082 , -0.083751 ,  0.0038915];
static C6: [f64; 3] = [-0.4803, -0.082676,  0.0030302];

/// Implements an approximation of the distribution of the W
/// statistic of the [Shapiro-Wilk test](https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test).
///
/// # References
///
/// Royston, P. (1992). Approximating the Shapiro-Wilk W-test for non-normality.
///     Statistics and Computing, 2(3), 117–119. <https://doi.org/10.1007/BF01891203>
///
/// Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality
///      (complete samples)†. Biometrika, 52(3–4), 591–611. <https://doi.org/10.1093/biomet/52.3-4.591>
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct ShapiroWilk {
    normal: Normal,
    n: usize,
    mean: f64,
    std_dev: f64
}

impl ShapiroWilk {
    /// Create a new approximation for the distribution of W for a given sample size `n`.
    ///
    /// # Examples
    ///
    /// ```
    /// use stattest::distribution::ShapiroWilk;
    ///
    /// let result = ShapiroWilk::new(25);
    /// assert!(result.is_ok());
    /// ```
    pub fn new(n: usize) -> Result<ShapiroWilk> {
        let float_n = n as f64;

        let (mean, std_dev) = if n <= 11 {
            (polynomial(float_n, &C3), polynomial(float_n, &C4).exp())
        } else {
            let log_n = float_n.ln();
            (polynomial(log_n, &C5), polynomial(log_n, &C6).exp())
        };

        let result = Normal::new(mean, std_dev);
        match result {
            Ok(normal) => Ok(ShapiroWilk { normal, n, mean, std_dev }),
            Err(error) => Err(error)
        }
    }

    /// Return the sample size associated with this distribution.
    pub fn n(&self) -> usize {
        self.n
    }
}

impl ::rand::distributions::Distribution<f64> for ShapiroWilk {
    fn sample<R: Rng + ?Sized>(&self, r: &mut R) -> f64 {
        ::rand::distributions::Distribution::sample(&self.normal, r)
    }
}

impl ContinuousCDF<f64, f64> for ShapiroWilk {
    fn cdf(&self, x: f64) -> f64 {
        self.normal.cdf(x)
    }

    fn inverse_cdf(&self, x: f64) -> f64 {
        self.normal.inverse_cdf(x)
    }
}

impl Min<f64> for ShapiroWilk {
    fn min(&self) -> f64 {
        self.normal.min()
    }
}

impl Max<f64> for ShapiroWilk {
    fn max(&self) -> f64 {
        self.normal.max()
    }
}

impl Distribution<f64> for ShapiroWilk {
    fn mean(&self) -> Option<f64> {
        self.normal.mean()
    }
    fn variance(&self) -> Option<f64> {
        self.normal.variance()
    }
    fn entropy(&self) -> Option<f64> {
        self.normal.entropy()
    }
    fn skewness(&self) -> Option<f64> {
        self.normal.skewness()
    }
}

impl Mode<Option<f64>> for ShapiroWilk {
    fn mode(&self) -> Option<f64> {
        self.normal.mode()
    }
}

impl Continuous<f64, f64> for ShapiroWilk {
    fn pdf(&self, x: f64) -> f64 {
        self.normal.pdf(x)
    }

    fn ln_pdf(&self, x: f64) -> f64 {
        self.normal.ln_pdf(x)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn n_25() {
        let distribution = super::ShapiroWilk::new(25).unwrap();
        assert_eq!(distribution.mean, -3.3245620722111475);
        assert_eq!(distribution.std_dev, 0.4891787150186814);
    }
}
