use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use statrs::function::factorial::binomial;
use statrs::statistics::*;
use statrs::Result;
use rand::Rng;

/// Implements an approximation of the distribution of the W
/// statistic of the [Wilcoxon signed rank test](https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test).
///
/// # References
///
/// Wilcoxon, F. (1945). Individual Comparisons by Ranking Methods.
///     Biometrics Bulletin, 1(6), 80–83. <https://doi.org/10.2307/3001968>
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct SignedRank {
    approximation: Approximation,
    n: usize
}

#[derive(Debug, Copy, Clone, PartialEq)]
enum Approximation {
    Normal(Normal, f64),
    Exact
}

impl SignedRank {
    /// Determine whether to use an approximation or an exact distribution of W for a given sample size `n`
    /// and given a number of paired observations with a non-zero difference `non_zero`.
    ///
    /// # Examples
    ///
    /// ```
    /// use stattest::distribution::SignedRank;
    ///
    /// let result = SignedRank::new(25);
    /// assert!(result.is_ok());
    /// ```
    pub fn new(n: usize, non_zero: usize) -> Result<SignedRank> {
        if non_zero >= 20 {
            SignedRank::approximate(n, non_zero)
        } else {
            SignedRank::exact(n)
        }
    }

    /// Create a new approximation for the distribution of W for a given sample size `n`
    /// and given a number of paired observations with a non-zero difference `non_zero`.
    ///
    /// # Examples
    ///
    /// ```
    /// use stattest::distribution::SignedRank;
    ///
    /// let result = SignedRank::new(25);
    /// assert!(result.is_ok());
    /// ```
    pub fn approximate(n: usize, non_zero: usize) -> Result<SignedRank> {
        let mean = (non_zero * (non_zero + 1)) as f64 / 4.0;
        let var = mean * (2 * non_zero + 1) as f64 / 6.0;
        let normal = Normal::new(mean, var.sqrt())?;

        let rank_sum = (n * (n + 1)) as f64 / 2.0;
        let scale = normal.cdf(rank_sum) - normal.cdf(0.0);

        Ok(SignedRank {
            n,
            approximation: Approximation::Normal(normal, scale)
        })
    }

    /// Create a new exact distribution of W for a given sample size `n`.
    ///
    /// # Examples
    ///
    /// ```
    /// use stattest::distribution::SignedRank;
    ///
    /// let result = SignedRank::new(25);
    /// assert!(result.is_ok());
    /// ```
    pub fn exact(n: usize) -> Result<SignedRank> {
        Ok(SignedRank {
            n,
            approximation: Approximation::Exact
        })
    }
}

fn partitions(number: usize, parts: usize, size: usize) -> usize {
    if number == parts {
        1
    } else if number == 0 || parts == 0 || size == 0 {
        0
    } else {
        (1..=std::cmp::min(number, size))
            .map(|highest_part| partitions(number - highest_part, parts - 1, highest_part))
            .sum()
    }
}

impl ::rand::distributions::Distribution<f64> for SignedRank {
    fn sample<R: Rng + ?Sized>(&self, r: &mut R) -> f64 {
        match self.approximation {
            Approximation::Normal(normal, _) => ::rand::distributions::Distribution::sample(&normal, r),
            Approximation::Exact => {
                todo!();
            }
        }
    }
}

impl ContinuousCDF<f64, f64> for SignedRank {
    fn cdf(&self, x: f64) -> f64 {
        match self.approximation {
            Approximation::Normal(normal, scale) => 2.0 * normal.cdf(x) / scale,
            Approximation::Exact => {
                let r = x.ceil() as usize;

                let mut sum = 1;

                for n in 1..=self.n {
                    let n_choose_2 = binomial(n as u64, 2) as usize;
                    if n_choose_2 > r { continue }
                    for i in n..=(r - n_choose_2) {
                        sum += partitions(i, n, self.n - n + 1);
                    }
                }

                sum as f64 / 2_f64.powi(self.n as i32 - 1)
            }
        }
    }

    fn inverse_cdf(&self, x: f64) -> f64 {
        match self.approximation {
            Approximation::Normal(normal, _) => normal.inverse_cdf(x),
            Approximation::Exact => {
                todo!();
            }
        }
    }
}

/*
// Solve 1 + 2 + ... + k < r for k
// 1 + 2 + ... + k = k(k + 1) / 2
// k(k + 1) / 2 < r
// k² + k - 2r < 0
// let end = (0.5 * (8.0 * x + 1.0).sqrt() + 0.5) as usize;
// let end = (0.5 * (4.0 * x + 1.0).sqrt() + 0.5) as usize;

//     // Loop over all totals smaller than or equal to r
//     for n in 1..=r {
//         // For each total, we want the number of ways in which it can be constructed
//         // from the ranks, i.e. the number of unequal partitions of n. The number of
//         // partitions of r with no part larger than j is equal to the number of j-part
//         // partitions of r - (j choose 2).
//
//         let mut tmp = 0;
//         for j in 1..=self.n {
//             // if j > n {
//             //     break
//             // }
//             let j_choose_2 = binomial(j as u64, 2) as usize;
//             if j_choose_2 > n {
//                 break
//             }
//             println!("total={} {}-part partitions of {}: {}", r, j, n, partitions(n - j_choose_2, j));
//             tmp += partitions(n - j_choose_2, j);
//         }
//         sum += tmp;
//         // println!("{} {} {}", r, n, tmp);
//
//         // let n_choose_2 = binomial(std::cmp::min(n, self.n) as u64, 2) as usize;
//         // if n_choose_2 > r {
//         //     break
//         // }
//         // // sum += partitions(r - n_choose_2, n);
//         // // println!("{} {}", self.n, r - n_choose_2);
//         // for i in n..=(r - n_choose_2) {
//         //     // println!("{} {} {}", i, n, partitions(i, n));
//         //     sum += partitions(i, n);
//         // }
//     }
//
// sum as f64
*/

impl Min<f64> for SignedRank {
    fn min(&self) -> f64 {
        0.0
    }
}

impl Max<f64> for SignedRank {
    fn max(&self) -> f64 {
        (self.n * (self.n + 1) / 2) as f64
    }
}

impl Distribution<f64> for SignedRank {
    fn mean(&self) -> Option<f64> {
        match self.approximation {
            Approximation::Normal(normal, _) => normal.mean(),
            Approximation::Exact => {
                todo!();
            }
        }
    }
    fn variance(&self) -> Option<f64> {
        match self.approximation {
            Approximation::Normal(normal, _) => normal.variance(),
            Approximation::Exact => {
                todo!();
            }
        }
    }
    fn entropy(&self) -> Option<f64> {
        match self.approximation {
            Approximation::Normal(normal, _) => normal.entropy(),
            Approximation::Exact => {
                todo!();
            }
        }
    }
    fn skewness(&self) -> Option<f64> {
        match self.approximation {
            Approximation::Normal(normal, _) => normal.skewness(),
            Approximation::Exact => {
                todo!();
            }
        }
    }
}

impl Mode<Option<f64>> for SignedRank {
    fn mode(&self) -> Option<f64> {
        match self.approximation {
            Approximation::Normal(normal, _) => normal.mode(),
            Approximation::Exact => {
                todo!();
            }
        }
    }
}

impl Continuous<f64, f64> for SignedRank {
    fn pdf(&self, x: f64) -> f64 {
        match self.approximation {
            Approximation::Normal(normal, _) => normal.pdf(x),
            Approximation::Exact => {
                todo!();
            }
        }
    }

    fn ln_pdf(&self, x: f64) -> f64 {
        match self.approximation {
            Approximation::Normal(normal, _) => normal.ln_pdf(x),
            Approximation::Exact => {
                todo!();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use statrs::distribution::ContinuousCDF;
    use statrs::statistics::Distribution;

    #[test]
    fn n_20() {
        let distribution = super::SignedRank::new(20, 20).unwrap();
        assert_eq!(distribution.mean(), Some(105.0));
        assert_eq!(distribution.std_dev(), Some(26.78619047195775));
        assert_eq!(distribution.cdf(50.0), 0.0400473452471253);
    }

    #[test]
    fn n_30() {
        let distribution = super::SignedRank::new(30, 30).unwrap();
        assert_eq!(distribution.mean(), Some(232.5));
        assert_eq!(distribution.std_dev(), Some(48.61841215013094));
        assert_eq!(distribution.cdf(150.0), 0.08971799337794222);
    }

    #[test]
    fn n_20_exact() {
        let distribution = super::SignedRank::exact(20).unwrap();
        assert_eq!(distribution.cdf(50.0), 0.039989471435546875);
    }

    #[test]
    fn n_11() {
        let distribution = super::SignedRank::exact(11).unwrap();
        assert_eq!(distribution.cdf(11.0), 0.0537109375);
        assert_eq!(distribution.cdf(7.0), 0.0185546875);
        assert_eq!(distribution.cdf(5.0), 0.009765625);
    }

    #[test]
    fn n_10() {
        let distribution = super::SignedRank::exact(10).unwrap();
        assert_eq!(distribution.cdf(8.0), 0.048828125);
        assert_eq!(distribution.cdf(5.0), 0.01953125);
        assert_eq!(distribution.cdf(3.0), 0.009765625);
    }

    #[test]
    fn n_9() {
        let distribution = super::SignedRank::exact(9).unwrap();
        assert_eq!(distribution.cdf(6.0), 0.0546875);
        assert_eq!(distribution.cdf(3.0), 0.01953125);
        assert_eq!(distribution.cdf(2.0), 0.01171875);
    }

    #[test]
    fn n_8() {
        let distribution = super::SignedRank::exact(8).unwrap();
        assert_eq!(distribution.cdf(4.0), 0.0546875);
        assert_eq!(distribution.cdf(3.0), 0.0390625);
        // assert_eq!(distribution.cdf(2.5), 0.0390625);
        assert_eq!(distribution.cdf(2.0), 0.0234375);
    }

    #[test]
    fn partition() {
        assert_eq!(super::partitions(7, 3, 5), 4);
    }
}
