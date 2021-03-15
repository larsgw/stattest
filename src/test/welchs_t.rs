use statrs::distribution::{StudentsT, ContinuousCDF};
use statrs::statistics::Statistics;
use crate::statistics::StatisticsExt;

/// Implements [Welch's t-test](https://en.wikipedia.org/wiki/Welch's_t-test) (Welch, 1947).
///
/// # References
///
/// Welch, B. L. (1947). The generalisation of student’s problems when several different population
///     variances are involved. Biometrika, 34(1–2), 28–35. <https://doi.org/10.1093/biomet/34.1-2.28>
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct WelchsTTest  {
    df: f64,
    estimate: f64,
    effect_size: f64,
    p_value: f64
}

impl WelchsTTest  {
    /// Run Welch's two-sample t-test on samples `x` and `y`.
    pub fn new (x: &Vec<f64>, y: &Vec<f64>) -> statrs::Result<WelchsTTest > {
        let var_x = x.variance();
        let var_y = y.variance();
        let var_x_n = var_x / x.n();
        let var_y_n = var_y / y.n();
        let linear_combination = var_x_n + var_y_n;

        let df = linear_combination.powi(2) / (
            var_x_n.powi(2) / x.df() +
            var_y_n.powi(2) / y.df()
        );

        let mean_difference = x.mean() - y.mean();
        let effect_size = mean_difference.abs() / ((var_x + var_y) / 2.0).sqrt();
        let t = mean_difference / linear_combination.sqrt();

        let t_distribution = StudentsT::new(0.0, 1.0, df)?;
        let p_value = 2.0 * t_distribution.cdf(-t.abs());

        Ok(WelchsTTest  {
            df,
            effect_size,
            estimate: t,
            p_value: p_value
        })
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn students_t() {
        let x = vec!(134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0);
        let y = vec!(70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0);
        let test = super::WelchsTTest ::new(&x, &y).unwrap();
        assert_eq!(test.df, 13.081702113268564);
        assert_eq!(test.estimate, 1.9107001042454415);
        assert_eq!(test.effect_size, 0.904358069450997);
        assert_eq!(test.p_value, 0.0782070409214568);
    }

    #[test]
    fn reverse() {
        let x = vec!(134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0);
        let y = vec!(70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0);
        let test = super::WelchsTTest ::new(&y, &x).unwrap();
        assert_eq!(test.df, 13.081702113268564);
        assert_eq!(test.estimate, -1.9107001042454415);
        assert_eq!(test.effect_size, 0.904358069450997);
        assert_eq!(test.p_value, 0.0782070409214568);
    }
}
