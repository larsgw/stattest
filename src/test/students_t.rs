use statrs::statistics::Statistics;
use statrs::distribution::{StudentsT, ContinuousCDF};
use crate::statistics::StatisticsExt;

/// Implements [Student's t-test](https://en.wikipedia.org/wiki/Student%27s_t-test).
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct StudentsTTest {
    df: f64,
    estimate: f64,
    p_value: f64
}

impl StudentsTTest {
    /// Run Student's two-sample t-test on samples `x` and `y`.
    pub fn new (x: &Vec<f64>, y: &Vec<f64>) -> statrs::Result<StudentsTTest> {
        let n_x = x.n();
        let n_y = y.n();
        let df = n_x + n_y;

        let var_p = x.pooled_variance(y);

        let t = (x.mean() - y.mean()) / (var_p * (n_x.recip() + n_y.recip())).sqrt();

        let t_distribution = StudentsT::new(0.0, 1.0, df)?;
        let p_value = 2.0 * t_distribution.cdf(-t.abs());

        Ok(StudentsTTest {
            df,
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
        let test = super::StudentsTTest::new(&x, &y).unwrap();
        assert_eq!(test.estimate, 1.89143639744233);
        assert_eq!(test.p_value, 0.073911127032672);
    }

    #[test]
    fn reverse() {
        let x = vec!(134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0);
        let y = vec!(70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0);
        let test = super::StudentsTTest::new(&y, &x).unwrap();
        assert_eq!(test.estimate, -1.89143639744233);
        assert_eq!(test.p_value, 0.073911127032672);
    }
}
