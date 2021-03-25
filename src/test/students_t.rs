use statrs::distribution::{StudentsT, ContinuousCDF};
use statrs::statistics::Statistics;
use crate::statistics::StatisticsExt;

/// Implements [Student's t-test](https://en.wikipedia.org/wiki/Student%27s_t-test).
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct StudentsTTest {
    df: f64,
    estimate: f64,
    effect_size: f64,
    p_value: f64
}

impl StudentsTTest {
    /// Run Student's two-sample t-test on samples `x` and `y`.
    pub fn independent (x: &Vec<f64>, y: &Vec<f64>) -> statrs::Result<StudentsTTest> {
        let n_x = x.n();
        let n_y = y.n();
        let df = n_x + n_y;

        let effect_size = (x.mean() - y.mean()) / x.pooled_variance(y).sqrt();
        let t = effect_size / (n_x.recip() + n_y.recip()).sqrt();

        let t_distribution = StudentsT::new(0.0, 1.0, df)?;
        let p_value = 2.0 * t_distribution.cdf(-t.abs());

        Ok(StudentsTTest {
            df,
            effect_size: effect_size.abs(),
            estimate: t,
            p_value: p_value
        })
    }

    /// Run paired Student's t-test on samples `x` and `y`.
    pub fn paired (x: &Vec<f64>, y: &Vec<f64>) -> statrs::Result<StudentsTTest > {
        let d: Vec<_> = x.iter().zip(y).map(|(x, y)| x - y).collect();
        let df = x.df();
        let effect_size = (&d).mean() / (&d).std_dev();
        let t = effect_size * x.n().sqrt();

        let t_distribution = StudentsT::new(0.0, 1.0, df)?;
        let p_value = 2.0 * t_distribution.cdf(-t.abs());

        Ok(StudentsTTest {
            df,
            effect_size: effect_size.abs(),
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
        let test = super::StudentsTTest::independent(&x, &y).unwrap();
        assert_eq!(test.estimate, 1.8914363974423305);
        assert_eq!(test.effect_size, 0.8995574392432595);
        assert_eq!(test.p_value, 0.073911127032672);
    }

    #[test]
    fn reverse() {
        let x = vec!(134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0);
        let y = vec!(70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0);
        let test = super::StudentsTTest::independent(&y, &x).unwrap();
        assert_eq!(test.estimate, -1.8914363974423305);
        assert_eq!(test.effect_size, 0.8995574392432595);
        assert_eq!(test.p_value, 0.073911127032672);
    }

    #[test]
    fn paired() {
        let x = vec!(8.0, 6.0, 5.5, 11.0, 8.5, 5.0, 6.0, 6.0);
        let y = vec!(8.5, 9.0, 6.5, 10.5, 9.0, 7.0, 6.5, 7.0);
        let test = super::StudentsTTest::paired(&x, &y).unwrap();
        assert_eq!(test.estimate, -2.645751311064591);
        assert_eq!(test.effect_size, 0.9354143466934856);
        assert_eq!(test.p_value, 0.03314550026377362);
    }
}
