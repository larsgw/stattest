use statrs::statistics::Statistics;
use crate::statistics::StatisticsExt;
use std::borrow::Borrow;
use std::f64;

impl<T> StatisticsExt<f64> for T
where
    T: IntoIterator + Clone,
    T::Item: Borrow<f64>,
{
    fn n(self) -> f64 {
        self.into_iter().count() as f64
    }

    fn df(self) -> f64 {
        self.n() - 1.0
    }

    fn pooled_variance(self, other: Self) -> f64 {
        let df_x = self.clone().df();
        let df_y = other.clone().df();
        let var_x = self.clone().variance();
        let var_y = other.clone().variance();

        (df_x * var_x + df_y * var_y) / (df_x + df_y)
    }

    fn pooled_std_dev(self, other: Self) -> f64 {
        self.pooled_variance(other).sqrt()
    }

    fn variance_ratio(self, other: Self) -> f64 {
        self.variance() / other.variance()
    }
}

#[cfg(test)]
mod tests {
    #[allow(dead_code)]
    fn round(number: f64, digits: Option<i32>) -> f64 {
        match digits {
            Some(digits) => {
                let factor = 10_f64.powi(digits);
                (number * factor).round() / factor
            }
            None => number.round()
        }
    }

    #[test]
    fn pooled_variance() {
        let x = vec!(134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0);
        let y = vec!(70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0);
        assert_eq!(round(super::StatisticsExt::pooled_variance(&x, &y), Some(3)), 446.118);
    }

    #[test]
    fn pooled_std_dev() {
        let x = vec!(134.0, 146.0, 104.0, 119.0, 124.0, 161.0, 107.0, 83.0, 113.0, 129.0, 97.0, 123.0);
        let y = vec!(70.0, 118.0, 101.0, 85.0, 107.0, 132.0, 94.0);
        assert_eq!(round(super::StatisticsExt::pooled_std_dev(&x, &y), Some(3)), 21.121);
    }
}
