/// Provides additional statistical utilities, complementing [statrs::statistics::Statistics].
pub trait StatisticsExt<T> {
    /// Returns the amount of observations.
    fn n(self) -> T;
    /// Returns the degrees of freedom of the data.
    fn df(self) -> T;

    /// Returns the pooled variance.
    fn pooled_variance(self, other: Self) -> T;
    /// Returns the pooled standard deviation.
    fn pooled_std_dev(self, other: Self) -> T;

    /// Returns the ratio between two variances.
    fn variance_ratio(self, other: Self) -> T;
}
