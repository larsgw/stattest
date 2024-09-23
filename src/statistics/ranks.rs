use core::cmp::Ordering;
use core::convert::TryFrom;

use crate::traits::abs::Abs;

/// Trait defining how many times a value occurs.
pub trait Occurrence: Copy {
    /// Convert the weight to a floating point number.
    fn to_occurrence(&self) -> usize;
}

impl Occurrence for () {
    fn to_occurrence(&self) -> usize {
        1
    }
}

impl Occurrence for u8 {
    fn to_occurrence(&self) -> usize {
        usize::from(*self)
    }
}

impl Occurrence for &u8 {
    fn to_occurrence(&self) -> usize {
        usize::from(**self)
    }
}

impl Occurrence for u16 {
    fn to_occurrence(&self) -> usize {
        usize::from(*self)
    }
}

impl Occurrence for &u16 {
    fn to_occurrence(&self) -> usize {
        usize::from(**self)
    }
}

impl Occurrence for u32 {
    fn to_occurrence(&self) -> usize {
        usize::try_from(*self).unwrap()
    }
}

impl Occurrence for &u32 {
    fn to_occurrence(&self) -> usize {
        usize::try_from(**self).unwrap()
    }
}

impl Occurrence for usize {
    fn to_occurrence(&self) -> usize {
        *self
    }
}

impl Occurrence for &usize {
    fn to_occurrence(&self) -> usize {
        **self
    }
}

#[derive(Debug, Copy, Clone)]
/// A tuple of a value and its occurrences.
pub(crate) struct WeightedTuple<T, O> {
    /// The value.
    pub(crate) value: T,
    /// The occurrences.
    pub(crate) occurrences: O,
}

impl<T, W> From<(T, W)> for WeightedTuple<T, W> {
    fn from((value, occurrences): (T, W)) -> Self {
        WeightedTuple { value, occurrences }
    }
}

impl<T, W> PartialEq for WeightedTuple<T, W>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<T, W> PartialOrd for WeightedTuple<T, W>
where
    T: PartialOrd,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.value.partial_cmp(&other.value)
    }
}

impl<T, W> Abs for WeightedTuple<T, W>
where
    T: Abs,
{
    type Output = WeightedTuple<T::Output, W>;
    fn abs(self) -> Self::Output {
        WeightedTuple {
            value: self.value.abs(),
            occurrences: self.occurrences,
        }
    }
}

impl<T: Copy, W> Occurrence for WeightedTuple<T, W>
where
    W: Occurrence,
{
    fn to_occurrence(&self) -> usize {
        self.occurrences.to_occurrence()
    }
}

/// For the ranking of groups of variables.
pub trait Ranks<T> {
    /// Returns a vector of ranks.
    fn ranks(self) -> (Vec<T>, usize);
}

impl<T> Ranks<f64> for T
where
    T: IntoIterator,
    T::Item: PartialOrd + Copy + Default,
{
    #[inline]
    fn ranks(self) -> (Vec<f64>, usize) {
        let mut observations: Vec<(usize, T::Item)> = self.into_iter().enumerate().collect();
        observations.sort_unstable_by(|(_, a), (_, b)| {
            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut resolved_ties =
            ResolveTies::from(observations.iter().map(|(_, value)| WeightedTuple::from((*value, ()))));
        let mut ranks = vec![0.0; observations.len()];

        for ((rank, _), old_index) in
            (&mut resolved_ties).zip(observations.iter().map(|(index, _)| *index))
        {
            ranks[old_index] = rank;
        }

        (ranks, resolved_ties.tie_correction())
    }
}

pub(crate) struct ResolveTies<I, F>
where
    I: Iterator,
{
    iter: I,
    index: usize,
    resolved: f64,
    tie_correction: usize,
    current_normalized_item: Option<I::Item>,
    normalize: F,
}

impl<I, F> ResolveTies<I, F>
where
    I: Iterator,
{
    #[inline]
    pub(crate) fn new(iter: I, normalize: F) -> Self {
        ResolveTies {
            iter,
            index: 0,
            resolved: 0.0,
            current_normalized_item: None,
            tie_correction: 0,
            normalize,
        }
    }

    #[inline]
    pub(crate) fn tie_correction(&self) -> usize {
        self.tie_correction
    }
}

impl<I, J> From<J> for ResolveTies<I, fn(I::Item) -> I::Item>
where
    J: IntoIterator<IntoIter = I>,
    I: Iterator,
    I::Item: Copy,
{
    #[inline]
    fn from(iter: J) -> Self {
        ResolveTies::new(iter.into_iter(), |x| x)
    }
}

impl<I, F> Iterator for ResolveTies<I, F>
where
    I: Iterator + Clone,
    I::Item: Copy + PartialEq + Occurrence,
    F: Fn(I::Item) -> I::Item,
{
    type Item = (f64, I::Item);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|item| {
            let normalized_item = (self.normalize)(item);
            if self.current_normalized_item != Some(normalized_item) {
                self.current_normalized_item = Some(normalized_item);
                let count: usize = 1 + self
                    .iter
                    .clone()
                    .map(&self.normalize)
                    .take_while(|x| *x == normalized_item)
                    .map(|x| x.to_occurrence())
                    .sum::<usize>();
                self.resolved = (1.0 + count as f64) / 2.0 + self.index as f64;
                self.index += count;
                self.tie_correction += count.pow(3) - count;
            }
            (self.resolved, item)
        })
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let (low, _high) = self.iter.size_hint();
        (low, None)
    }
}
