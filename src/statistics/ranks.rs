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

        let mut resolved_ties = ResolveTies::from(observations.iter().map(|(_, value)| *value));
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
    I::Item: Copy + PartialEq,
    F: Fn(I::Item) -> I::Item,
{
    type Item = (f64, I::Item);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|item| {
            let normalized_item = (self.normalize)(item);
            if self.current_normalized_item != Some(normalized_item) {
                self.current_normalized_item = Some(normalized_item);
                let count = 1 + self
                    .iter
                    .clone()
                    .map(&self.normalize)
                    .take_while(|x| *x == normalized_item)
                    .count();
                self.resolved = (1.0 + count as f64) / 2.0 + self.index as f64;
                self.tie_correction += count.pow(3) - count;
            }
            self.index += 1;
            (self.resolved, item)
        })
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let (low, _high) = self.iter.size_hint();
        (low, None)
    }
}
