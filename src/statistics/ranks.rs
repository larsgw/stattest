use std::iter::FusedIterator;

/// For the ranking of groups of variables.
pub trait Ranks<T> {
    /// Returns a vector of ranks.
    fn ranks(self) -> (Vec<T>, usize);
}

impl<T> Ranks<f64> for T
where
    T: IntoIterator,
    T::Item: PartialOrd,
{
    fn ranks(self) -> (Vec<f64>, usize) {
        let mut observations: Vec<(usize, T::Item)> = self.into_iter().enumerate().collect();
        observations.sort_unstable_by(|(_, a), (_, b)| {
            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
        });

        let (sorted_indices, sorted_values): (Vec<_>, Vec<_>) = observations.into_iter().unzip();
        let groups = DedupWithCount::new(sorted_values.iter());

        let mut resolved_ties = ResolveTies::new(groups);
        let mut ranks = vec![0.0; sorted_values.len()];

        for (rank, old_index) in (&mut resolved_ties).zip(sorted_indices) {
            ranks[old_index] = rank;
        }

        (ranks, resolved_ties.tie_correction())
    }
}

struct ResolveTies<I, T>
where
    I: Iterator<Item = (usize, T)>,
{
    iter: I,
    index: usize,
    left: usize,
    resolved: f64,
    tie_correction: usize,
}

impl<I, T> ResolveTies<I, T>
where
    I: Iterator<Item = (usize, T)>,
{
    fn new(iter: I) -> Self {
        ResolveTies {
            iter,
            index: 0,
            left: 0,
            resolved: 0.0,
            tie_correction: 0,
        }
    }

    fn tie_correction(&self) -> usize {
        self.tie_correction
    }
}

impl<I, T> Iterator for ResolveTies<I, T>
where
    I: Iterator<Item = (usize, T)>,
{
    type Item = f64;

    fn next(&mut self) -> Option<f64> {
        if self.left > 0 {
            self.left -= 1;
            self.index += 1;
            Some(self.resolved)
        } else {
            self.iter.next().map(|(count, _)| {
                self.resolved = ((1.0 + count as f64) / 2.0) + self.index as f64;
                self.left = count - 1;
                self.index += 1;
                self.tie_correction += count.pow(3) - count;
                self.resolved
            })
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (low, _high) = self.iter.size_hint();
        (low, None)
    }
}

impl<I, T> FusedIterator for ResolveTies<I, T> where I: Iterator<Item = (usize, T)> {}

struct DedupWithCount<I>
where
    I: Iterator,
    I::Item: PartialEq,
{
    iter: I,
    current_value: Option<I::Item>,
}

impl<I> DedupWithCount<I>
where
    I: Iterator,
    I::Item: PartialEq,
{
    fn new(mut iter: I) -> Self {
        DedupWithCount {
            current_value: iter.next(),
            iter,
        }
    }
}

impl<I> Iterator for DedupWithCount<I>
where
    I: Iterator,
    I::Item: PartialEq + Copy,
{
    type Item = (usize, I::Item);

    fn next(&mut self) -> Option<(usize, I::Item)> {
        self.current_value.map(|current_value| {
            let mut count = 1;
            self.current_value = self.iter.next();
            while self
                .current_value
                .map(|next_value| next_value == current_value)
                .unwrap_or(false)
            {
                count += 1;
                self.current_value = self.iter.next();
            }
            (count, current_value)
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (low, high) = self.iter.size_hint();
        (if low == 0 { 0 } else { 1 }, high)
    }
}

impl<I: Iterator> FusedIterator for DedupWithCount<I> where I::Item: PartialEq + Copy {}
