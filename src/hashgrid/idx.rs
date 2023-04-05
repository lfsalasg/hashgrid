/// A trait for the elements capable of indexing the `Hashgrid`. The `flatten` method should take a high dimensional
/// index of a cell (an slice for example) and return an index for the position of that `HashCell` given the size of 
/// the grid. The `deflate` method tries to convert an index of a cell to a high dimensional index using the size of
/// the grid.
/// 
/// For readiness, if the dimensionality fails, it would panic instead of returning a `Result` type. 
pub trait Idx {
    fn flatten<const N: usize>(&self, grid:[usize; N]) -> usize;
    fn deflate<const N: usize>(&self, grid:[usize;N]) -> [usize; N];
}

impl Idx for usize {
    fn flatten<const N: usize>(&self, _grid:[usize; N]) -> usize {
        *self
    }

    fn deflate<const N: usize>(&self, grid:[usize;N]) -> [usize; N] {
        if *self > grid.iter().fold(1, |acc, x| acc * x) {
            panic!("Deflating to something larger than the grid's size")
        }
        let mut indices:[usize; N] = [0; N];
        let mut remainder = *self;

        for i in (0..grid.len()).rev() {
            let idx = remainder % grid[i];
            indices[i] = idx;
            remainder /= grid[i];
        }

        indices 
    }
}