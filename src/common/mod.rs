//! Exposes common structs, traits and types, useful to work with the 
//! hashgrid crate.
mod point;
mod unittest;

pub use point::{Point, Point2D, Point3D};

#[cfg(feature = "double-precision")]
pub type Float = f64;

#[cfg(not(feature = "double-precision"))]
pub type Float = f32;

/// The Cardinality trait obligates all elements inside a hashgrid to  have
/// a cardinal location in the N-space, where N is the dimension of the grid.
/// It also has some useful space and vector related functions (not the fastest though)
pub trait Cardinality<const N: usize> {
    fn coord(&self) -> Point<N>;

    fn set_coord(&mut self, coord:Point<N>);
}

/// A trait for the elements capable of indexing the `Hashgrid`. The `flatten` method should take a high dimensional
/// index of a cell (an slice for example) and return an index for the position of that `HashCell` given the size of 
/// the grid. The `deflate` method tries to convert an index of a cell to a high dimensional index using the size of
/// the grid.
/// 
/// For readiness, if the dimensionality fails, it will panic instead of returning a `Result` type. 
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

impl<const M: usize> Idx for [usize; M] {
    fn flatten<const N: usize> (&self, grid:[usize; N]) -> usize {
        let mut index = 0;
        
        for i in 0..self.len() {
            if self[i] >= grid[i] {
                panic!("Index is {} but size in the {}-dimension is {}", self[i], i + 1, grid[i])
            }
            index += self[i] * grid[i+1..]
                .iter()
                .fold(1, |acc, x| acc * x);
        }

        index   
    }

    fn deflate<const N: usize>(&self, grid:[usize;N]) -> [usize; N] {
        if self.len() != grid.len() {
            panic!("Index length should be {} but found {}", grid.len(), self.len())
        }
        let mut out = [0; N];
        for i in 0..self.len() {
            out[i] = self[i]
        }
        out
    }
}