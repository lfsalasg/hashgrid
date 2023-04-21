mod cocurrent;
mod unittest;

use crate::core::{HashGrid, ReadGrid, WriteGrid};
use crate::common::{Cardinality, Point};
pub use crate::dynamic::cocurrent::MultiThreaded;

#[derive(PartialEq)]
enum IsoHashGridState {
    Committed,
    Initialized
}

/// The IsoHashGrid is a struct build upon the basic `HashGrid` struct and stands
/// for *Isolated HashGrid*. The grid provides two states of a `HashGrid`: A *present*
/// or *current* view of the grid and a future, mutable view of it. It reproduces a 
/// snapshot isolation strategy, where all read transactions are performed over the
/// present view of the grid while all write transactions occur on the future view. 
/// 
pub struct IsoHashGrid<const N:usize, E: Clone + Cardinality<N>> {
    present: HashGrid<N, E>,
    future: HashGrid<N, E>,
    state: IsoHashGridState
}

impl<const N:usize, E:Clone + Cardinality<N>> IsoHashGrid<N, E> {
    /// Creates an IsoHashGrid from a `HashGrid` instance
    pub fn from(grid:HashGrid<N, E>) -> IsoHashGrid<N, E> {
        Self {
            present: grid.clone(),
            future: grid,
            state: IsoHashGridState::Initialized
        }
    }

    /// Updates the *present* view with the information of the *future* view.
    /// Currently, it implements a naive `clone` strategy to clone the data
    pub fn commit(&mut self){
        self.present = self.future.clone();

        self.state = IsoHashGridState::Committed;
    }

    /// Removes any changes in the *future* view, making it identical to the 
    /// *present* view. Currently, it implements a naive `clone` strategy to
    /// clone the data.
    pub fn rollback(&mut self) {
        self.future = self.present.clone();
    }
}

impl <const N:usize, E: Clone + Cardinality<N>> ReadGrid<N,E> for IsoHashGrid<N, E> {
    fn cell_anchor<I: crate::common::Idx>(&self, coord:I) -> Point<N> {
        self.present.cell_anchor(coord)
    }
    
    fn cell_center<I: crate::common::Idx>(&self, coord:I) -> Point<N> {
        self.present.cell_center(coord)
    }

    fn get_cells(&self) -> &[crate::HashCell<N, E>] {
        self.present.get_cells()
    }

    fn get_cells_coords(&self) -> Vec<[usize; N]> {
        self.present.get_cells_coords()
    }

    fn get_cells_index(&self) -> Vec<usize> {
        self.present.get_cells_index()
    }

    fn bounding_cell_coord(&self, coord:Point<N>) -> Result<[usize; N], crate::HashGridError> {
        self.present.bounding_cell_coord(coord)
    }

    fn get_dwellers<I: crate::common::Idx>(&self, coord:I) -> &[E] {
        self.present.get_dwellers(coord)
    }

    fn get_neighbors<I: crate::common::Idx>(&self, coord:I) -> Vec<(&crate::HashCell<N, E>, [isize; N])> {
        self.present.get_neighbors(coord)
    }

    fn get_neighbors_coords<I: crate::common::Idx>(&self, coord:I) -> Vec<([usize; N], [isize; N])> {
        self.present.get_neighbors_coords(coord)
    }

    fn get_neighbors_dwellers<I: crate::common::Idx>(&self , coord:I) -> Vec<&E> {
        self.present.get_neighbors_dwellers(coord)
    }

    fn get_all_dwellers(&self) -> Vec<&E> {
        self.present.get_all_dwellers()
    }

    fn population(&self) -> usize {
        self.present.population()
    }

    fn shape(&self) -> [usize; N] {
        self.present.shape()
    }

    fn size(&self) -> usize {
        self.present.size()
    }
}

impl <const N:usize, E: Clone + Cardinality<N>> WriteGrid<N,E> for IsoHashGrid<N, E> {
    
    fn get_mut_cells(&mut self) -> &mut [crate::HashCell<N, E>] {
        self.future.get_mut_cells()
    }

    fn get_mut_dwellers<I: crate::common::Idx>(&mut self, coord:I) -> &mut [E] {
        self.future.get_mut_dwellers(coord)
    }

    fn set_dwellers<I: crate::common::Idx>(&mut self, coord:I, dwellers:Vec<E>) {
        self.future.set_dwellers(coord, dwellers)
    }
    
    fn add_dweller<I: crate::common::Idx>(&mut self, coord:I, dweller:E) {
        self.future.add_dweller(coord, dweller)
    }

    fn drop_dweller<I: crate::common::Idx>(&mut self, indx:usize, coord:I) {
        self.future.drop_dweller(indx, coord)
    }

    fn move_dweller<I: crate::common::Idx>(&mut self, indx:usize, from:I, to:I) {
        self.future.move_dweller(indx, from, to)
    }

    fn update_neighbors<I: crate::common::Idx>(&mut self, coord:I, periodic_images:[crate::PeriodicImage;N]) {
        self.future.update_neighbors(coord, periodic_images)
    }
}