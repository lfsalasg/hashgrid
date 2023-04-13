mod cocurrent;
mod unittest;

use crate::hashgrid::{HashGrid, ReadGrid, WriteGrid};
use crate::common::Cardinality;
use crate::dynamic::cocurrent::MultiThreaded;

#[derive(PartialEq)]
enum IsoHashgridState {
    Committed,
    Initialized
}

pub struct IsoHashgrid<const N:usize, E: Clone + Cardinality<N>> {
    images:[HashGrid<N, E>; 2],
    p: usize,
    f: usize,
    state: IsoHashgridState
}

impl<const N:usize, E:Clone + Cardinality<N>> IsoHashgrid<N, E> {
    pub fn from(grid:HashGrid<N, E>) -> IsoHashgrid<N, E> {
        Self {
            images: [grid.clone(), grid],
            p: 0,
            f: 1,
            state: IsoHashgridState::Initialized
        }
    }

    pub fn commit(&mut self){
        if self.p == 0 {
            self.p = 1;
            self.f = 0
        }else {
            self.p = 0;
            self.f = 1
        }

        self.state = IsoHashgridState::Committed;
    }

    pub fn rollback(&mut self) {
        self.images[self.f] = self.images[self.p].clone()
    }
}

impl <const N:usize, E: Clone + Cardinality<N>> ReadGrid<N,E> for IsoHashgrid<N, E> {
    fn cell_center<I: crate::common::Idx>(&self, coord:I) -> [crate::common::Float; N] {
        self.images[self.p].cell_center(coord)
    }

    fn get_cells(&self) -> &[crate::HashCell<N, E>] {
        self.images[self.p].get_cells()
    }

    fn get_cells_coords(&self) -> Vec<[usize; N]> {
        self.images[self.p].get_cells_coords()
    }

    fn get_cells_index(&self) -> Vec<usize> {
        self.images[self.p].get_cells_index()
    }

    fn get_bounding_cell(&self, coord:[crate::common::Float; N]) -> Result<[usize; N], crate::HashGridError> {
        self.images[self.p].get_bounding_cell(coord)
    }

    fn get_dwellers<I: crate::common::Idx>(&self, coord:I) -> &[E] {
        self.images[self.p].get_dwellers(coord)
    }

    fn get_neighbors<I: crate::common::Idx>(&self, coord:I) -> Vec<(&crate::HashCell<N, E>, [isize; N])> {
        self.images[self.p].get_neighbors(coord)
    }

    fn get_neighbors_coords<I: crate::common::Idx>(&self, coord:I) -> Vec<([usize; N], [isize; N])> {
        self.images[self.p].get_neighbors_coords(coord)
    }

    fn get_neighbors_dwellers<I: crate::common::Idx>(&self , coord:I) -> Vec<&E> {
        self.images[self.p].get_neighbors_dwellers(coord)
    }

    fn get_all_dwellers(&self) -> Vec<&E> {
        self.images[self.p].get_all_dwellers()
    }

    fn population(&self) -> usize {
        self.images[self.p].population()
    }

    fn shape(&self) -> [usize; N] {
        self.images[self.p].shape()
    }

    fn size(&self) -> usize {
        self.images[self.p].size()
    }
}

impl <const N:usize, E: Clone + Cardinality<N>> WriteGrid<N,E> for IsoHashgrid<N, E> {
    
    fn get_mut_cells(&mut self) -> &mut [crate::HashCell<N, E>] {
        self.images[self.f].get_mut_cells()
    }

    fn get_mut_dwellers<I: crate::common::Idx>(&mut self, coord:I) -> &mut [E] {
        self.images[self.f].get_mut_dwellers(coord)
    }

    fn set_dwellers<I: crate::common::Idx>(&mut self, coord:I, dwellers:Vec<E>) {
        self.images[self.f].set_dwellers(coord, dwellers)
    }
    
    fn add_dweller<I: crate::common::Idx>(&mut self, coord:I, dweller:E) {
        self.images[self.f].add_dweller(coord, dweller)
    }

    fn drop_dweller<I: crate::common::Idx>(&mut self, indx:usize, coord:I) {
        self.images[self.f].drop_dweller(indx, coord)
    }

    fn move_dweller<I: crate::common::Idx>(&mut self, indx:usize, from:I, to:I) {
        self.images[self.f].move_dweller(indx, from, to)
    }

    fn update_neighbors<I: crate::common::Idx>(&mut self, coord:I, periodic_images:[crate::PeriodicImage;N]) {
        self.images[self.f].update_neighbors(coord, periodic_images)
    }
}