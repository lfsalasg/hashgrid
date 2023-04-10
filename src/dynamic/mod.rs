use std::ops::Deref;

use crate::hashgrid::{HashGrid, PeriodicImage};
use crate::common::{Cardinality, Point};

struct DynamiHashgrid<const N:usize, E: Clone + Cardinality<N>> {
    images:(HashGrid<N, E>, HashGrid<N, E>),
    pointer: usize
}

impl<const N:usize, E:Clone + Cardinality<N>> DynamiHashgrid<N, E> {
    pub fn generate_uniform_grid(grid:[usize; N], periodicity:[PeriodicImage; N], dims:Point<N>) -> Self {
        let grid = HashGrid::generate_uniform_grid(grid, periodicity, dims);
        let ref_grid = &grid;
        
        Self {
            images: (grid.clone(), grid),
            pointer: 0
        }
    }

    pub fn commit(&mut self){
        if self.pointer == 0 {
            self.pointer = 1
        }else {
            self.pointer = 0
        }
    }
}

impl<const N: usize, E: Clone + Cardinality<N>> Deref for DynamiHashgrid<N, E> {
    type Target = HashGrid<N,E>;

    fn deref(&self) -> &Self::Target {
        if self.pointer == 0 {
            return &self.images.0
        }
        else if self.pointer == 1 {
            return &self.images.1;
        }
        else {
            panic!("Pointer is pointing to a non-existent image")
        }
    }
}