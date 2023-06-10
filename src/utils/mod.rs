//! Utilities for working with the HashGrid family.
mod unittest;

use rand::{thread_rng, Rng};

use crate::{HashGrid, ReadGrid, WriteGrid};
use crate::common::{Cardinality, Idx, Point};

pub fn relocate_dwellers<const N: usize, E: Cardinality<N> + Clone> (grid:&mut HashGrid<N, E>) {
    let mut changes = Vec::new();
    for (indx, cell) in grid.get_cells().iter().enumerate() {
        let current_index = indx.deflate(grid.shape());
        for (dw_index, dweller) in cell.get_dwellers().iter().enumerate() {
            match grid.bounding_cell_coord(dweller.coord()) {
                Ok(new_index) => {
                    if new_index != current_index {
                        changes.push((dw_index, current_index, new_index))
                    }
                },

                Err(e) => {
                    panic!("{:?}", e)
                }   
            }
        }
    }

    for (dw_index, from, to) in changes {
        grid.move_dweller(dw_index, from, to)
    }
}

pub fn populate_randomly<const N:usize, E:Clone + Cardinality<N>>(grid: &mut HashGrid<N, E>, mut elements:Vec<E>) {
    let mut rng = thread_rng();
    while let Some(mut element) = elements.pop() {
        let movement:Point<N> = grid.dims()
            .iter()
            .map(|x| rng.gen_range(0.0..*x))
            .collect();

        element.set_coord(element.coord() + movement);

        let cell_coord = grid.bounding_cell_coord(element.coord()).expect("Runtime error: Unexpected error found");
        grid[cell_coord].add_dweller(element);
        
    } 
}