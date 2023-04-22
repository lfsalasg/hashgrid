//! Utilities for working with the HashGrid family.
mod unittest;

use crate::{HashGrid, ReadGrid, WriteGrid};
use crate::common::{Cardinality, Idx};

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