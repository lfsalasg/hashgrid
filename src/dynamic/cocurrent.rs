use std::sync::Arc;

use crossbeam::scope;

use crate::common::Cardinality;
use crate::dynamic::IsoHashgrid;
use crate::hashgrid::{HashCell, ReadGrid, WriteGrid};

pub trait MultiThreaded<const N:usize, E: Clone + Cardinality<N>> {
    fn split_task<T, F: Fn(&mut HashCell<N, E>) -> T + Sync>(&mut self, n_threads:usize, f: F);

    fn split_and_move<F: Fn(&HashCell<N, E>) -> Vec<(E, usize, usize)>>(&mut self, n_threads:usize, f: F);
}

impl<const N:usize, E:Clone + Cardinality<N> + Sync + Send> MultiThreaded <N, E> for IsoHashgrid<N, E> {
    fn split_task<T, F: Fn(&mut HashCell<N, E>) -> T + Sync>(&mut self, n_threads:usize, f: F) {
        let chunk_size = (self.get_cells_index().len() +  n_threads - 1) / n_threads;
    
        scope(|s| {
            for chunk in self.get_mut_cells().chunks_mut(chunk_size) {
                s.spawn(|_| {
                    for i in chunk {
                        f(i);
                    }
                });
            }
            /*
            for _ in 0..n_threads {
                let start = chunk_start;
                let end = chunk_end;
                let g = &f;
                let mut image = &self.images[self.f];
                s.spawn(move |_| {
                    for i in start..end {
                        g(&mut image[i]);
                    }
                });
                chunk_start = chunk_end;
                chunk_end = (chunk_end + chunk_size).min(self.get_cells_index().len());
                
            }
            */   
        }).unwrap();
        
    }

    fn split_and_move<F: Fn(&HashCell<N, E>) -> Vec<(E, usize, usize)>>(&mut self, n_threads:usize, f: F) {
        for cell in self.get_cells() {
            let _movements = f(&cell);
        }
    }
}