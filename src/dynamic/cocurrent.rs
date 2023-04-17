use std::sync::Arc;

use crossbeam::scope;

use crate::common::Cardinality;
use crate::dynamic::{IsoHashGrid, IsoHashGridState};
use crate::hashgrid::{HashGrid, HashCell, ReadGrid, WriteGrid};

/// The MultiThreaded trait provides a set of methods to read and mutate
/// ÌsoHashGrid` instances cocurrently, providing a [serializable isolation
/// level](https://en.wikipedia.org/wiki/Isolation_(database_systems))
pub trait MultiThreaded<const N:usize, E: Clone + Cardinality<N>> {
    fn split_task<T, F>(&mut self, n_threads:usize, f: F)
        where
            F: Fn(&mut HashCell<N, E>, &Arc<HashGrid<N, E>>) -> T + Sync;
}


impl<const N:usize, E:Clone + Cardinality<N> + Sync + Send> MultiThreaded <N, E> for IsoHashGrid<N, E> {
    
    /// Applies the function `f` over each cell in the grid by parallelizing the
    /// task in `n_threads`. Notice that `split_task` might mutate the *future* state
    /// of the IsoHashGrid but it requires to `commit` the changes.
    fn split_task<T, F>(&mut self, n_threads:usize, f: F)
        where
            F: Fn(&mut HashCell<N, E>, &Arc<HashGrid<N, E>>) -> T + Sync
     {
        
        let chunk_size = (self.get_cells_index().len() +  n_threads - 1) / n_threads;
        let shared_present = Arc::new(self.present.clone());
        scope(|s| {
            for chunk in self.get_mut_cells().chunks_mut(chunk_size) {
                let g = &f;
                let local_present = shared_present.clone();
                s.spawn(move |_| {
                    for i in chunk {
                        g(i, &local_present);
                    }
                });
            }
        }).unwrap();
        
    }
}