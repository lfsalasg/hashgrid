use std::sync::Arc;
use crate::common::Cardinality;
use crate::hashgrid::HashGrid;

type AtomicHashGrid<const N:usize, E> = Arc<HashGrid<N, E>>;

struct AIsoHashGrid<const N:usize, E: Clone + Cardinality<N>> {
    images:[AtomicHashGrid<N, E>; 2],
    p: usize,
    f: usize
}