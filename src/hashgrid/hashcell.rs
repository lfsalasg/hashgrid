use crate::hashgrid::PeriodicImage;

/// A representation of a cell contained inside a N-dimensional hash grid. It stores elements of type `E` and 
/// can interact directly (with some limitations) or through the `HashGrid` struct.
#[derive(Clone)]
pub struct HashCell<const N:usize, E: Clone> {
    pub dwellers: Vec<E>,
    pub neighbors: Vec<usize>,
    pub periodicity: [PeriodicImage; N] 
}

impl<const N:usize, E: Clone> HashCell<N, E> {
    pub fn new (periodicity:[PeriodicImage; N]) -> Self{
        Self {
            dwellers: Vec::new(),
            neighbors: Vec::new(),
            periodicity
        }
    }

    pub fn get_dwellers(&self) -> &[E] {
        self.dwellers.as_slice()
    }

    pub fn get_mut_dwellers(&mut self) -> &mut [E] {
        self.dwellers.as_mut_slice()
    }

    pub fn set_dwellers(&mut self, dwellers:Vec<E>) {
        self.dwellers = dwellers
    }

    pub fn add_dweller(&mut self, dweller:E) {
        self.dwellers.push(dweller)
    }

    pub fn drop_dweller(&mut self, indx:usize) {
        self.dwellers.remove(indx);
    }

    pub fn population(&self) -> usize {
        self.dwellers.len()
    }
    
}