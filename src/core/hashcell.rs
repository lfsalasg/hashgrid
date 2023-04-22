use serde::ser::{Serialize, Serializer, SerializeStruct};
use serde::de::{Deserialize, Deserializer};

use crate::common::Cardinality;

/// A representation of a cell contained inside a N-dimensional hash grid. It stores elements of type `E` and 
/// can interact directly (with some limitations) or through the `HashGrid` struct.
#[derive(Clone, Debug)]
pub struct HashCell<const N:usize, E: Clone + Cardinality<N>> {
    pub dwellers: Vec<E>,
    pub neighbors: Vec<(usize, [isize; N])>
}

impl<const N:usize, E: Clone + Cardinality<N>> HashCell<N, E> {
    pub fn new () -> Self{
        Self {
            dwellers: Vec::new(),
            neighbors: Vec::new()
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

    pub fn drop_dweller(&mut self, indx:usize) -> E{
        self.dwellers.remove(indx)
    }

    pub fn population(&self) -> usize {
        self.dwellers.len()
    }
    
}


impl<const N: usize, E: Clone + Serialize + Cardinality<N>> Serialize for HashCell<N, E> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let neighbors:Vec<(usize, Vec<isize>)> = self.neighbors.iter()
                        .map(|x| (x.0, x.1.to_vec()))
                        .collect();
        let mut state = serializer.serialize_struct("HashCell", 3)?;
        state.serialize_field("dwellers", &self.dwellers)?;
        state.serialize_field("neighbors", &neighbors)?;
        state.end()
    }
}


impl<'de, const N: usize, E: Clone + Deserialize<'de> + Cardinality<N>> Deserialize<'de> for HashCell<N, E> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(serde::Deserialize)]
        struct HashCellHelper<E: Clone> {
            dwellers: Vec<E>,
            neighbors: Vec<(usize, Vec<isize>)>,
        }

        let helper = HashCellHelper::deserialize(deserializer)?;
        let neighbors:Vec<(usize, [isize; N])> = helper.neighbors.iter()
                        .map(|x| (x.0, x.1.clone().try_into().unwrap()))
                        .collect();
        Ok(Self {
            dwellers: helper.dwellers,
            neighbors,
        })
    }
}