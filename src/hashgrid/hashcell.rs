use serde::ser::{Serialize, Serializer, SerializeStruct};
use serde::de::{Deserialize, Deserializer};
use crate::hashgrid::PeriodicImage;

/// A representation of a cell contained inside a N-dimensional hash grid. It stores elements of type `E` and 
/// can interact directly (with some limitations) or through the `HashGrid` struct.
#[derive(Clone, Debug)]
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


impl<const N: usize, E: Clone + Serialize> Serialize for HashCell<N, E> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("HashCell", 3)?;
        state.serialize_field("dwellers", &self.dwellers)?;
        state.serialize_field("neighbors", &self.neighbors)?;
        state.serialize_field("periodicity", &self.periodicity.to_vec())?;
        state.end()
    }
}


impl<'de, const N: usize, E: Clone + Deserialize<'de>> Deserialize<'de> for HashCell<N, E> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(serde::Deserialize)]
        struct HashCellHelper<E: Clone> {
            dwellers: Vec<E>,
            neighbors: Vec<usize>,
            periodicity: Vec<PeriodicImage>
        }

        let helper = HashCellHelper::deserialize(deserializer)?;
        Ok(Self {
            dwellers: helper.dwellers,
            neighbors: helper.neighbors,
            periodicity: helper.periodicity.try_into().unwrap()
        })
    }
}