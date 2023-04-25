//! Contains the basic structs `HashCell` and HashGrid that provide
//! a spatially defined structure to store objects that implement the 
//! Cardinality trait.
//!

mod hashcell;
mod unittest;

use std::ops::{Index, IndexMut};

use serde::{Serialize, Deserialize};
use serde::ser::{Serializer, SerializeStruct};
use serde::de::{Deserializer};

pub use crate::core::hashcell::HashCell;

use crate::common::{Cardinality, Float, Idx, Point};


/// Custom catcheable runtime errors when working with hash grids

#[derive(Debug)]
pub enum HashGridError {
    MismatchedGridSize(String, usize, usize), // Error message, expected cells found cells
    OutOfBounds(String, usize, Float), // Error message, dimension, grid size on this dimension
    WrongDimensionality(String)
}

/// Standard definitions on how to apply the periodic image policies. `LEFT`
/// creates images for the n - 1 face of the cell and `RIGHT` to the n + 1
/// face of the cell.  
#[derive(Clone, Copy, Debug, PartialEq,Serialize, Deserialize)]
pub enum PeriodicImage {
    NONE,
    BOTH,
    LEFT,
    RIGHT
}

/// Common methods for the `HashGrid` and other `HashGrid` based structs to
/// read their content at the grid or cell level.

pub trait ReadGrid <const N: usize, E: Clone + Cardinality<N>> {
    // CELL LEVEL METHODS

    /// Returns a reference of a slice from the cells composing the grid
    fn get_cells(&self) -> &[HashCell<N, E>];

    /// Returns a vector with the indexes of cells
    fn get_cells_index(&self) -> Vec<usize>;

    /// Returns a vector of the coordinates of all cells in the grid
    fn get_cells_coords(&self) -> Vec<[usize; N]>;

    /// Returns an instance of `Point<N>` that represents the anchor of a cell. The anchor is the 
    /// vertex formed by the smallest point bounded by the cell on each dimension.
    fn cell_anchor<I: Idx>(&self, coord:I) -> Point<N>;

    /// Returns an instance of `Point<N>` that represents the center of a cell.
    fn cell_center<I: Idx>(&self, coord:I) -> Point<N>;

    /// Returns the coordinates of the cell that geometrically contains the coordinates `coord`. 
    /// Notice that an element can be *registered* to one cell but can be geometrically located
    /// inside another one.  
    fn bounding_cell_coord(&self, coord:Point<N>) -> Result<[usize; N], HashGridError>;

    // DWELLERS LEVEL METHODS
    
    /// Returns a reference to a slice of all elements inside a cell with coordenates `coord`
    fn get_dwellers<I: Idx>(&self, coord:I) -> &[E];

    /// Returns a vector of tuples. The first element of the tuple is a reference to a neighboring
    /// cell and the second element is the periodic image where the neighbor is located. If the
    /// neighboring cell is located in the real grid, the second element is `[0; N]`
    fn get_neighbors<I: Idx>(&self, coord:I) -> Vec<(&HashCell<N, E>, [isize; N])>;

    /// Returns a vector of tuples. The first element of the tuple is the coordinates of a neighboring
    /// cell and the second element is the periodic image where the neighbor is located. If the
    /// neighboring cell is located in the real grid, the second element is `[0; N]`
    fn get_neighbors_coords<I: Idx>(&self, coord:I) -> Vec<([usize; N], [isize; N])>;

    /// Returns a vector of references to elements registered in the grid. This vector is a 
    /// collection of the elements neighboring the cell with coordinates `coord`.
    fn get_neighbors_dwellers<I: Idx>(&self , coord:I) -> Vec<&E>;

    /// Returns a vector of references to all the elements registered in the grid
    fn get_all_dwellers(&self) -> Vec<&E>;

    // GRID LEVEL METHODS

    /// Shorthand for counting the total number of elements registered in the grid
    fn population(&self) -> usize;

    /// Shorthand for counting the total cells inside the grid
    fn size(&self) -> usize;

    /// Returns a slice with the shape *i.e.* number of cells on each dimension, of the grid.
    fn shape(&self) -> [usize; N];

}

/// Common methods for the `HashGrid` and other `HashGrid` based structs to
/// mutate their content at the grid or cell level.
pub trait WriteGrid <const N: usize, E: Clone + Cardinality<N>> {
    /// Returns a mutable reference of a slice with the cells inside the grid
    fn get_mut_cells(&mut self) -> &mut [HashCell<N, E>];

    /// Returns a mutable reference of the elements registered inside a cell with coordinates
    /// `coord`
    fn get_mut_dwellers<I: Idx>(&mut self, coord:I) -> &mut [E];

    /// Register the elements of a cell with coordinates `coord` using the `dwellers` vector.
    /// **Notice** this method will overwrite any other registered dweller.
    fn set_dwellers<I: Idx>(&mut self, coord:I, dwellers:Vec<E>);   

    /// Register a new `dweller` to the cell with coordinates `coord`
    fn add_dweller<I: Idx>(&mut self, coord:I, dweller:E);

    /// Drops the element with index `indx` registered in the cell with coordinates `coord` and
    /// returns the element. Internally, it performs a `Vec<E>::remove()` operation so it could be
    ///  slow with *O(n)* as the worst case scenario.
    fn drop_dweller<I: Idx>(&mut self, indx:usize, coord:I) -> E;

    /// Drops a dweller from cell `from` and register it into `to`
    fn move_dweller<I: Idx>(&mut self, indx:usize, from:I, to:I);
}

/// An N-dimensional, agnostic grid that provides an interface to interact with its cells and the registered
/// elements. The `HashGrid` struct is defined over `N` dimensions of size `dims` and contain cells of uniform
/// where the elements of type `E` are registered.
#[derive(Clone, Debug)]
pub struct HashGrid<const N:usize, E: Clone + Cardinality<N>> 
{
    grid: [usize; N],
    cells:Vec<HashCell<N, E>>,
    periodicity: [PeriodicImage; N],
    pub dims: Point<N>
}

impl<const N: usize, E: Clone + Cardinality<N>> HashGrid<N, E> {

    /// Creates a uniform grid in the N-dimensional space with the same boundaries and 
    /// periodic conditions
    pub fn generate_uniform_grid(
        grid:[usize; N], 
        periodicity:[PeriodicImage; N], 
        dims:Point<N>) -> Self{
        
        let e_size = grid.iter().fold(1, |acc, x| acc * x);
        let mut hashgrid = HashGrid {
            grid,
            cells: vec![HashCell::<N, E>::new(); e_size],
            periodicity,
            dims
        };
        for i in 0..hashgrid.cells.len() {
            hashgrid.cells[i].neighbors = hashgrid.list_combinations(hashgrid.ndim_from_1dim(i), periodicity)
                .iter()
                .map(|x| (hashgrid.ndim_to_1dim(x.0), x.1))
                .collect()
        }
        hashgrid
    }

    /// Creates a grid in the N-dimensional space starting from a collection of `HashCell`. The 
    /// number of cells should be equal to the expected size of the grid. It returns the `MismatchedSize` error
    /// if the cells passed do not correspond to the dimensions of the Hashgrid
    pub fn generate_from_cells(
        grid:[usize; N], 
        cells:Vec<HashCell<N, E>>,
        periodicity:[PeriodicImage; N], 
        dims:Point<N>) -> Result<HashGrid<N, E>, HashGridError>{
        
        let expected_length = grid.iter().fold(1, |acc, x| acc * x);
        if cells.len() !=  expected_length {
            return Err(HashGridError::MismatchedGridSize(
                format!("Expected number of cells was {} but  {} were found", expected_length, cells.len()),
                expected_length,
                cells.len()
            ));
        }

        let mut hashgrid = HashGrid {
            grid,
            cells,
            periodicity,
            dims
        };

        for i in 0..hashgrid.cells.len() {
            hashgrid.cells[i].neighbors = hashgrid.list_combinations(hashgrid.ndim_from_1dim(i), hashgrid.periodicity)
                .iter()
                .map(|x| (hashgrid.ndim_to_1dim(x.0), x.1))
                .collect()
        }

        Ok(hashgrid)
    }
    
    fn ndim_to_1dim(&self, coord:[usize; N]) -> usize {
        let mut index = 0;
        
        for i in 0..self.grid.len() {
            if coord[i] >= self.grid[i] {
                panic!("Index is {} but size in the {}-dimension is {}", coord[i], i + 1, self.grid[i])
            }
            index += coord[i] * self.grid[i+1..]
                .iter().fold(1, |acc, x| acc * x);
        }

        index
    }

    fn ndim_from_1dim(&self, indx:usize) -> [usize; N] {
        let mut indices:[usize; N] = [0; N];
        let mut remainder = indx;

        for i in (0..self.grid.len()).rev() {
            let idx = remainder % self.grid[i];
            indices[i] = idx;
            remainder /= self.grid[i];
        }

        indices
    }

    /// Lists all n-dimensional indexes for the periodic images of a central cell located at `coord`.
    /// The definition of what is considered a *neighbor* depends on the value of `periodic_image`
    /// where the position *i* of the array corresponds to the ith dimension and the value indicates
    /// how periodicity should be applied on this dimension (face)   
    fn list_combinations<I: Idx>(&self, coord:I, 
        periodic_images:[PeriodicImage;N]) -> Vec<([usize; N], [isize; N])> {

        let cell = coord.deflate(self.grid);
        let mut all_combs = Vec::new();
        self.list_combinations_helper(0, cell, [0; N],&mut all_combs, &periodic_images);
        all_combs.remove(0);
        all_combs
    }

    fn list_combinations_helper(&self, i: usize, 
        comb: [usize; N],
        grid: [isize; N], 
        all_combs: &mut Vec<([usize; N], [isize; N])>, 
        periodic_images:&[PeriodicImage;N]) {
        /*
        let translations:Vec<isize> = match periodic_images[i] {
            PeriodicImage::NONE => {vec![0]},
            PeriodicImage::LEFT => {vec![0, -1]},
            PeriodicImage::RIGHT => {vec![0, 1]},
            PeriodicImage::BOTH => {vec![0, -1, 1]}
        };
        */
        let translations = [0, -1, 1];
        if i == N - 1 {
            for k in translations { // Modifies the last dimension
                let mut cell = comb.clone();
                let mut grid = grid.clone();

                let dim = cell[i]  as isize + k;
                
                if dim >= 0 && dim < self.grid[i] as isize{
                    cell[i] = dim as usize;
                    grid[i] = 0;
                    all_combs.push((cell, grid))
                }
                else if dim < 0 && (periodic_images[i] == PeriodicImage::LEFT || periodic_images[i] == PeriodicImage::BOTH) {
                    cell[i] = self.grid[i] - 1;
                    grid[i] = -1;
                    all_combs.push((cell, grid)); 
                }else if dim >= self.grid[i] as isize &&  (periodic_images[i] == PeriodicImage::RIGHT || periodic_images[i] == PeriodicImage::BOTH) {
                    cell[i] = 0;
                    grid[i] = 1;
                    all_combs.push((cell, grid)); 
                }                
            }
            return
        }
        else {
            for k in translations { // Modifies the kth dimension
                let mut cell = comb.clone();
                let mut grid = grid.clone();
                let dim = cell[i]  as isize + k;
                if dim >= 0 && dim < self.grid[i] as isize{
                    cell[i] = dim as usize;
                    grid[i] = 0;

                    self.list_combinations_helper(i + 1, cell, grid, all_combs, periodic_images)
                }
                else if dim < 0 && (periodic_images[i] == PeriodicImage::LEFT || periodic_images[i] == PeriodicImage::BOTH) {
                    cell[i] = self.grid[i] - 1;
                    grid[i] = -1;

                    self.list_combinations_helper(i + 1, cell, grid, all_combs, periodic_images)
                }else if dim >= self.grid[i] as isize &&  (periodic_images[i] == PeriodicImage::RIGHT || periodic_images[i] == PeriodicImage::BOTH) {
                    cell[i] = 0;
                    grid[i] = 1;

                    self.list_combinations_helper(i + 1, cell, grid, all_combs, periodic_images)
                }     
            }
        }
    }
}

impl <const N: usize, E: Clone + Cardinality<N>> ReadGrid<N, E> for HashGrid<N, E> {
    fn get_cells(&self) -> &[HashCell<N, E>] {
        self.cells.as_slice()
    }

    fn get_cells_index(&self) -> Vec<usize> {
        let v:Vec<usize> = (0..self.cells.len()).collect();
        v
    }

    fn get_cells_coords(&self) -> Vec<[usize; N]> {
        let coords:Vec<[usize; N]> = (0..self.cells.len())
                                        .into_iter()
                                        .map(|x| self.ndim_from_1dim(x))
                                        .collect();
        coords
    }
    
    // Returns a referece of the elements registered under a cell with coordinates `coord` in the
    /// N-dimensional space    
    fn get_dwellers<I: Idx>(&self, coord:I) -> &[E] {
        let indx = coord.flatten(self.grid);
        //let indx = self.ndim_to_1dim(coord);
        self.cells[indx].get_dwellers()
    }

    fn get_neighbors<I: Idx>(&self, coord:I) -> Vec<(&HashCell<N, E>, [isize; N])> {
        let mut neighbors = Vec::new();
        let indx = coord.flatten(self.grid);
        for (cell_index, grid) in self.cells[indx].neighbors.iter() {
            neighbors.push((&self.cells[*cell_index], *grid))
        } 

        neighbors
    }

    fn get_neighbors_coords<I: Idx>(&self, coord:I) -> Vec<([usize; N], [isize; N])> {
        let mut neighbors = Vec::new();
        let indx = coord.flatten(self.grid);
        for (cell_index, grid) in self.cells[indx].neighbors.iter() {
            neighbors.push((self.ndim_from_1dim(*cell_index), *grid))
        } 

        neighbors
    }

    fn get_neighbors_dwellers<I: Idx>(&self , coord:I) -> Vec<&E> {
        let mut neighbors = Vec::new();
        let indx = coord.flatten(self.grid);
        for (cell_index, _) in self.cells[indx].neighbors.iter() {
            neighbors.extend(self.cells[*cell_index].get_dwellers())
        } 

        neighbors
    }   

    fn get_all_dwellers(&self) -> Vec<&E> {
        let mut pop:Vec<&E> = Vec::new();

        for cell in self.cells.iter() {
            pop.extend(cell.get_dwellers())
        }

        pop
    }

    fn population(&self) -> usize {
        self.cells.iter().fold(0, |acc, x| acc + x.population())
    }

    fn size(&self) -> usize {
        self.cells.len()
    }

    fn shape(&self) -> [usize; N] {
        self.grid
    }

    fn cell_anchor<I: Idx>(&self, coord:I) -> Point<N> {
        let cell = coord.deflate(self.grid);
        let cell_size = self.dims / self.grid;

        let anchor = cell_size * cell;
        
        anchor
    }

    fn cell_center<I: Idx>(&self, coord:I) -> Point<N> {
        let cell = coord.deflate(self.grid);
        let cell_size = self.dims / self.grid;

        let anchor = cell_size * cell;
        let center = anchor + cell_size / 2.0;

        center
    }

    fn bounding_cell_coord(&self, coord:Point<N>) -> Result<[usize; N], HashGridError> {
        let mut cell:[usize; N] = [0; N];
        //if coord >= Point::from_scalar(0.0) && coord <= self.dims {
        //    println!("{:?}, {:?}", coord, self.dims);
        //    let rel_position = coord / self.dims * self.grid;
        //    for dim in 0..N {
        //        cell[dim] = rel_position[dim] as usize;
        //    }
//
        //    return Ok(cell)
//
        //}

        for dim in 0..N {
            
            if coord[dim] >= 0.0 && coord[dim] <= self.dims[dim] {
                cell[dim] = (coord[dim] / self.dims[dim] * self.grid[dim] as Float) as usize
            }else if coord[dim] < 0.0 {
                if self.periodicity[dim] == PeriodicImage::LEFT || self.periodicity[dim] == PeriodicImage::BOTH {
                    let rel_cell = coord[dim] / self.dims[dim];  
                    cell[dim] = ((self.dims[dim] - (rel_cell - rel_cell.ceil())) *  self.grid[dim] as Float) as usize;
                    cell[dim] = ((self.dims[dim] + coord[dim] % self.dims[dim]) / self.dims[dim] * self.grid[dim] as Float) as usize;
                }else {
                    return Err(HashGridError::OutOfBounds(
                        format!("In dimension {} object is {} from origin but grid cannot have negative values and it does not implement Periodicity", dim, coord[dim]), 
                        dim, coord[dim]))
                }
            } else {
                if self.periodicity[dim] == PeriodicImage::RIGHT || self.periodicity[dim] == PeriodicImage::BOTH {
                    let rel_cell = coord[dim] / self.dims[dim];                    
                    cell[dim] = ((rel_cell - rel_cell.floor()) * self.grid[dim] as Float) as usize;

                }else {
                    return Err(HashGridError::OutOfBounds(
                        format!("In dimension {} object is {} from origin but grid has size {} and it does not implement Periodicity", dim, coord[dim], self.dims[dim]), 
                        dim, coord[dim]))
                }
            }
        }

        Ok(cell)
    }
}

impl<const N: usize, E: Clone + Cardinality<N>> WriteGrid<N, E> for HashGrid<N, E> {
    /// Returns a mutable reference of a slice from the cells composing the GRID
    fn get_mut_cells(&mut self) -> &mut [HashCell<N, E>] {
        self.cells.as_mut_slice()
    }

    fn get_mut_dwellers<I: Idx>(&mut self, coord:I) -> &mut [E] {
        let indx = coord.flatten(self.grid);
        self.cells[indx].get_mut_dwellers()
    }

    fn set_dwellers<I: Idx>(&mut self, coord:I, dwellers:Vec<E>) {
        let indx = coord.flatten(self.grid);
        self.cells[indx].set_dwellers(dwellers)
    }

    fn add_dweller<I: Idx>(&mut self, coord:I, dweller:E) {
        let indx = coord.flatten(self.grid);
        self.cells[indx].add_dweller(dweller)    
    }

    fn drop_dweller<I: Idx>(&mut self, indx:usize, coord:I) -> E {
        let cell_index = coord.flatten(self.grid);
        self.cells[cell_index].drop_dweller(indx)
    }

    fn move_dweller<I: Idx>(&mut self, indx:usize, from:I, to:I) {
        let from_indx = from.flatten(self.grid);
        let to_indx = to.flatten(self.grid);
        let dw = self.cells[from_indx].drop_dweller(indx);
        self.cells[to_indx].add_dweller(dw);
    }

}

impl<const N: usize, E: Clone + Cardinality<N>, I: Idx> Index<I> for HashGrid<N, E>{
    type Output = HashCell<N, E>;
    fn index(&self, index: I) -> &Self::Output {

        &self.cells[index.flatten(self.grid)]
    }
}

impl< const N: usize, E: Clone + Cardinality<N>, I: Idx> IndexMut<I> for HashGrid<N, E>{
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        &mut self.cells[index.flatten(self.grid)]
    }
}

impl<const N: usize, E: Clone + Serialize + Cardinality<N>> Serialize for HashGrid<N, E> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("HashGrid", 3)?;
        state.serialize_field("grid", &self.grid.to_vec())?;
        state.serialize_field("cells", &self.cells)?;
        state.serialize_field("periodicity", &self.periodicity.to_vec())?;
        state.serialize_field("dims", &self.dims.to_vec())?;
        state.end()
    }
}

impl<'de, const N: usize, E: Clone + Deserialize<'de> + Cardinality<N>> Deserialize<'de> for HashGrid<N, E> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(serde::Deserialize)]
        struct HashGridHelper<const N:usize, E: Clone + Cardinality<N>> {
            grid: Vec<usize>,
            cells: Vec<HashCell<N, E>>,
            periodicity: Vec<PeriodicImage>,
            dims: Vec<Float>
        }

        let helper = HashGridHelper::deserialize(deserializer)?;
        Ok(Self {
            grid: helper.grid.try_into().unwrap(),
            cells: helper.cells,
            periodicity: helper.periodicity.try_into().unwrap(),
            dims: Point::from_vec(helper.dims).unwrap()
        })
    }
}