mod unittest;
mod hashcell;

use std::ops::{Index, IndexMut};

pub use crate::hashgrid::hashcell::HashCell;

#[cfg(feature = "double-precision")]
type Float = f64;

#[cfg(not(feature = "double-precision"))]
type Float = f32;

#[derive(Debug)]
pub enum HashGridError {
    MismatchedSize(String),
    OutOfBounds(String)
}

#[derive(Clone, Copy)]
pub enum PeriodicImage {
    NONE,
    BOTH,
    LEFT,
    RIGHT
}

/// An N-dimensional, agnostic grid that provides an interface to interact with its cells and the registered
/// elements. The `HashGrid` struct is defined over `N` dimensions of size `[T; N]` and contain cells of uniform
/// size `dims`where the elements of type `E` are registered.
pub struct HashGrid<const N:usize, E: Clone> 
{
    grid: [usize; N],
    cells:Vec<HashCell<N, E>>,
    pub dims: [Float; N]
}

impl<const N: usize, E: Clone> HashGrid<N, E> {

    /// Creates a uniform grid in the N-dimensional space with the same boundaries and 
    /// periodic conditions
    pub fn generate_uniform_grid(grid:[usize; N], periodicity:[PeriodicImage; N], dims:[Float; N]) -> Self{
        let e_size = grid.iter().fold(1, |acc, x| acc * x);
        let mut hashgrid = HashGrid {
            grid,
            cells: vec![HashCell::<N, E>::new(periodicity); e_size],
            dims
        };
        for i in 0..hashgrid.cells.len() {
            hashgrid.cells[i].neighbors = hashgrid.list_combinations(hashgrid.ndim_from_1dim(i), periodicity)
                .iter()
                .map(|x| hashgrid.ndim_to_1dim(*x))
                .collect()
        }
        hashgrid
    }

    /// Creates a grid in the N-dimensional space starting from a collection of `HashCell`. The 
    /// number of cells should be equal to the expected size of the grid. It returns the `MismatchedSize` error
    /// if the cells passed do not correspond to the dimensions of the Hashgrid
    pub fn generate_from_cells(grid:[usize; N], cells:Vec<HashCell<N, E>>, dims:[Float; N]) -> Result<HashGrid<N, E>, HashGridError>{
        let expected_length = grid.iter().fold(1, |acc, x| acc * x);
        if cells.len() !=  expected_length {
            return Err(HashGridError::MismatchedSize(
                format!("Expected number of cells was {} but  {} were found", expected_length, cells.len())
            ));
        }

        let mut hashgrid = HashGrid {
            grid,
            cells,
            dims
        };

        for i in 0..hashgrid.cells.len() {
            hashgrid.cells[i].neighbors = hashgrid.list_combinations(hashgrid.ndim_from_1dim(i), hashgrid.cells[i].periodicity)
                .iter()
                .map(|x| hashgrid.ndim_to_1dim(*x))
                .collect()
        }

        Ok(hashgrid)
    }

    // CELL-LEVEL EXPLORATION

    /// Returns a reference of a slice from the cells composing the GRID
    pub fn get_cells(&self) -> &[HashCell<N, E>] {
        self.cells.as_slice()
    }

    /// Returns a mutable reference of a slice from the cells composing the GRID
    pub fn get_mut_cells(&mut self) -> &mut [HashCell<N, E>] {
        self.cells.as_mut_slice()
    }

    pub fn get_cells_index(&self) -> Vec<usize> {
        let v:Vec<usize> = (0..self.cells.len()).collect();
        v
    }

    pub fn get_cells_coords(&self) -> Vec<[usize; N]> {
        let coords:Vec<[usize; N]> = (0..self.cells.len())
                                        .into_iter()
                                        .map(|x| self.ndim_from_1dim(x))
                                        .collect();
        coords
    }

    /// Returns a referece of the elements registered under a cell with coordinates `coord` in the
    /// N-dimensional space
    pub fn get_dwellers(&self, coord:[usize; N]) -> &[E] {
        let indx = self.ndim_to_1dim(coord);
        self.cells[indx].get_dwellers()
    }

    /// Returns a mutable referece of the elements registered under a cell with coordinates `coord` in the
    /// N-dimensional space
    pub fn get_mut_dwellers(&mut self, coord:[usize; N]) -> &mut [E] {
        let indx = self.ndim_to_1dim(coord);
        self.cells[indx].get_mut_dwellers()
    }

    /// Sets the dwellers of a certain cell. It overwrites any previous registered dweller
    pub fn set_dwellers(&mut self, cell:[usize; N], dwellers:Vec<E>) {
        let indx = self.ndim_to_1dim(cell);
        self.cells[indx].set_dwellers(dwellers)
    }

    pub fn add_dweller(&mut self, cell:[usize; N], dweller:E) {
        let indx = self.ndim_to_1dim(cell);
        self.cells[indx].add_dweller(dweller)    
    }

    pub fn drop_dweller(&mut self, indx:usize, cell:[usize; N]) {
        let cell_index = self.ndim_to_1dim(cell);
        self.cells[cell_index].drop_dweller(indx)
    }

    pub fn move_dweller(&mut self, indx:usize, from:[usize; N], to:[usize; N]) {
        let from_indx = self.ndim_to_1dim(from);
        let to_indx = self.ndim_to_1dim(to);
        let dw = self.cells[from_indx].dwellers[indx].clone();
        self.cells[from_indx].drop_dweller(indx);
        self.cells[to_indx].add_dweller(dw);
    }

    pub fn get_neighbors(&self, coord:[usize; N]) -> Vec<&HashCell<N, E>> {
        let mut neighbors = Vec::new();
        let indx = self.ndim_to_1dim(coord);
        for cell_index in self.cells[indx].neighbors.iter() {
            neighbors.push(&self.cells[*cell_index])
        } 

        neighbors
    }

    pub fn get_neighbors_coords(&self, coord:[usize; N]) -> Vec<[usize; N]> {
        let mut neighbors = Vec::new();
        let indx = self.ndim_to_1dim(coord);
        for cell_index in self.cells[indx].neighbors.iter() {
            neighbors.push(self.ndim_from_1dim(*cell_index))
        } 

        neighbors
    }

    pub fn get_neighbors_dwellers(&self , coord:[usize; N]) -> Vec<&E> {
        let mut neighbors = Vec::new();
        let indx = self.ndim_to_1dim(coord);
        for cell_index in self.cells[indx].neighbors.iter() {
            neighbors.extend(self.cells[*cell_index].get_dwellers())
        } 

        neighbors
    }   

    pub fn get_all_dwellers(&self) -> Vec<&E> {
        let mut pop:Vec<&E> = Vec::new();

        for cell in self.cells.iter() {
            pop.extend(cell.get_dwellers())
        }

        pop
    }

    pub fn update_neighbors(&mut self, cell:[usize; N], periodic_images:[PeriodicImage;N]) {
        let cell_index = self.ndim_to_1dim(cell);
        self.cells[cell_index].neighbors = self.list_combinations(cell, periodic_images)
                .iter()
                .map(|x| self.ndim_to_1dim(*x))
                .collect();
    }

    pub fn population(&self) -> usize {
        self.cells.iter().fold(0, |acc, x| acc + x.population())
    }

    pub fn size(&self) -> usize {
        self.cells.len()
    }

    pub fn shape(&self) -> [usize; N] {
        self.grid
    }

    pub fn cell_center(&self, cell:[usize; N]) -> [Float; N] {
        let mut center = [0.0; N];
        for dim in 0..cell.len() {
            center[dim] = (self.dims[dim] / self.grid[dim] as Float) * (cell[dim] as Float + 0.5)

        }

        center.map(|x| x.into())
    }

    pub fn get_bounding_cell(&self, coord:[Float; N]) -> Result<[usize; N], HashGridError> {
        let mut cell:[usize; N] = [0; N];
        for c in 0..coord.len() {
            let tmp = (coord[c] / self.dims[c]).floor() as usize;
            if tmp > self.grid[c] {
                return Err(HashGridError::OutOfBounds(format!("Coordinate in position {} is out of the grid {}", c, self.dims[c])));
            }
            cell[c] = tmp
        }

        Ok(cell)
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

    fn list_combinations(&self, cell:[usize; N], periodic_images:[PeriodicImage;N]) -> Vec<[usize; N]> {
        let mut all_combs = Vec::new();
        self.list_combinations_helper(0, cell, &mut all_combs, &periodic_images);
        all_combs.remove(0);
        all_combs
    }

    fn list_combinations_helper(&self, i: usize, comb: [usize; N], all_combs: &mut Vec<[usize; N]>, periodic_images:&[PeriodicImage;N]) {
        let translations:Vec<isize> = match periodic_images[i] {
            PeriodicImage::NONE => {vec![0]},
            PeriodicImage::LEFT => {vec![0, -1]},
            PeriodicImage::RIGHT => {vec![0, 1]},
            PeriodicImage::BOTH => {vec![0, -1, 1]}
        };
        if i == N - 1 {
            for k in translations {
                let mut cell = comb.clone();

                let dim = cell[i]  as isize + k;
                if dim < 0 {
                    cell[i] = self.grid[i] - 1;
                }else if dim >= self.grid[i] as isize {
                    cell[i] = 0;
                }else {
                    cell[i] = dim as usize
                }
                all_combs.push(cell) 
            }
            return
        }
        else {
            for k in translations {
                let mut cell = comb.clone();

                let dim = cell[i]  as isize + k;
                if dim < 0 {
                    cell[i] = self.grid[i] - 1;
                }else if dim >= self.grid[i] as isize {
                    cell[i] = 0;
                }else {
                    cell[i] = dim as usize
                }

                self.list_combinations_helper(i + 1, cell, all_combs, periodic_images)
            }
        }
    }
}

impl<const N: usize, E: Clone> Index<usize> for HashGrid<N, E>{
    type Output = HashCell<N, E>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.cells[index]
    }
}

impl< const N: usize, E: Clone> Index<[usize; N]> for HashGrid<N, E>{
    type Output = HashCell<N, E>;
    fn index(&self, index: [usize; N]) -> &Self::Output {
        let indx = self.ndim_to_1dim(index);
        &self.cells[indx]
    }
}


impl< const N: usize, E: Clone> IndexMut<usize> for HashGrid<N, E>{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.cells[index]
    }
}

impl< const N: usize, E: Clone> IndexMut<[usize; N]> for HashGrid<N, E>{
    fn index_mut(&mut self, index: [usize; N]) -> &mut Self::Output {
        let indx = self.ndim_to_1dim(index);
        &mut self.cells[indx]
    }
}