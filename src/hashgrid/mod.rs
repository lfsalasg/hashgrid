mod unittest;

use std::ops::{Index};

pub enum HashGridError {
    MismatchedSize(String)
}

#[derive(Clone, Copy)]
pub enum PeriodicImage {
    NONE,
    BOTH,
    LEFT,
    RIGHT
}

pub struct HashGrid<T: Into<f64> + Clone, const N:usize, E: Clone> {
    grid: [usize; N],
    cells:Vec<HashCell<T, N, E>>,
}

#[derive(Clone)]
pub struct HashCell<T: Into<f64> + Clone, const N:usize, E: Clone> {
    dwellers: Vec<E>,
    neighbors: Vec<usize>,
    boundaries:[T; N], 
    periodicity: [PeriodicImage; N] 
}

impl<T: Into<f64> + Clone, const N: usize, E: Clone> HashGrid<T, N, E> {
    /// Creates a uniform grid in the N-dimensional space with the same boundaries and 
    /// periodic conditions
    pub fn generate_uniform_grid(grid:[usize; N], boundaries:[T; N], periodicity:[PeriodicImage; N]) -> Self{
        let e_size = grid.iter().fold(1, |acc, x| acc * x);
        let mut hashgrid = HashGrid {
            grid,
            cells: vec![HashCell::<T, N, E>::new(boundaries, periodicity); e_size],
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
    /// number of cells should be equal to the expected size of the grid.
    pub fn generate_from_cells(grid:[usize; N], cells:Vec<HashCell<T, N, E>>) -> Result<HashGrid<T, N, E>, HashGridError>{
        let expected_length = grid.iter().fold(1, |acc, x| acc * x);
        if cells.len() !=  expected_length {
            return Err(HashGridError::MismatchedSize(
                format!("Expected number of cells was {} but  {} were found", expected_length, cells.len())
            ));
        }

        let mut hashgrid = HashGrid {
            grid,
            cells,
        };

        for i in 0..hashgrid.cells.len() {
            hashgrid.cells[i].neighbors = hashgrid.list_combinations(hashgrid.ndim_from_1dim(i), hashgrid.cells[i].periodicity)
                .iter()
                .map(|x| hashgrid.ndim_to_1dim(*x))
                .collect()
        }

        Ok(hashgrid)
    }

    /// Returns a referece of the elements registered under a cell with coordinates `coord` in the
    /// N-dimensional space
    pub fn get_dwellers(&self, coord:[usize; N]) -> &[E] {
        let indx = self.ndim_to_1dim(coord);
        self.cells[indx].dwellers.as_slice()
    }

    /// Returns a mutable referece of the elements registered under a cell with coordinates `coord` in the
    /// N-dimensional space
    pub fn get_dwellers_mut(&mut self, coord:[usize; N]) -> &mut [E] {
        let indx = self.ndim_to_1dim(coord);
        self.cells[indx].dwellers.as_mut_slice()
    }

    pub fn set_dwellers(&mut self, cell:[usize; N], dwellers:Vec<E>) {
        let indx = self.ndim_to_1dim(cell);
        self.cells[indx].dwellers = dwellers
    }

    pub fn get_neighbors(&self, coord:[usize; N]) -> Vec<&HashCell<T, N, E>> {
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
            neighbors.extend(&self.cells[*cell_index].dwellers)
        } 

        neighbors
    }

    pub fn move_dweller(&mut self, indx:usize, from:[usize; N], to:[usize; N]) {
        let from_indx = self.ndim_to_1dim(from);
        let to_indx = self.ndim_to_1dim(to);
        let dw = self.cells[from_indx].dwellers[indx].clone();
        self.cells[from_indx].dwellers.remove(indx);
        self.cells[to_indx].dwellers.push(dw);
    }

    pub fn add_dweller(&mut self, d:E, cell:[usize; N]) {
        let indx = self.ndim_to_1dim(cell);
        self.cells[indx].dwellers.push(d);
    }

    pub fn drop_dweller(&mut self, indx:usize, cell:[usize; N]) {
        let cell_index = self.ndim_to_1dim(cell);
        self.cells[cell_index].dwellers.remove(indx);
    }

    pub fn update_neighbors(&mut self, cell:[usize; N], periodic_images:[PeriodicImage;N]) {
        let cell_index = self.ndim_to_1dim(cell);
        self.cells[cell_index].neighbors = self.list_combinations(cell, periodic_images)
                .iter()
                .map(|x| self.ndim_to_1dim(*x))
                .collect();
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

impl<T: Into<f64> + Clone, const N: usize, E: Clone> Index<usize> for HashGrid<T, N, E>{
    type Output = HashCell<T, N, E>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.cells[index]
    }
}

impl<T: Into<f64> + Clone, const N: usize, E: Clone> Index<[usize; N]> for HashGrid<T, N, E>{
    type Output = HashCell<T, N, E>;
    fn index(&self, index: [usize; N]) -> &Self::Output {
        let indx = self.ndim_to_1dim(index);
        &self.cells[indx]
    }
}

impl<T: Into<f64> + Clone, const N:usize, E: Clone> HashCell<T, N, E> {
    pub fn new (boundaries:[T; N], periodicity:[PeriodicImage; N]) -> Self{
        Self {
            dwellers: Vec::new(),
            neighbors: Vec::new(),
            boundaries,
            periodicity
        }
    }

    
}