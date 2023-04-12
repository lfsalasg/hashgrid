use std::time::Instant;
use std::sync::{Arc, Mutex};
use std::thread;

use hashgrid::{HashGrid, PeriodicImage, ReadGrid};
use hashgrid::common::{Cardinality, Point3D};


#[derive(Clone)]
struct LJSphere {
    x:f32,
    y:f32,
    z:f32
}

impl Cardinality<3> for LJSphere {
    fn coord(&self) -> Point3D {
        Point3D::new([self.x, self.y, self.z])
    }
}

fn pairwise_energy(epsilon:f32, squared_sig:f32, squared_distance:f32) -> f32 {
    4.0 * epsilon * ((squared_sig / squared_distance).powi(6) - (squared_sig / squared_distance).powi(3))
}

fn brute_force(p_per_side:u32, sep:f32, cutoff:f32) {
    let start = Instant::now();
    let mut particles = Vec::with_capacity(512);
    let squared_cutoff = cutoff.powi(2);
    let dim = p_per_side as f32 * sep; 
    let mut images = Vec::with_capacity(27);

    for i in 0..p_per_side {
        for j in 0..p_per_side {
            for k in 0..p_per_side {
                let sphere = LJSphere{x:i as f32 * sep, y:j as f32 * sep, z:k as f32 * sep};
                particles.push(sphere)
            }
        }
    }

    for i in -1..2 {
        for j in -1..2 {
            for k in -1..2 {
                images.push([i as f32 *dim, j as f32*dim, k as f32 *dim])
            }
        }
    }

    let mut total_energy = 0.0;

    let mut interaction_counter = 0;
    for dim in images {
        let disp = Point3D::new(dim);
        for i in 0..particles.len() {
            for j in (0..i).chain(i+1..particles.len()) {
                let d = particles[i].coord().squared_distance(&(particles[j].coord() + disp));
                if d > squared_cutoff {
                    continue;
                }
                interaction_counter += 1; 
                total_energy += pairwise_energy(5.0, 3.34, d)
            }
        }
    }
    
    println!("Total number of interactions: {}", interaction_counter);
    println!("Energy of the system: {}. It took {:?}", total_energy / 2.0, Instant::now() - start);
}

fn with_cells(p_per_side:u32, sep:f32, cutoff:f32) {
    let start = Instant::now();
    let num_grids = (p_per_side as f32 * sep / cutoff).round() as usize;
    let squared_cutoff = cutoff.powi(2);
    let dim = p_per_side as f32 * sep; 

    let mut grid: HashGrid<3, LJSphere> = HashGrid::generate_uniform_grid(
        [num_grids, num_grids, num_grids], 
        [PeriodicImage::BOTH; 3], 
        Point3D::new([cutoff, cutoff, cutoff]));

    for i in 0..p_per_side {
        for j in 0..p_per_side {
            for k in 0..p_per_side {
                let sphere = LJSphere{x:i as f32 * sep, y:j as f32 * sep, z:k as f32 * sep};
                let coord = grid.get_bounding_cell([sphere.x, sphere.y, sphere.z]).unwrap();
                grid[coord].add_dweller(sphere);
            }
        }
    }

    println!("Total number of cells: {}", grid.get_cells().len());
    let mut total_energy = 0.0;
    let mut interaction_counter = 0;
    for indx in grid.get_cells_index() {
        let dwellers = grid[indx].get_dwellers();
        let neighbors = grid.get_neighbors(indx);
        for i in 0..dwellers.len() {
            for j in (0..i).chain(i+1..dwellers.len()) {
                let d = dwellers[i].coord().squared_distance(&dwellers[j].coord());
                if d > squared_cutoff {
                    continue;
                }
                interaction_counter += 1;
                total_energy += pairwise_energy(5.0, 3.34, d);
            }

            for (cell, pi) in neighbors.iter() {
                let disp = Point3D::from_scalar(dim) * *pi;
                for j in cell.get_dwellers() {
                    let d = dwellers[i].coord().squared_distance(&(j.coord() + disp));
                    if d > squared_cutoff {
                        continue;
                    }
                    interaction_counter += 1;
                    total_energy += pairwise_energy(5.0, 3.34, d);
                }                 
            }
        }
    }
    println!("Total number of interactions: {}", interaction_counter);
    println!("Energy of the system: {}. It took {:?}", total_energy / 2.0, Instant::now() - start);
}

fn with_cells_in_parallel(p_per_side:u32, sep:f32, cutoff:f32, n_threads:usize) {
    let start = Instant::now();
    let num_grids = (p_per_side as f32 * sep / cutoff).round() as usize;
    let squared_cutoff = cutoff.powi(2);
    let dim = p_per_side as f32 * sep; 

    let mut grid: HashGrid<3, LJSphere> = HashGrid::generate_uniform_grid(
        [num_grids, num_grids, num_grids], 
        [PeriodicImage::BOTH; 3], 
        Point3D::new([cutoff, cutoff, cutoff]));

    for i in 0..p_per_side {
        for j in 0..p_per_side {
            for k in 0..p_per_side {
                let sphere = LJSphere{x:i as f32 * sep, y:j as f32 * sep, z:k as f32 * sep};
                let coord = grid.get_bounding_cell([sphere.x, sphere.y, sphere.z]).unwrap();
                grid[coord].add_dweller(sphere);
            }
        }
    }

    let chunk_size = grid.get_cells_index().len() / n_threads;

    let shared_grid = Arc::new(grid);
    let mut handlers = Vec::new();
    let mut chunk_start = 0;
    let mut chunk_end = chunk_start + chunk_size;
    let total_energy = Arc::new(Mutex::new(0.0));
    for _ in 0..n_threads {
        let grid_pointer = shared_grid.clone();
        let total_energy = Arc::clone(&total_energy);
        let handler = thread::spawn(move || {
            for j in chunk_start..chunk_end {
                let dwellers = &grid_pointer[j].get_dwellers();
                let neighbors = &grid_pointer.get_neighbors(j);
                let mut cell_energy:f32 = 0.0;
                for i in 0..dwellers.len() {
                    for j in (0..i).chain(i+1..dwellers.len()) {
                        let d = dwellers[i].coord().squared_distance(&dwellers[j].coord());
                        if d > squared_cutoff {
                            continue;
                        }
                        cell_energy += pairwise_energy(5.0, 3.34, d);
                    }
        
                    for (cell, pi) in neighbors.iter() {
                        let disp = Point3D::from_scalar(dim) * *pi;
                        for j in cell.get_dwellers() {
                            let d = dwellers[i].coord().squared_distance(&(j.coord() + disp));
                            if d > squared_cutoff {
                                continue;
                            }
                            cell_energy += pairwise_energy(5.0, 3.34, d);
                        }                 
                    }
                }
                let mut total_energy = total_energy.lock().unwrap();
                *total_energy += cell_energy;
                //println!("{}", cell_energy)
            }
        });
        handlers.push(handler);

        chunk_start = chunk_end;
        chunk_end = (chunk_end + chunk_size).min(shared_grid.get_cells_index().len())
    }

    for handle in handlers {
        handle.join().unwrap()
    }

    println!("Total energy is {}. It took {:?}", *total_energy.lock().unwrap() / 2.0, Instant::now() - start);
}


fn main() {
    println!("Running brute force test");
    brute_force(8, 4.5, 12.0);
    println!("Running hashgrid test (no parallelization)");
    with_cells(8, 4.5, 12.0);
    println!("Running hashgrid test (with parallelization)");
    with_cells_in_parallel(8, 4.5, 12.0, 3)
}