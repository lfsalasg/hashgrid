use std::time::Instant;
use hashgrid::{HashGrid, PeriodicImage, common};

#[derive(Clone)]
struct LJSphere {
    x:f32,
    y:f32,
    z:f32
}

impl LJSphere {
    
    fn squared_distance(&self, b:&LJSphere) -> f32{
        (self.x - b.x).powi(2) + (self.y - b.y).powi(2) + (self.z - b.z).powi(2)
    }
}

impl common::Cardinality<3> for LJSphere {
    fn coord(&self) -> [common::Float; 3] {
        [self.x, self.y, self.z]
    }
}

fn pairwise_energy(epsilon:f32, squared_sig:f32, squared_distance:f32) -> f32 {
    4.0 * epsilon * ((squared_sig / squared_distance).powi(6) - (squared_sig / squared_distance).powi(3))
}

fn brute_force(p_per_side:u32, sep:f32, cutoff:f32) {
    let start = Instant::now();
    let mut particles = Vec::with_capacity(343);
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

    for i in -1..1 {
        for j in -1..1 {
            for k in -1..1 {
                images.push([i as f32 *dim, j as f32*dim, k as f32 *dim])
            }
        }
    }

    let mut total_energy = 0.0;

    for dim in images {
        for i in 0..particles.len() {
            for j in (0..i).chain(i+1..particles.len()) {
                let d = particles[i].squared_distance(&particles[j]);
                if d > squared_cutoff {
                    continue;
                }
                total_energy += pairwise_energy(5.0, 3.34, d)
            }
        }
    }
    

    println!("Energy of the system: {}. It took {:?}", total_energy / 2.0, Instant::now() - start);
}

fn with_cells(p_per_side:u32, sep:f32, cutoff:f32) {
    let start = Instant::now();
    let num_grids = (p_per_side as f32 * sep / cutoff).round() as usize;
    let squared_cutoff = cutoff.powi(2);

    let mut grid: HashGrid<3, LJSphere> = HashGrid::generate_uniform_grid([num_grids, num_grids, num_grids], [PeriodicImage::BOTH; 3], [12.0, 12.0, 12.0]);

    for i in 0..p_per_side {
        for j in 0..p_per_side {
            for k in 0..p_per_side {
                let sphere = LJSphere{x:i as f32 * sep, y:j as f32 * sep, z:k as f32 * sep};
                let coord = grid.get_bounding_cell([sphere.x, sphere.y, sphere.z]).unwrap();
                grid[coord].add_dweller(sphere);
            }
        }
    }
    
    println!("Total number of cells: {}", grid[[0, 0, 1]].get_dwellers().len());
    let mut total_energy = 0.0;
    for indx in grid.get_cells_index() {
        let dwellers = grid[indx].get_dwellers();
        let neighbors = grid.get_neighbors(indx);
        for i in 0..dwellers.len() {
            for j in (0..i).chain(i+1..dwellers.len()) {
                let d = dwellers[i].squared_distance(&dwellers[j]);
                if d > squared_cutoff {
                    continue;
                }
                total_energy += pairwise_energy(5.0, 3.34, d);
            }

            for (cell, pi) in neighbors {
                for j in cell. 
                let d = dwellers[i].squared_distance(&neighbors[j]);
                if d > squared_cutoff {
                    continue;
                }
                total_energy += pairwise_energy(5.0, 3.34, d);
            }
        }
    }
    println!("Energy of the system: {}. It took {:?}", total_energy / 2.0, Instant::now() - start);
}

fn main() {
    println!("Running brute force test");
    brute_force(8, 4.5, 12.0);
    println!("Running hashgrid test (no parallelization)");
    with_cells(8, 4.5, 12.0)
}