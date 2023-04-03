use rand::{thread_rng, Rng};

use hashgrid::{HashGrid, HashGridError,PeriodicImage};

#[derive(Clone)]
struct LJSphere {
    x:f32,
    y:f32,
    z:f32
}

impl LJSphere {
    
    fn squared_distance(&self, b:&LJSphere) -> f32{
        (self.x - b.x).powi(2) + (self.x - b.x).powi(2) + (self.x - b.x).powi(2)
    }
}

fn pairwise_energy(epsilon:f32, squared_sig:f32, squared_distance:f32) -> f32 {
    4.0 * epsilon * ((squared_sig / squared_distance).powi(6) - (squared_sig / squared_distance).powi(3))
}

#[test]
fn lj_energy() {
    let mut rng = thread_rng();
    let mut grid: HashGrid<3, LJSphere> = HashGrid::generate_uniform_grid([3,3,3], [PeriodicImage::BOTH; 3], [12.0, 12.0, 12.0]);
    
    for _ in 0..500 {
        let sphere = LJSphere { x: rng.gen_range(0.0..12.0), y: rng.gen_range(0.0..12.0), z: rng.gen_range(0.0..12.0) };
        
        let coord = grid.get_bounding_cell([sphere.x, sphere.y, sphere.z]).unwrap();
        grid[coord].add_dweller(sphere);
    }

    let mut total_energy = 0.0;
    for indx in grid.get_cells_index() {
        let dwellers = grid[indx].get_dwellers();
        for i in 0..dwellers.len() {
            for j in i+1..dwellers.len() {
                let d = dwellers[i].squared_distance(&dwellers[j]);
                total_energy += pairwise_energy(32.0, 0.13, d)
            }
        }
    }
}