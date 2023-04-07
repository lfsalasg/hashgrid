use hashgrid::{HashGrid, PeriodicImage};

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

fn pairwise_energy(epsilon:f32, squared_sig:f32, squared_distance:f32) -> f32 {
    4.0 * epsilon * ((squared_sig / squared_distance).powi(6) - (squared_sig / squared_distance).powi(3))
}

fn main() {
    let mut grid: HashGrid<3, LJSphere> = HashGrid::generate_uniform_grid([3,3,3], [PeriodicImage::BOTH; 3], [12.0, 12.0, 12.0]);
    
    for i in 0..8 {
        for j in 0..8 {
            for k in 0..8 {
                let sphere = LJSphere{x:i as f32 * 4.5, y:j as f32 * 4.5, z:k as f32 * 4.5};
                let coord = grid.get_bounding_cell([sphere.x, sphere.y, sphere.z]).unwrap();
                grid[coord].add_dweller(sphere);
            }
        }
    }
    
    println!("Total number of cells: {}",grid[[0, 0, 1]].get_dwellers().len());
    let mut total_energy = 0.0;
    for indx in grid.get_cells_index() {
        let dwellers = grid[indx].get_dwellers();
        let neighbors = grid.get_neighbors_dwellers(indx);
        for i in 0..dwellers.len() {
            for j in (0..i).chain(i+1..dwellers.len()) {
                let d = dwellers[i].squared_distance(&dwellers[j]);
                total_energy += pairwise_energy(5.0, 3.34, d);
            }

            for j in 0..neighbors.len() {
                let d = dwellers[i].squared_distance(&neighbors[j]);
                total_energy += pairwise_energy(5.0, 3.34, d);
            }
        }
    }
    println!("Energy of the system: {}", total_energy / 2.0);
}