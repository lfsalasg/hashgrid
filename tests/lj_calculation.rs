use rand::{thread_rng, Rng};

use hashgrid::{HashGrid, PeriodicImage};

#[derive(Clone)]
struct LJSphere {
    x:f32,
    y:f32,
    z:f32
}

#[test]
fn lj_energy() {
    let mut rng = thread_rng();
    let mut grid: HashGrid<3, LJSphere> = HashGrid::generate_uniform_grid([3,3,3], [PeriodicImage::BOTH; 3], [12.0, 12.0, 12.0]);
    for cell in grid.get_mut_cells() {
        for _ in 0..5 {
            cell.add_dweller(LJSphere { x: rng.gen_range(0.0..12.0), y: rng.gen_range(0.0..12.0), z: rng.gen_range(0.0..12.0) })
        }
    }
}