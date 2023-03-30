use hashgrid::{HashGrid, PeriodicImage};

#[test]

fn basic_usage() {
    let mut grid:HashGrid<f32, 3, u16> = HashGrid::generate_uniform_grid(
        [3, 3, 3], 
        [1.0, 1.0, 1.0], 
        [PeriodicImage::BOTH; 3]
    );

    let mut l = 0;
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                grid.set_dwellers([i, j, k], vec![l, l+1, l+2]);
                l += 3;
            }
        }
    }
}