#[cfg(test)]
use std::sync::Arc;
use std::thread;

use hashgrid::{HashGrid, PeriodicImage, ReadGrid, WriteGrid};
use hashgrid::common::{Point2D};

#[test]
fn cocurrent_read() {
    let mut grid:HashGrid<2, Point2D> = HashGrid::generate_uniform_grid(
        [2, 2], 
        [PeriodicImage::BOTH, PeriodicImage::BOTH], 
        Point2D::new([1.0, 1.0]));


    let mut l = 0;
    for i in 0..2 {
        for j in 0..2 {
            grid.set_dwellers([i, j], vec![
                Point2D::from_scalar(l as f32), 
                Point2D::from_scalar((l+1) as f32), 
                Point2D::from_scalar((l+2) as f32)]);
            l += 3;
        }
    }

    let shared_grid = Arc::new(grid);

    let mut handlers = Vec::new();

    for i in 0..3 {
        let grid_pointer = shared_grid.clone();
        let handler = thread::spawn(move || {
            let cell = &grid_pointer[i];
            let sum = cell.dwellers.iter().fold(0.0, |acc, x| acc + x.norm());
            println!("{}", sum)
        });
        handlers.push(handler)
    }

    for handle in handlers {
        handle.join().unwrap()
    }
}