use rand::{thread_rng, Rng};

use hashgrid::{HashGrid, PeriodicImage, ReadGrid, WriteGrid};
use hashgrid::common::Point3D;

fn populate_grid() -> HashGrid<3, Point3D>{
    // Define a 3-dimensional grid that stores `Point3D` elements and is divided in 3 by 3 by 3 cell 
    // of edge 10.0. Periodicity is set for all faces of the grid
    
    let mut grid:HashGrid<3, Point3D> = HashGrid::generate_uniform_grid(
        [3, 3, 3],  
        [PeriodicImage::NONE; 3],
        Point3D::from_scalar(10.0)
    );

    // Populate the grid with Point3D elements at random location. The element is registered with the add_dweller
    // method in the cell that contains it
    let mut rng = thread_rng();
    for _ in 0..100 {
        let dweller = Point3D::new([
            rng.gen_range(0.0..10.0),
            rng.gen_range(0.0..10.0),
            rng.gen_range(0.0..10.0)
        ]);

        grid.add_dweller(grid.bounding_cell_coord(dweller).unwrap(), dweller)
    }

    assert_eq!(grid.population(), 100);

    grid
}

fn count_average_neighbors() {
    // Determine the average number of elements in a radius of 2.0 for each element
    let mut average_neighbors:f32 = 0.0;

    let grid = populate_grid();
    for cell_index in grid.get_cells_index() {
        for (i, dweller) in grid.get_dwellers(cell_index).iter().enumerate() {
            let mut true_neighbors = 0;
            for neighbor in grid.get_neighbors_dwellers(cell_index) {
                if dweller.distance(&neighbor) < 3.0 {
                    true_neighbors += 1; 
                }
            }
            average_neighbors = (i as f32 * average_neighbors + true_neighbors as f32) / (i + 1) as f32
            
        }
    }

    println!("Average neighbors are: {}", average_neighbors)
}

fn main () {
    count_average_neighbors();
}