use hashgrid::{HashGrid, PeriodicImage};

fn simple_grid() -> HashGrid<3, u16>{
    // Create and populate a grid
    let mut grid:HashGrid<3, u16> = HashGrid::generate_uniform_grid(
        [3, 3, 3],  
        [PeriodicImage::BOTH; 3],
        [1.0, 1.0, 1.0]
    );

    let mut l = 0;
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                grid.set_dwellers([i, j, k], vec![l, l+1, l+2]);
                l += 3;
            }
        }
    }

    grid
}

#[test]
fn test_populate() {
    // Create and populate a grid
    let grid = simple_grid();

    assert_eq!(grid.size(), 27);
    assert_eq!(grid.population(), 27 * 3);
    let copy_of_elements:Vec<_> = grid.get_all_dwellers().iter().map(|x| **x).collect();
    assert_eq!(copy_of_elements, (0..81).collect::<Vec<_>>())
}

#[test]
fn test_modify_elements() {
    let mut grid = simple_grid();

    grid.add_dweller([1, 0, 2], 66);
    assert_eq!(grid[[1, 0, 2]].get_dwellers()[3], 66);
    
    grid[[1, 0, 2]].drop_dweller(0);
    assert_eq!(grid[[1, 0, 2]].get_dwellers().len(), 3);

    grid.move_dweller(0, [1, 1, 1], [2, 2, 2]);
    assert_eq!(grid[[1, 1, 1]].get_dwellers().len(), 2);
    assert_eq!(grid[[2, 2, 2]].get_dwellers().len(), 4);
}