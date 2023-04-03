#[cfg(test)]
mod test {
    use crate::hashgrid::{HashGrid, PeriodicImage};

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
    fn test_ndim_to_1dim() {
        let hashgrid:HashGrid<3, usize> = HashGrid::generate_uniform_grid([2, 3, 1], [PeriodicImage::BOTH; 3], [1.0, 1.0, 1.0]);
        assert_eq!(hashgrid.ndim_to_1dim([0, 0, 0]), 0);
        assert_eq!(hashgrid.ndim_to_1dim([0, 1, 0]), 1);
        assert_eq!(hashgrid.ndim_to_1dim([0, 2, 0]), 2);
        assert_eq!(hashgrid.ndim_to_1dim([1, 0, 0]), 3);
        assert_eq!(hashgrid.ndim_to_1dim([1, 1, 0]), 4);

        assert_eq!(hashgrid.ndim_from_1dim(3), [1, 0, 0]);
        assert_eq!(hashgrid.ndim_from_1dim(2), [0, 2, 0]);

        
    }

    #[test]
    fn test_set_and_get_dwellers() {
        let mut hashgrid = HashGrid::generate_uniform_grid([2, 2, 2], [PeriodicImage::BOTH; 3], [3.0, 3.0, 3.0]);
        hashgrid.set_dwellers([0,0,0], vec![0, 1, 2]);
        let dwellers = hashgrid.get_dwellers([0,0,0]);
        assert_eq!(dwellers, [0, 1, 2])
    }

    #[test]
    fn test_list_combination() {
        let hashgrid:HashGrid<2, usize> = HashGrid::generate_uniform_grid([3,3], [PeriodicImage::BOTH; 2], [3.0, 3.0]);
        let result = hashgrid.list_combinations([1, 1], [PeriodicImage::BOTH, PeriodicImage::BOTH]);
        assert_eq!(result.len(), 8);
        assert_eq!(result[0], [1, 0]);
        assert_eq!(result[result.len()-1], [2, 2]);

        let hashgrid:HashGrid<2, usize> = HashGrid::generate_uniform_grid([3,3], [PeriodicImage::BOTH; 2], [3.0, 3.0]);
        let result = hashgrid.list_combinations([2, 2], [PeriodicImage::BOTH, PeriodicImage::BOTH]);
        assert_eq!(result.len(), 8);
        assert_eq!(result[0], [2, 1]);
        assert_eq!(result[result.len()-1], [0, 0])

    }

    #[test]
    fn test_simple_cycle() {
        let mut hashgrid:HashGrid<2, usize> = HashGrid::generate_uniform_grid([3, 3], [PeriodicImage::BOTH; 2], [1.0, 1.0]);
        hashgrid.set_dwellers([1,1], vec![1, 2, 3]);
        hashgrid.set_dwellers([1,0], vec![4, 5, 6]);
        let cell = &hashgrid[[1, 1]];
        for neigbhor in &cell.neighbors {
            println!("{:?}", hashgrid.ndim_from_1dim(*neigbhor))
        }
        
        println!("{:?}", hashgrid.get_neighbors_dwellers([1,1]));

        let dwellers = hashgrid.get_mut_dwellers([1,1]);
        dwellers[2] = 8;

        println!("{:?}", hashgrid.get_dwellers([1, 1]));
    }

    #[test]
    fn test_iterate_elements () {
        let grid = simple_grid();
        let mut k:u16 = 0;
        for cell in grid.get_cells() {
            for e in cell.get_dwellers() {
            assert_eq!(*e, k);
            k += 1;
            }
        }

        k = 0;
        for e in grid.get_all_dwellers() {
            assert_eq!(*e, k);
            k += 1;
        }
    }

    #[test]
    fn test_bounding_cell() {
        let grid = simple_grid();
        assert_eq!(grid.get_bounding_cell([1.5, 0.0, 2.1]).unwrap(), [1, 0, 2])
        
    }
}