#[cfg(test)]
mod test {
    use crate::hashgrid::{HashGrid, PeriodicImage};

    #[test]
    fn test_ndim_to_1dim() {
        let hashgrid:HashGrid<f32, 3, usize> = HashGrid::generate_uniform_grid([2, 3, 1], [PeriodicImage::BOTH; 3], [1.0, 1.0, 1.0]);
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
        let hashgrid:HashGrid<f32, 2, usize> = HashGrid::generate_uniform_grid([3,3], [PeriodicImage::BOTH; 2], [3.0, 3.0]);
        let result = hashgrid.list_combinations([1, 1], [PeriodicImage::BOTH, PeriodicImage::BOTH]);
        assert_eq!(result.len(), 8);
        assert_eq!(result[0], [1, 0]);
        assert_eq!(result[result.len()-1], [2, 2]);

        let hashgrid:HashGrid<f32, 2, usize> = HashGrid::generate_uniform_grid([3,3], [PeriodicImage::BOTH; 2], [3.0, 3.0]);
        let result = hashgrid.list_combinations([2, 2], [PeriodicImage::BOTH, PeriodicImage::BOTH]);
        assert_eq!(result.len(), 8);
        assert_eq!(result[0], [2, 1]);
        assert_eq!(result[result.len()-1], [0, 0])

    }

    #[test]
    fn test_simple_cycle() {
        let mut hashgrid:HashGrid<f32, 2, usize> = HashGrid::generate_uniform_grid([3, 3], [PeriodicImage::BOTH; 2], [1.0, 1.0]);
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
    fn test_into_iter() {
        let hashgrid:HashGrid<f32, 2, usize> = HashGrid::generate_uniform_grid(
            [3, 3], [PeriodicImage::BOTH; 2], [1.0, 1.0]);
        
        assert!(hashgrid.into_iter().all(|x| x.periodicity.len() == 2));
    }
}