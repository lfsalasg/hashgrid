#[cfg(test)]
mod test {
    use serde_json;

    use crate::common::{Point2D,Point3D};
    use crate::hashgrid::{HashGrid, HashCell, PeriodicImage, ReadGrid, WriteGrid};

    fn simple_grid() -> HashGrid<3, Point3D>{
        // Create and populate a grid
        let mut grid:HashGrid<3, Point3D> = HashGrid::generate_uniform_grid(
            [3, 3, 3],  
            [PeriodicImage::BOTH; 3],
            Point3D::new([1.0, 1.0, 1.0])
        );
    
        let mut l = 0;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let p = Point3D::new([0.0, 0.0, 0.0]);
                    grid.set_dwellers([i, j, k], vec![p + l as f32, p + (l+1) as f32, p + (l+2) as f32]);
                    l += 3;
                }
            }
        }
    
        grid
    }

    #[test]
    fn test_ndim_to_1dim() {
        let hashgrid:HashGrid<3, Point3D> = HashGrid::generate_uniform_grid([2, 3, 1], [PeriodicImage::BOTH; 3], Point3D::new([1.0, 1.0, 1.0]));
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
        let mut hashgrid = HashGrid::generate_uniform_grid([2, 2, 2], [PeriodicImage::BOTH; 3], Point3D::new([3.0, 3.0, 3.0]));
        hashgrid.set_dwellers([0,0,0], 
            vec![
                Point3D::new([0.0, 0.0, 0.0]), 
                Point3D::new([1.0, 1.0, 1.0]), 
                Point3D::new([2.0, 2.0, 2.0])]
            );
        let dwellers = hashgrid.get_dwellers([0,0,0]);
        assert_eq!(dwellers, [
            Point3D::new([0.0, 0.0, 0.0]), 
            Point3D::new([1.0, 1.0, 1.0]), 
            Point3D::new([2.0, 2.0, 2.0])
            ])
    }

    #[test]
    fn test_list_combination() {
        let hashgrid:HashGrid<2, Point2D> = HashGrid::generate_uniform_grid([3,3], [PeriodicImage::BOTH; 2], Point2D::new([3.0, 3.0]));
        let result = hashgrid.list_combinations([1, 1], [PeriodicImage::BOTH, PeriodicImage::BOTH]);
        assert_eq!(result.len(), 8);
        assert_eq!(result[0].0, [1, 0]);
        assert_eq!(result[0].1, [0, 0]);
        assert_eq!(result[result.len()-1].0, [2, 2]);

        let hashgrid:HashGrid<2, Point2D> = HashGrid::generate_uniform_grid([3,3], [PeriodicImage::BOTH; 2], Point2D::new([3.0, 3.0]));
        let result = hashgrid.list_combinations([2, 2], [PeriodicImage::BOTH, PeriodicImage::BOTH]);
        assert_eq!(result.len(), 8);
        assert_eq!(result[0].0, [2, 1]);
        assert_eq!(result[result.len()-1].0, [0, 0]);
        assert_eq!(result[result.len()-1].1, [1, 1], "{:?}", result[result.len()-1].1)

    }

    #[test]
    fn test_simple_cycle() {
        let mut hashgrid:HashGrid<2, Point2D> = HashGrid::generate_uniform_grid([3, 3], [PeriodicImage::BOTH; 2], Point2D::new([1.0, 1.0]));
        hashgrid.set_dwellers([1, 1], vec![
            Point2D::new([1.0, 1.0]),
            Point2D::new([2.0, 2.0]),
            Point2D::new([3.0, 3.0]),
        ]);

        hashgrid.set_dwellers([1, 0], vec![
            Point2D::new([4.0, 4.0]),
            Point2D::new([5.0, 5.0]),
            Point2D::new([6.0, 6.0]),
        ]);

        let cell = &hashgrid[[1, 1]];
        for neigbhor in &cell.neighbors {
            println!("{:?}", hashgrid.ndim_from_1dim(neigbhor.0))
        }
        
        println!("{:?}", hashgrid.get_neighbors_dwellers([1,1]));

        let dwellers = hashgrid.get_mut_dwellers([1,1]);
        dwellers[2] = Point2D::new([8.0, 8.0]);

        println!("{:?}", hashgrid.get_dwellers([1, 1]));
    }

    #[test]
    fn test_iterate_elements () {
        let grid = simple_grid();
        let mut k:u16 = 0;
        for cell in grid.get_cells() {
            for e in cell.get_dwellers() {
            assert_eq!(*e, Point3D::new([k as f32, k as f32, k as f32]));
            k += 1;
            }
        }

        k = 0;
        for e in grid.get_all_dwellers() {
            assert_eq!(*e, Point3D::new([k as f32, k as f32, k as f32]));
            k += 1;
        }
    }

    #[test]
    fn test_bounding_cell() {
        let grid:HashGrid<3, Point3D> = HashGrid::generate_uniform_grid(
            [3, 3, 3], 
            [PeriodicImage::BOTH, PeriodicImage::BOTH, PeriodicImage::NONE], 
            Point3D::from_scalar(1.0)
        );

        assert_eq!(grid.bounding_cell_coord(Point3D::from_scalar(0.1)).unwrap(), [0, 0, 0]);
        assert_eq!(grid.bounding_cell_coord(Point3D::from_scalar(0.5)).unwrap(), [1, 1, 1]);
        assert_eq!(grid.bounding_cell_coord(Point3D::new([1.5, 0.1, 0.1])).unwrap(), [1, 0, 0]);
        assert_eq!(grid.bounding_cell_coord(Point3D::new([-1.1, 0.1, 0.1])).unwrap(), [2, 0, 0]);
        match grid.bounding_cell_coord(Point3D::new([0.5, 0.5, 1.3])) {
            Ok(_) => {
                panic!("Should panic")
            }
            Err(_) => {}
        }
    }

    #[test]
    fn test_serialize_deserialize_hashcell() {
        // create a HashCell struct
        let hashcell = HashCell::<3, Point3D> {
            dwellers: vec![
                Point3D::new([1.0, 1.0, 1.0]),
                Point3D::new([2.0, 2.0, 2.0]),
                Point3D::new([3.0, 3.0, 3.0]),
            ],
            neighbors: vec![(0, [0,0,0]), (2, [0,0,0])],
        };

        // serialize the HashCell struct to JSON
        let json = serde_json::to_string(&hashcell).unwrap();
        println!("{}", json);
        // deserialize it back to 
        let new_cell: HashCell<3, Point3D> = serde_json::from_str(&json).unwrap();
        assert_eq!(new_cell.dwellers, hashcell.dwellers)
    }

    #[test]
    fn test_serialize_desearialize_hashgrid() {
        let grid:HashGrid<2, Point2D> = HashGrid::generate_uniform_grid(
            [3, 3], 
            [PeriodicImage::BOTH, PeriodicImage::BOTH], 
            Point2D::new([1.0, 1.0])
        );

        let json = serde_json::to_string(&grid).unwrap();
        println!("{}", json);
        let new_grid: HashGrid<2, Point2D> = serde_json::from_str(&json).unwrap();
        assert_eq!(grid.cells.len(), new_grid.cells.len())
    }

    #[test]
    fn test_anchor_and_center() {
        let grid = simple_grid();
        assert_eq!(grid.cell_anchor([0, 0, 0]), Point3D::from_scalar(0.0));
        assert_eq!(grid.cell_anchor([1, 0, 0]), Point3D::new([1.0/3.0, 0.0, 0.0]));
        
        assert_eq!(grid.cell_center([0, 0, 0]), Point3D::from_scalar(1.0/6.0));
        assert_eq!(grid.cell_center([1, 0, 0]), Point3D::new([1.0/2.0, 1.0/6.0, 1.0/6.0]))
    }
}