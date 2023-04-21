#[cfg(test)]
mod test {
    use crate::hashgrid::{HashGrid, PeriodicImage, ReadGrid};
    use crate::utils::relocate_dwellers;
    use crate::common::{Point2D, Cardinality};

    #[test]
    fn test_relocate_dwellers() {
        let mut grid:HashGrid<2, Point2D> = HashGrid::generate_uniform_grid(
            [3, 3], 
            [PeriodicImage::BOTH; 2], 
            Point2D::from_scalar(3.0)
        );

        grid[[1, 1]].add_dweller(Point2D::new([1.5, 1.5]));
        grid[[0, 1]].add_dweller(Point2D::new([0.3, 0.1]));
        grid[[0, 0]].add_dweller(Point2D::new([3.7, 2.5]));

        relocate_dwellers(&mut grid);

        let dwellers = grid.get_all_dwellers();
        assert_eq!(grid.bounding_cell_coord(dwellers[0].coord()).unwrap(), [0, 0]);
        assert_eq!(grid.bounding_cell_coord(dwellers[1].coord()).unwrap(), [0, 2]);
        assert_eq!(grid.bounding_cell_coord(dwellers[2].coord()).unwrap(), [1, 1]);
    }
}