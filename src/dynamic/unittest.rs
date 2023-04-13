#[cfg(test)]
mod test {
    use crate::hashgrid::{HashGrid, HashCell, PeriodicImage, WriteGrid, ReadGrid};
    use crate::dynamic::{IsoHashgrid, MultiThreaded};
    use crate::common::{Point2D, Cardinality};
    

    #[test]
    fn test_isogrid() {
        let mut grid:HashGrid<2, Point2D> = HashGrid::generate_uniform_grid(
            [3, 3], [PeriodicImage::NONE; 2], 
            Point2D::from_scalar(1.0));

        
        let mut l = 0;
        for i in 0..3 {
            for j in 0..3 {
                let p = Point2D::new([0.0, 0.0]);
                grid.set_dwellers([i, j], vec![p + l as f32, p + (l+1) as f32, p + (l+2) as f32]);
                l += 3;
            }
        } 

        let mut isogrid = IsoHashgrid::from(grid);

        isogrid.get_mut_cells()[0].dwellers[0] = Point2D::new([3.0, 0.0]);

        assert_ne!(isogrid.images[isogrid.f].get_cells()[0].dwellers[0], 
                    isogrid.images[isogrid.p].get_cells()[0].dwellers[0]);

        isogrid.commit();

        assert_eq!(isogrid.images[isogrid.p].get_cells()[0].dwellers[0], Point2D::new([3.0, 0.0]))

    }

    fn shrink(cell: &mut HashCell<2, Point2D>) {
        //let sum = cell.get_dwellers().iter().fold(0.0, |acc, x| acc + x.coord().norm());
        //println!("{}", sum);
        cell.get_mut_dwellers().iter_mut().for_each(|x| *x = 0.7 * *x);
        //sleep(Duration::from_millis(500))
        
    }

    #[test]
    fn test_split() {
        let mut grid:HashGrid<2, Point2D> = HashGrid::generate_uniform_grid(
            [3, 3], [PeriodicImage::NONE; 2], 
            Point2D::from_scalar(1.0));

        
        let mut l = 0;
        for i in 0..3 {
            for j in 0..3 {
                let p = Point2D::new([0.0, 0.0]);
                grid.set_dwellers([i, j], vec![p + l as f32, p + (l+1) as f32, p + (l+2) as f32]);
                l += 3;
            }
        } 
        let element = grid[0].get_dwellers()[2].clone();
        let mut isogrid = IsoHashgrid::from(grid);

        isogrid.split_task(2, |x| shrink(x));
        isogrid.commit();

        assert_eq!(isogrid.get_cells()[0].get_dwellers()[2], element * 0.7);
    }
}