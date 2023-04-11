#[cfg(test)]
mod test {
    use crate::hashgrid::{HashGrid, PeriodicImage, WriteGrid, ReadGrid};
    use crate::dynamic::IsoHashgrid;
    use crate::common::Point2D;

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
}