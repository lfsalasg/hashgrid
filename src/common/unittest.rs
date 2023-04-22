#[cfg(test)]
mod test {
    use crate::common::Point2D;
    #[test]
    fn point_operations() {
        let a = Point2D::new([1.0, 2.0]);
        let b = Point2D::new([5.0, -1.0]);

        assert_eq!(a + b, Point2D::new([6.0, 1.0]));
        assert_eq!(a - b, Point2D::new([-4.0, 3.0]));
        assert_eq!(a * b, Point2D::new([5.0, -2.0]));
        assert_eq!(a / b, Point2D::new([0.2, -2.0]));
        assert_eq!(a.distance(&b), 5.0)
    }
}