use std::ops::{Add, Mul, Div, Index, IndexMut, Deref};

use serde::ser::{Serialize, Serializer};
use serde::de::{Deserialize, Deserializer};

use crate::common::{Cardinality, Float};
use crate::hashgrid::HashGridError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<const N: usize>([Float; N]);

pub type Point2D = Point<2>;
pub type Point3D = Point<3>;

impl<const N:usize> Point<N> {
    pub fn new(coord:[Float; N]) -> Self {
        Self(coord)
    }
    
    pub fn from_scalar(scalar:Float) -> Self {
        Self([scalar; N])
    }

    pub fn squared_distance(&self, rhs:&Point<N>) -> Float {
        self.0.iter()
        .zip(rhs.0.iter())
        .map(|(&x, &y)| (x - y) * (x - y))
        .sum()
    }

    pub fn distance(&self, rhs:&Point<N>) -> Float {
        self.squared_distance(rhs).sqrt()
    }

    pub fn norm(&self) -> Float{
        self.distance(&Point::from_scalar(0.0))
    }

    pub fn to_vec(&self) -> Vec<Float> {
        self.0.to_vec()
    }

    pub fn from_vec(vec:Vec<Float>) -> Result<Self, HashGridError> {
        if vec.len() != N {
            return Err(HashGridError::WrongDimensionality(format!("A vector of size {} was expected but instead a vector of size {} was found", N, vec.len())))
        }

        Ok(Self(vec[0..N].try_into().unwrap()))
    }
}

impl<const N: usize> Cardinality<N> for Point<N> {
    fn coord(&self) -> Self {
        *self
    }
}

impl<const N: usize> Add for Point<N> {
    type Output = Point<N>;

    fn add(self, rhs: Point<N>) -> Point<N> {
        let mut result = [0.0; N];
        for i in 0..N {
            result[i] = self.0[i] + rhs.0[i];
        }
        Point(result)
    }
}

impl<const N: usize> Add<Float> for Point<N> {
    type Output = Self;
    fn add(self, rhs: Float) -> Self::Output {
        let mut out = self.clone();
        for i in 0..self.0.len() {
            out.0[i] += rhs
        }

        out
    }
}

impl<const N: usize> Add<Point<N>> for Float {
    type Output = Point<N>;

    fn add(self, rhs: Point<N>) -> Self::Output {
        let mut out = rhs.clone();
        for i in 0..rhs.0.len() {
            out.0[i] += self
        }

        out
    }
} 

impl<const N: usize> Mul<Point<N>> for Point<N> {
    type Output = Point<N>;

    fn mul(self, rhs: Point<N>) -> Point<N> {
        let mut result = [0.0; N];
        for i in 0..N {
            result[i] = self.0[i] * rhs.0[i];
        }
        Point(result)
    }
}

impl<const N: usize> Mul<Float> for Point<N> {
    type Output = Point<N>;

    fn mul(self, scalar: Float) -> Point<N> {
        let mut result = [0.0; N];
        for i in 0..N {
            result[i] = self.0[i] * scalar;
        }
        Point(result)
    }
}

impl<const N: usize> Mul<[isize; N]> for Point<N> {
    type Output = Point<N>;

    fn mul(self, rhs: [isize; N]) -> Point<N> {
        let mut result = [0.0; N];
        for i in 0..N {
            result[i] = self.0[i] * rhs[i] as Float;
        }
        Point(result)
    }
}

impl<const N: usize> Mul<Point<N>> for Float {
    type Output = Point<N>;

    fn mul(self, point: Point<N>) -> Point<N> {
        let mut result = [0.0; N];
        for i in 0..N {
            result[i] = point.0[i] * self;
        }
        Point(result)
    }
}

impl<const N: usize> Div<Float> for Point<N> {
    type Output = Point<N>;

    fn div(self, scalar: Float) -> Self::Output {
        self * (1.0 / scalar)
    }
}

impl<const N: usize> Div<Point<N>> for Point<N> {
    type Output = Point<N>;

    fn div(self, rhs: Point<N>) -> Point<N> {
        let mut result = [0.0; N];
        for i in 0..N {
            result[i] = self.0[i] / rhs.0[i];
        }
        Point(result)
    }
}

impl<const N: usize> Div<Point<N>> for Float {
    type Output = Point<N>;

    fn div(self, rhs: Point<N>) -> Self::Output {
        let mut result = [0.0; N];
        for i in 0..N {
            result[i] = self / rhs.0[i];
        }
        Point(result)
    }
}

impl<const N: usize> Index<usize> for Point<N> {
    type Output = Float;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<const N: usize> IndexMut<usize> for Point<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<const N: usize> Deref for Point<N> {
    type Target = [Float; N];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: usize> Serialize for Point<N> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        self.0.to_vec().serialize(serializer)
    }
}

impl<'de, const N: usize> Deserialize<'de> for Point<N> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let arr:Vec<Float> = Deserialize::deserialize(deserializer)?;
        Ok(Point(arr.try_into().unwrap()))
    }
}