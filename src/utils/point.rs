use std::ops::{Add, Mul, Div};

use serde::ser::{Serialize, Serializer};
use serde::de::{Deserialize, Deserializer};

use crate::common::{Cardinality, Float};

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
}

impl<const N: usize> Cardinality<N> for Point<N> {
    fn coord(&self) -> [Float; N] {
        self.0
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