# Hashgrid

The hashgrid library is an implementation of the hash grid data structure written in pure rust. It provides a N-dimensional, agnostic grid abstraction to store and manipulate data in contiguous uniform cells. Some features of this crate include:

- N-Dimensional: You can represent any space by simply defining the dimensionality. It becomes specially useful if the implementation requires grids in 2 and 3 dimensional spaces at the same time or if you are doing some high dimensional calculations (k-neighbors for example)

- Agnostic: Store any type of data that implements the `Clone` and `Cardinality` traits.

- Memory safe: Run cocurrent code without fear using the `IsoHashGrid` struct and the `MultiThreaded` trait.

- Periodic images: Specially useful for physics, the current implementation allows to create periodic images of the grid and move/locate element through the periodic image in a very simple and convenient manner.

# Use Hashgrid

Initializing a hashgrid is very simple. You can initialize an empty grid based
on the dimensions of it or by passing an initialized vector of `HashCell`. 

```
// Define a 3-dimensional grid that stores `u16` values and is divided in 3 by 3 by 3 cell 
// of edge 1.0 each. All cells can see their periodic images on all of 
// its faces.

use hashgrid::{HashGrid, PeriodicImage, ReadGrid, WriteGrid};
use hashgrid::common::Point3D;

let mut grid:HashGrid<3, Point3D> = HashGrid::generate_uniform_grid(
        [3, 3, 3],  
        [PeriodicImage::BOTH; 3],
        Point3D::new([1.0, 1.0, 1.0])
    );

    let mut l = 0;
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                grid.set_dwellers([i, j, k], vec![Point3D::from_scalar(l as f32), 
                Point3D::from_scalar((l+1) as f32), 
                Point3D::from_scalar((l+2) as f32)
            ]);
                l += 3;
            }
        }
    }
```