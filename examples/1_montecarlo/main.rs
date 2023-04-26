use plotters::prelude::*; //To use it the fontconfig library should be installed
                          // to install it run apt install libfontconfig1-dev
use rand::{Rng, thread_rng};
use rand_distr::{Normal, Distribution};

use hashgrid::{HashGrid, PeriodicImage, ReadGrid, WriteGrid};
use hashgrid::common::{Cardinality, Point3D};

fn pairwise_energy(epsilon:f32, squared_sig:f32, squared_distance:f32) -> f32 {
    4.0 * epsilon * ((squared_sig / squared_distance).powi(6) - (squared_sig / squared_distance).powi(3))
}

fn metropolis_criteria(new_energy:f32, old_energy:f32, temperature:f32) -> f32 {
    if new_energy < old_energy {
        1.0
    }else {
        (-1.0 * (new_energy - old_energy) / 1.0 / temperature).exp()
    }
}

fn graph(data:Vec<(f32, f32)>) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("examples/1_montecarlo/result.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Energy vs step", ("sans-serif", 50).into_font())
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f32..70000f32, -10000f32..10000f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            data,
            &RED,
        ))?
        .label("Average energy")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

#[derive(Clone, Copy, Debug)]
struct LJSphere {
    x:f32,
    y:f32,
    z:f32,
}

impl Cardinality<3> for LJSphere {
    fn coord(&self) -> Point3D {
        Point3D::new([self.x, self.y, self.z])
    }
}

fn main() {
    // Create the hashgrid and populate it with 300 spheres. We use a 27.0 A box
    // with periodic image 

    let cutoff2 = 81.0;
    let steps = 60000;

    let mut grid:HashGrid<3, LJSphere> = HashGrid::generate_uniform_grid(
        [3; 3], 
        [PeriodicImage::BOTH; 3], 
        Point3D::from_scalar(27.0)
    );

    let mut rng = thread_rng();
    let normal = Normal::new(2.7, 0.5).unwrap();

    for _ in 0..100 {
        let particle = LJSphere {
            x: rng.gen_range(0.0..27.0),
            y: rng.gen_range(0.0..27.0),
            z: rng.gen_range(0.0..27.0),
        };

        grid.add_dweller(grid.bounding_cell_coord(particle.coord()).unwrap(), particle)
    }

    // Calculate the intial energy for each element in the grid and the total grid energy
    let mut total_energy = 0.0;
    let mut average_energy = 0.0;
    let mut history:Vec<(f32, f32)> = Vec::with_capacity(steps / 100);

    for cell_index in grid.get_cells_index() {
        for dweller in grid.get_dwellers(cell_index) {
            for neighbor in grid.get_neighbors_dwellers(cell_index) {
                let d2 = dweller.coord().distance(&neighbor.coord());
                if d2 <= cutoff2 {
                    let energy = pairwise_energy(35.0, 13.69, d2);
                    total_energy += energy;
                }
            }
        }
    }

    total_energy /= 2.0;

    // Now we do the canonical Monte Carlo
    for step in 0..steps {
        if step % 100 == 0{
            history.push((step as f32, total_energy))
        }

        for cell_index in grid.get_cells_index() {
            for i in 0..grid[cell_index].population() {
                let translate = Point3D::new([
                    normal.sample(&mut rng),
                    normal.sample(&mut rng),
                    normal.sample(&mut rng)
                ]);

                let new_pos = grid[cell_index][i].coord() + translate;
                let mut current_particle_energy = 0.0;
                let mut proposed_particle_energy = 0.0;
                for neighbor in grid.get_neighbors_dwellers(cell_index) {
                    let current_d2 = neighbor.coord().distance(&grid[cell_index][i].coord());
                    let new_d2 = neighbor.coord().distance(&new_pos);
                    if current_d2 <= cutoff2 {
                        current_particle_energy += pairwise_energy(35.0, 13.69, current_d2);
                    }
                    if new_d2 <= cutoff2 {
                        
                        proposed_particle_energy += pairwise_energy(35.0, 13.69, new_d2) 
                    }
                }

                let new_energy = total_energy + proposed_particle_energy - current_particle_energy;
                let prob = metropolis_criteria(new_energy, total_energy, 200.0);
                if prob > rng.gen() {
                    total_energy = new_energy;
                    grid[cell_index][i].x = new_pos[0];
                    grid[cell_index][i].y = new_pos[1];
                    grid[cell_index][i].z = new_pos[2];
                }
            }
            
        }
        
        average_energy = (step as f32 * average_energy + total_energy) / (step + 1) as f32

    }

    graph(history).unwrap();
    println!("The average energy of the system is {}", average_energy)

}