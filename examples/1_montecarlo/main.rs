use plotters::prelude::*; //To use this the fontconfig library should be installed
                          // to install it run apt install libfontconfig1-dev
use rand::{Rng, thread_rng};
use rand_distr::{Normal, Distribution};
use clap::Parser;

use hashgrid::{HashGrid, HashCell, PeriodicImage, ReadGrid, WriteGrid};
use hashgrid::common::{Cardinality, Point3D};

fn pairwise_energy(epsilon:f32, squared_sig:f32, squared_distance:f32) -> f32 {
    4.0 * epsilon * ((squared_sig / squared_distance).powi(6) - (squared_sig / squared_distance).powi(3))
}

fn metropolis_criteria(new_energy:f32, old_energy:f32, temperature:f32) -> f32 {
    if new_energy <= old_energy {
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
        .build_cartesian_2d(0f32..70000f32, -20000f32..10000f32)?;

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

    fn set_coord(&mut self, coord:Point3D) {
        self.x = coord[0];
        self.y = coord[1];
        self.z = coord[2];
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long, default_value_t = 250.0)]
    temperature: f32,
    #[arg(short, long, default_value_t = 36.0)]
    box_size: f32,
    #[arg(short, long, default_value_t = 100)]
    n_particles: i32,
    #[arg(short, long, default_value_t = 12.0)]
    cutoff: f32,
    #[arg(short, long, default_value_t = 50000)]
    steps: usize
}

fn calc_cell_energy(cell:&HashCell<3, LJSphere>) -> f32 {
    let mut cell_energy = 0.0;
    for i in 0..cell.population() {
        for j in (0..i).chain(i+1..cell.population()) {
            let d2 = cell[i].coord().squared_distance(&cell[j].coord());
            cell_energy += pairwise_energy(158.5, 13.838, d2);
        }
    }

    cell_energy
}

fn calc_box_energy(grid: &HashGrid<3, LJSphere>, cutoff2:f32) -> f32 {
    let mut total_energy = 0.0;
    for cell_index in grid.get_cells_index() {
        total_energy += calc_cell_energy(&grid[cell_index]);
        for dweller in grid.get_dwellers(cell_index) {
            for (coord, image) in grid.get_neighbors_coords(cell_index) {
                for neighbor in grid.get_dwellers(coord) {
                    let d2 = dweller.coord().squared_distance(&(neighbor.coord() + grid.dims() * image));
                    if d2 <= cutoff2 {
                        let energy = pairwise_energy(158.5, 13.838, d2);
                        total_energy += energy;
                    }
                }
            }
        }
    }

    total_energy / 2.0
}

fn main() {
    // Create the hashgrid and populate it with 300 spheres. We use a 27.0 A box
    // with periodic image 
    let args = Args::parse();

    let cutoff2 = args.cutoff.powi(2);
    if args.cutoff > 2.0 * args.box_size {
        panic!("Irrecoverable error: The box size must be at least two times the cutoff")
    }
    
    let n_grids = (args.box_size / args.cutoff) as usize;
    println!("Creating hashgrid. Grid shape {} x {} x {}", n_grids, n_grids, n_grids);


    let mut grid:HashGrid<3, LJSphere> = HashGrid::generate_uniform_grid(
        [n_grids; 3], 
        [PeriodicImage::BOTH; 3], 
        Point3D::from_scalar(args.box_size)
    );

    let mut rng = thread_rng();
    let normal = Normal::new(args.box_size*0.05, 0.05).unwrap();

    for _ in 0..args.n_particles {
        let particle = LJSphere {
            x: rng.gen_range(0.0..args.box_size),
            y: rng.gen_range(0.0..args.box_size),
            z: rng.gen_range(0.0..args.box_size),
        };

        grid.add_dweller(grid.bounding_cell_coord(particle.coord()).unwrap(), particle)
    }

    // Calculate the intial energy for each element in the grid and the total grid energy
    let mut average_energy = 0.0;
    let mut history:Vec<(f32, f32)> = Vec::with_capacity(args.steps / 100);
    let mut movements = (0, 0); // total vs accepted

    let mut total_energy = calc_box_energy(&grid, cutoff2);
    println!("Initial energy: {}", total_energy);

    // Now we do the canonical Monte Carlo
    for step in 0..args.steps {
        if step % 50 == 0{
            history.push((step as f32, total_energy));
            total_energy = calc_box_energy(&grid, cutoff2)
        }
        
        for cell_index in grid.get_cells_index() {
            let mut elements_to_move = Vec::with_capacity(grid[cell_index].population());
            for i in 0..grid[cell_index].population() {
                movements.0 += 1;
                let translate = Point3D::new([
                    normal.sample(&mut rng),
                    normal.sample(&mut rng),
                    normal.sample(&mut rng)
                ]);


                let new_pos = grid[cell_index][i].coord() + translate;
                let mut current_particle_energy = 0.0;
                let mut proposed_particle_energy = 0.0;
                
                for j in (0..i).chain(i+1..grid[cell_index].population()) {
                    let current_d2 = grid[cell_index][i].coord().squared_distance(&grid[cell_index][j].coord());
                    let new_d2 = new_pos.coord().squared_distance(&grid[cell_index][j].coord());
                    if new_d2 <= cutoff2 {
                        proposed_particle_energy += pairwise_energy(158.5, 13.838, new_d2);
                    }
                    if current_d2 <= cutoff2 {
                        current_particle_energy += pairwise_energy(158.5, 13.838, current_d2);
                    }
                       
                }

                for (coord, image) in grid.get_neighbors_coords(cell_index) {
                    let displace = grid.dims() * image;
                    
                    for neighbor in grid.get_dwellers(coord) {
                        let current_d2 = grid[cell_index][i].coord().squared_distance(&(neighbor.coord() + displace));
                        let new_d2 = new_pos.coord().squared_distance(&(neighbor.coord() + displace));
                        
                        if current_d2 <= cutoff2 {
                            current_particle_energy += pairwise_energy(158.5, 13.838, current_d2);
                        }
                        if new_d2 <= cutoff2 {
                            proposed_particle_energy += pairwise_energy(158.5, 13.838, new_d2);   
                        }

                        
                    }

                }
                    
                let new_energy = total_energy + proposed_particle_energy - current_particle_energy;
               
                let prob = metropolis_criteria(new_energy, total_energy, args.temperature);
                if prob > rng.gen() {

                    grid[cell_index][i].set_coord(new_pos);
                    total_energy = new_energy;
                    elements_to_move.push(i);
                    movements.1 += 1;
                }
                
            }
            let elements = grid.purge(cell_index, &mut elements_to_move);
            grid.try_allocate(elements).expect("Fatal error: Elements out of bounds");
            
        }

        
        
        //average_energy = (step as f32 * average_energy + total_energy) / (step + 1) as f32
        average_energy = (average_energy + total_energy) / 2.0;

    }
    //println!("{:?}",history);
    graph(history).unwrap();
    
    //println!("{}", total_energy);

    //println!("{}", calc_box_energy(&grid, cutoff2));
    println!("The average energy of the system is {}", average_energy);
    println!("Number of accepted movements: {} / {}", movements.1, movements.0)

}