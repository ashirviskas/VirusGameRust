//! A simplified implementation of the classic game "Breakout"

use bevy::{
    core::FixedTimestep,
    math::{const_vec2, const_vec3},
    prelude::*,
    sprite::collide_aabb::{collide, Collision},
};
use rand::prelude::*;
use bevy_prototype_lyon::prelude::*;

// Defines the amount of time that should elapse between each physics step.
const TIME_STEP: f32 = 1.0 / 60.0;

const WORLD_SIZE: f32 = 800.0;
const WORLD_SCALE: f32 = 1.0;
const MAX_FOOD: usize = 1200;

const BOUNDARY_THICKNESS: f32 = 10.0;

const CELL_SIZE: Vec2 = const_vec2!([40., 40.]);
const FOOD_SIZE: f32 = 8.0;
const FOOD_ENERGY: f32 = 10.0;
const DEFAULT_CELL_ENERGY: f32 = 100.0;
// These values are exact
const GAP_BETWEEN_CELLS: f32 = 1.0;

const BACKGROUND_COLOR: Color = Color::rgb(0.9, 0.9, 0.9);
const WALL_COLOR: Color = Color::rgb(0.2, 0.2, 0.2);
const FOOD_COLOR: Color = Color::rgb(0.9, 0.1, 0.2);

fn main() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(ShapePlugin)
        .insert_resource(ClearColor(BACKGROUND_COLOR))
        .add_startup_system(setup)
        .add_event::<CollisionEvent>()
        .add_system_set(
            SystemSet::new()
                .with_run_criteria(FixedTimestep::step(TIME_STEP as f64))
                .with_system(check_for_collisions)
                // .with_system(move_paddle.before(check_for_collisions))
                .with_system(apply_velocity.before(check_for_collisions)),
        )
        .add_system(bevy::input::system::exit_on_esc_system)
        .add_system(energy_system)
        .add_system(energy_size_system)
        .add_system(food_dispenser)
        .run();
}

#[derive(Component)]
struct Paddle;

#[derive(Component)]
struct Particle;

#[derive(Component)]
struct Energy {
    value: f32,
}

#[derive(Component, Deref, DerefMut)]
struct Velocity(Vec2);

#[derive(Component)]
struct Collider;

#[derive(Default)]
struct CollisionEvent;

#[derive(Component)]
struct Looptarounder;

#[derive(Default)]
struct LooptaroundEvent;

#[derive(Component)]
struct Cell;

#[derive(Component)]
struct Wall;

#[derive(Component)]
struct Food;

/// Which side of the arena is this wall located on?
enum WallLocation {
    Left,
    Right,
    Bottom,
    Top,
}

#[derive(PartialEq)]
enum CellType {
    Empty,
    Wall,
    Cell,
}

impl CellType {
    fn to_color(&self) -> Color {
        match self {
            CellType::Empty => Color::rgb(1.0, 1.0, 1.0),
            CellType::Wall => WALL_COLOR,
            CellType::Cell => Color::rgb(0.0, 0.1, 1.0),
        }
    }
}
// Will check if the brick should be spawned. Used for fractal formation
fn should_spawn_type(x: usize, y: usize, world_size: usize) -> usize {
    // array of 4 integers:
    let funny = [0, 1, 1, 2];
    let mut xx = (x / 4) * 3;
    let mut yy = (y / 4) * 3;

    xx = xx + funny[(x % 4) as usize];
    yy = yy + funny[(y % 4) as usize];

    let mut result = 0;
    let mut i = 1;
    while i < world_size {
        if (xx/i)%3 == 1 && (yy/i)%3 == 1 {
            result = 1;
            let x_part = xx % i;
            let y_part = yy % i;
            let left = x_part == 0;
            let right = x_part == i - 1;
            let top = y_part == 0;
            let bottom = y_part == i - 1;
            if left || right || top || bottom {
                result = 2;
            }
        }
        i *= 3; 
    }
    return result;
}

fn energy_system(mut query: Query<&mut Energy>){
    for mut energy in query.iter_mut() {
        if energy.value > 0.0 {
            energy.value -= 0.05;
        // println!("Energy: {}", energy.value);
        }
    }
}

// for entities with energy, we will change their size depending on their energy
fn energy_size_system(mut query: Query<(&Energy, &mut Transform)>) {
    for (energy, mut transform) in query.iter_mut() {
        let scale_x = energy.value * CELL_SIZE.x / DEFAULT_CELL_ENERGY;
        let scale_y = energy.value * CELL_SIZE.y / DEFAULT_CELL_ENERGY;
        transform.scale = Vec3::new(scale_x, scale_y, 1.0);
    }
}


// Add the game's entities to our world
fn setup(mut commands: Commands, asset_server: Res<AssetServer>) {
    // Cameras
    commands.spawn_bundle(OrthographicCameraBundle::new_2d());
    commands.spawn_bundle(UiCameraBundle::default());

    // Cells
    // Negative scales result in flipped sprites / meshes,
    // which is definitely not what we want here
    assert!(CELL_SIZE.x > 0.0);
    assert!(CELL_SIZE.y > 0.0);

    // Given the space available, compute how many rows and columns of bricks we can fit
    let n_columns = (WORLD_SIZE * WORLD_SCALE / (CELL_SIZE.x + GAP_BETWEEN_CELLS)).floor() as usize;
    let n_rows = (WORLD_SIZE * WORLD_SCALE / (CELL_SIZE.y + GAP_BETWEEN_CELLS)).floor() as usize;
    let n_vertical_gaps = n_columns - 1;

    // In Bevy, the `translation` of an entity describes the center point,
    // not its bottom-left corner
    let offset_x = CELL_SIZE.x - (WORLD_SIZE * WORLD_SCALE / 2.);
    let offset_y = CELL_SIZE.y - WORLD_SIZE * WORLD_SCALE / 2.;

    for row in 0..n_rows {
        for column in 0..n_columns {
            let cell_idx = should_spawn_type(row, column, n_rows);
            let cell_type = match cell_idx {
                0 => CellType::Empty,
                1 => CellType::Wall,
                _ => CellType::Cell,

            };
            let cell_position = Vec2::new(
                offset_x + column as f32 * (CELL_SIZE.x + GAP_BETWEEN_CELLS),
                offset_y + row as f32 * (CELL_SIZE.y + GAP_BETWEEN_CELLS),
            );
                


            // cell
            if cell_type == CellType::Cell {
                commands
                    .spawn()
                    .insert(Cell)
                    .insert_bundle(SpriteBundle {
                        sprite: Sprite {
                            color: cell_type.to_color(),
                            ..default()
                        },
                        transform: Transform {
                            translation: cell_position.extend(0.0),
                            scale: Vec3::new(CELL_SIZE.x, CELL_SIZE.y, 1.0),
                            ..default()
                        },
                        ..default()
                    })
                    .insert(Collider)
                    .insert(Cell)
                    .insert(Energy{value: DEFAULT_CELL_ENERGY});
            } else if cell_type == CellType::Wall {
                commands
                    .spawn()
                    .insert_bundle(SpriteBundle {
                        sprite: Sprite {
                            color: cell_type.to_color(),
                            ..default()
                        },
                        transform: Transform {
                            translation: cell_position.extend(0.0),
                            scale: Vec3::new(CELL_SIZE.x, CELL_SIZE.y, 1.0),
                            ..default()
                        },
                        ..default()
                    })
                    .insert(Collider)
                    .insert(Wall);
            }
        }
    }
}

fn apply_velocity(mut query: Query<(&mut Transform, &Velocity)>) {
    for (mut transform, velocity) in query.iter_mut() {
        transform.translation.x += velocity.x * TIME_STEP;
        transform.translation.y += velocity.y * TIME_STEP;
    }
}

fn check_for_collisions(
    mut commands: Commands,
    mut particle_query: Query<(&mut Velocity, &Transform), With<Particle>>,
    collider_query: Query<&Transform, (With<Collider>, Without<Particle>)>,
    mut collision_events: EventWriter<CollisionEvent>,

) {
    if particle_query.iter().count() == 0 {
        return;
    }
    // let (mut ball_velocity, ball_transform) = particle_query.single_mut();
    // let ball_size = ball_transform.scale.truncate();
    let particle_size: Vec2 = const_vec2!([FOOD_SIZE, FOOD_SIZE]);
    for (mut particle_velocity, particle_transform) in particle_query.iter_mut() {
    // check collision with walls
        for (transform) in collider_query.iter() {
            let collision = collide(
                particle_transform.translation,
                particle_size,
                transform.translation,
                transform.scale.truncate(),
            );
            if let Some(collision) = collision {
                // Sends a collision event so that other systems can react to the collision
                collision_events.send_default();

                // if maybe_food.is_some() {
                //     if maybe_cell.is_some(){
                //         let mut cell = maybe_cell.as_mut().unwrap();
                //         // cell.value += FOOD_ENERGY;
                //         commands.entity(collider_entity).despawn();
                    
                //     }
                // }


                // reflect the ball when it collides
                let mut collide_x = false;
                let mut collide_y = false;

                // only reflect if the ball's velocity is going in the opposite direction of the
                // collision
                match collision {
                    Collision::Left => collide_x = particle_velocity.x > 0.0,
                    Collision::Right => collide_x = particle_velocity.x < 0.0,
                    Collision::Top => collide_y = particle_velocity.y < 0.0,
                    Collision::Bottom => collide_y = particle_velocity.y > 0.0,
                    Collision::Inside => { /* do nothing */ }
                }

                if collide_x {
                    particle_velocity.x = -particle_velocity.x;
                }

                // reflect velocity on the y-axis if we hit something on the y-axis
                if collide_y {
                    particle_velocity.y = -particle_velocity.y;
                }
                if collide_y || collide_x {
                    break;
                }
            }
        }
    }
}

fn food_dispenser(mut commands: Commands, food_query: Query<(&Particle, &Food)>, cells_query: Query<&Transform, (With<Collider>, Without<Food>)>) {
    // Check if we have enough food
    let center_vec = Vec2::new(0., 0.);
    let food_size_vec = Vec2::new(FOOD_SIZE, FOOD_SIZE);
    if food_query.iter().count() < MAX_FOOD {
        let mut random_food_position: Vec2 = Vec2::new(0., 0.);
        while true {
            let mut found_empty_spot = true;
            random_food_position = Vec2::new(random::<f32>() * WORLD_SIZE * WORLD_SCALE - WORLD_SIZE * WORLD_SCALE / 2.0, random::<f32>() * WORLD_SIZE * WORLD_SCALE - WORLD_SIZE * WORLD_SCALE / 2.0);
            for cell_transform in cells_query.iter() {
                let collision = collide(
                    random_food_position.extend(1.0),
                    food_size_vec,
                    cell_transform.translation,
                    cell_transform.scale.truncate(),
                );
                if let Some(collision) = collision {
                    // printing collision details
                    // println!("Collision: {:?}", random_food_position);
                    found_empty_spot = false;
                    break;
                }          
            }
            if found_empty_spot {
                // println!("Found empty spot: {:?}", random_food_position);
                break;
            }
        }
        let shape = shapes::Circle{
            radius: FOOD_SIZE / 2.0,
            center: center_vec,
        };
        let random_direction = Vec2::new(random::<f32>() * 50.0 - 25.0, random::<f32>() * 50.0 - 25.0);

        commands
        .spawn_bundle(GeometryBuilder::build_as(
            &shape,
            DrawMode::Outlined {
                fill_mode: FillMode::color(FOOD_COLOR),
                outline_mode: StrokeMode::new(Color::BLACK, 0.),
            },
            Transform {
                translation: random_food_position.extend(0.0),
                scale: Vec3::new(1.0, 1.0, 1.0),
                ..default()
            },
        ))
        .insert(Food)
        .insert(Collider)
        .insert(Velocity(random_direction))
        .insert(Particle);
    }
}