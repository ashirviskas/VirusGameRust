//! A simplified implementation of the classic game "Breakout"

use bevy::{
    core::FixedTimestep,
    math::{const_vec2, const_vec3},
    prelude::*,
    sprite::collide_aabb::{collide, Collision},
};
use bevy_prototype_lyon::prelude::*;
use rand::prelude::*;

// Defines the amount of time that should elapse between each physics step.
const TIME_STEP: f32 = 1.0 / 60.0;

const WORLD_SIZE: f32 = 800.0;
const WORLD_SCALE: f32 = 1.0;
const MAX_FOOD: usize = 1200;

const CELL_SIZE: Vec2 = const_vec2!([20., 20.]);
const FOOD_SIZE: f32 = 8.0;
const FOOD_ENERGY: f32 = 30.0;
const DEFAULT_CELL_ENERGY: f32 = 100.0;
const ENERGY_RESOLUTION: f32 = 0.1;
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
                .with_system(apply_velocity.before(check_for_collisions))
                .with_system(loop_around),
        )
        .add_system(bevy::input::system::exit_on_esc_system)
        .add_system(energy_system)
        .add_system(energy_color_system)
        .add_system(food_dispenser)
        .run();
}

#[derive(Component)]
struct Paddle;

#[derive(Component)]
struct Particle;

#[derive(Component)]
struct Energy {
    value: f32,      // 0.0 - 100.0
    max_energy: f32, // maximum amount of energy
    waste: f32, // generated waste that needs to be cleared. If 0.0, then the cell can have max_energy, if > 0.0, then the cell can have max_energy - waste
}

impl Energy {
    fn new(max_energy: f32) -> Energy {
        Energy {
            value: max_energy,
            max_energy,
            waste: 0.0,
        }
    }
    // Add energy
    fn add_energy(&mut self, amount: f32) -> bool {
        if self.value + self.waste + ENERGY_RESOLUTION >= self.max_energy {
            false
        } else if self.value + self.waste + amount > self.max_energy {
            self.value = self.max_energy - self.waste;
            true
        } else {
            self.value += amount;
            true
        }
    }
    // Removes energy if available and adds half of it to waste. If energy is not available, returns false.
    fn remove_energy(&mut self, amount: f32) -> bool {
        if self.value >= amount {
            self.value -= amount;
            self.waste += amount / 1.5;
            true
        } else {
            false
        }
    }
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

#[derive(Component)]
struct CodonHealth {
    value: f32,
    max_value: f32,
}

enum CodonType {
    // OG first part
    None,
    Digest,
    Remove,
    Repair,
    MoveHand,
    Read,
    Write,
    // OG second part
    Food,
    Waste,
    Wall,
    WeakLoc,
    Inward,
    Outward,
    RGL,
    // Additions
    Energy,  // Can read energy levels
    LogicIf, // Will execute next codons depending on the last read state (if read num > 0.5)
}

#[derive(Component)]
struct Codon {
    type_: CodonType,
    health: CodonHealth,
    value_a: f32,
    value_b: f32,
}

#[derive(Component)]
struct Genome {
    codons: Vec<Codon>,
}

#[derive(Component)]
struct GenomeExecutor {
    genome: Genome,
    current_codon: usize,
    hand_speed: f32,
    hand_position: usize,
}

#[derive(Component, PartialEq)]
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
        if (xx / i) % 3 == 1 && (yy / i) % 3 == 1 {
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

fn energy_system(mut commands: Commands, mut query: Query<(Entity, &mut Energy)>) {
    for (entity, mut energy) in query.iter_mut() {
        energy.remove_energy(0.05);
        if energy.value - ENERGY_RESOLUTION <= 0.0 {
            // TODO: Handle despawning in the cell
            commands.entity(entity).despawn();
        }
    }
}

// for entities with energy, we will change their size depending on their energy
fn energy_color_system(mut query: Query<(&Energy, &mut Sprite, &CellType)>) {
    let red = Color::rgb(1.0, 0.0, 0.0);
    for (energy, mut sprite, cell_type) in query.iter_mut() {
        let cell_color = cell_type.to_color();
        let scale = energy.value / DEFAULT_CELL_ENERGY;
        sprite.color = red * (1.0 - scale) + cell_color * scale;
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
                    .insert(cell_type)
                    .insert(Energy {
                        value: DEFAULT_CELL_ENERGY,
                        max_energy: DEFAULT_CELL_ENERGY,
                        waste: 0.0,
                    });
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

fn loop_around(mut query: Query<(&mut Transform, &Velocity)>) {
    for (mut transform, velocity) in query.iter_mut() {
        if transform.translation.x < -WORLD_SIZE * WORLD_SCALE / 2.0 {
            transform.translation.x += WORLD_SIZE * WORLD_SCALE;
        } else if transform.translation.x > WORLD_SIZE * WORLD_SCALE / 2.0 {
            transform.translation.x -= WORLD_SIZE * WORLD_SCALE;
        }

        if transform.translation.y < -WORLD_SIZE * WORLD_SCALE / 2.0 {
            transform.translation.y += WORLD_SIZE * WORLD_SCALE;
        } else if transform.translation.y > WORLD_SIZE * WORLD_SCALE / 2.0 {
            transform.translation.y -= WORLD_SIZE * WORLD_SCALE;
        }
    }
}

fn check_for_collisions(
    mut commands: Commands,
    mut particle_query: Query<(Entity, &mut Velocity, &Transform, Option<&Food>), With<Particle>>,
    mut collider_query: Query<
        (&Transform, Option<&mut Energy>),
        (With<Collider>, Without<Particle>),
    >,
    mut collision_events: EventWriter<CollisionEvent>,
) {
    if particle_query.iter().count() == 0 {
        return;
    }
    // let (mut ball_velocity, ball_transform) = particle_query.single_mut();
    // let ball_size = ball_transform.scale.truncate();
    let particle_size: Vec2 = const_vec2!([FOOD_SIZE, FOOD_SIZE]);
    for (particle_entity, mut particle_velocity, particle_transform, maybe_food) in
        particle_query.iter_mut()
    {
        // check collision with walls
        for (transform, mut maybe_energy) in collider_query.iter_mut() {
            let collision = collide(
                particle_transform.translation,
                particle_size,
                transform.translation,
                transform.scale.truncate(),
            );
            if let Some(collision) = collision {
                // Sends a collision event so that other systems can react to the collision
                collision_events.send_default();

                if maybe_food.is_some() {
                    if maybe_energy.is_some() {
                        let cell_energy = maybe_energy.as_mut().unwrap();
                        if cell_energy.add_energy(FOOD_ENERGY) {
                            commands.entity(particle_entity).despawn();
                            break;
                        }
                    }
                }

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

fn food_dispenser(
    mut commands: Commands,
    food_query: Query<(&Particle, &Food)>,
    cells_query: Query<&Transform, (With<Collider>, Without<Food>)>,
) {
    // Check if we have enough food
    let center_vec = Vec2::new(0., 0.);
    let food_size_vec = Vec2::new(FOOD_SIZE, FOOD_SIZE);
    if food_query.iter().count() < MAX_FOOD {
        let mut random_food_position: Vec2 = Vec2::new(0., 0.);
        while true {
            let mut found_empty_spot = true;
            random_food_position = Vec2::new(
                random::<f32>() * WORLD_SIZE * WORLD_SCALE - WORLD_SIZE * WORLD_SCALE / 2.0,
                random::<f32>() * WORLD_SIZE * WORLD_SCALE - WORLD_SIZE * WORLD_SCALE / 2.0,
            );
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
        let shape = shapes::Circle {
            radius: FOOD_SIZE / 2.0,
            center: center_vec,
        };
        let random_direction =
            Vec2::new(random::<f32>() * 50.0 - 25.0, random::<f32>() * 50.0 - 25.0);

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
