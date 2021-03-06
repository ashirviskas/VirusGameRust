//! A simplified implementation of the classic game "Breakout"

use bevy::diagnostic::{FrameTimeDiagnosticsPlugin, LogDiagnosticsPlugin};
use bevy::{
    core::FixedTimestep,
    input::mouse::{MouseMotion, MouseWheel},
    math::{const_vec2, const_vec3},
    prelude::*,
    sprite::collide_aabb::{collide, Collision},
    tasks::AsyncComputeTaskPool,
    utils::tracing::dispatcher,
};
use bevy_prototype_lyon::prelude::FillMode;
use bevy_prototype_lyon::prelude::*;
use bevy_rapier2d::prelude::*;
use rand::prelude::*;
use std::mem;
use std::ops::Deref;

const PI: f32 = 3.141592653589793;
// Defines the amount of time that should elapse between each physics step.
const TIME_STEP: f32 = 1.0 / 60.0;

const WORLD_SIZE: f32 = 800.0;
const WORLD_SCALE: f32 = 1.0;

const MAX_FOOD: usize = 400;
const MAX_WASTE: usize = 200;

const FOOD_SPAWN_TIMER: f32 = 0.1;
const FOOD_SPAWN_RATE: usize = 2;

const WASTE_DISPAWN_TIMER: f32 = 0.6;
const WASTE_DESPAWN_PER_TIMER: usize = 1;

const CELL_Z_LAYER: f32 = 0.0;
const CELL_WALL_Z_LAYER: f32 = 0.1;
const CODON_READER_Z_LAYER: f32 = 0.15;
const CODON_Z_LAYER: f32 = 0.2;
const UGO_Z_LAYER: f32 = 0.1;
const FOOD_Z_LAYER: f32 = 0.3;
const WASTE_Z_LAYER: f32 = 0.3;

const CELL_SIZE: Vec2 = const_vec2!([35., 35.]);
const DEFAULT_CELL_WALL_THICKNESS: f32 = 0.125;
const GAP_BETWEEN_CELLS: f32 = 0.;
const FOOD_SIZE: f32 = 8.0;
const FOOD_ENERGY: f32 = 30.0;
const DEFAULT_CELL_ENERGY: f32 = 100.0;
const ENERGY_RESOLUTION: f32 = 0.1;
const CELL_ENERGY_LOSS_RATE: f32 = 0.05;
const MAX_CELL_CODONS: usize = 45;

const CODON_SIZE: f32 = 0.1;
const CODON_RADIUS: f32 = 0.25;
const CODON_ENERGY_COST_PER_EXECUTION: f32 = 0.03; // how much energy will be taken out of the cell for executing a codon
const CODON_HEALTH_COST_PER_EXECUTION: f32 = 0.03; // how much health will be taken out of the codon for executing

const CODON_EXECUTION_RATE: f32 = 0.08; // Seconds between each codon execution for cell
const CODON_EXECUTION_RATE_VARIATION: f32 = 0.8; // How much the execution rate can vary from the base rate by %

const CODON_READ_MUTATION_RATE: f32 = 0.001; // Misread chance when reading codon
const MAX_CODON_IDX: i32 = 60;
const MIN_CODON_IDX: i32 = -60;

const UGO_SIZE: f32 = 32.0;

const BACKGROUND_COLOR: Color = Color::rgb(0.9, 0.9, 0.9);
const WALL_COLOR: Color = Color::rgb(0.2, 0.2, 0.2);
const FOOD_COLOR: Color = Color::rgb(0.9, 0.1, 0.2);
const WASTE_COLOR: Color = Color::rgb(0.5, 0.5, 0.2);
const UGO_COLOR: Color = Color::rgb(0.15, 0.1, 0.1);
const CELL_WALL_COLOR: Color = Color::rgb(0.5, 0.5, 0.5);

const GENOME_READER_COLOR: Color = Color::rgba(1.0, 1.0, 1.0, 0.8);

fn main() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(ShapePlugin)
        .add_plugin(RapierPhysicsPlugin::<NoUserData>::pixels_per_meter(100.0))
        .add_plugin(RapierDebugRenderPlugin::default())
        .add_plugin(LogDiagnosticsPlugin::default())
        .add_plugin(FrameTimeDiagnosticsPlugin::default())
        .insert_resource(ClearColor(BACKGROUND_COLOR))
        .add_startup_system(setup)
        .add_event::<CollisionEvent>()
        .add_system_set(
            SystemSet::new()
                .with_run_criteria(FixedTimestep::step(TIME_STEP as f64))
                .with_system(bevy::transform::transform_propagate_system)
                .with_system(check_for_collisions)
                .with_system(apply_velocity.before(check_for_collisions))
                .with_system(energy_system)
                .with_system(loop_around),
        )
        .add_system(bevy::input::system::exit_on_esc_system)
        .add_system(pan_orbit_camera)
        .add_system(energy_color_system)
        .insert_resource(WasteDespawnTimer(Timer::from_seconds(
            WASTE_DISPAWN_TIMER,
            true,
        )))
        .insert_resource(FoodSpawnTimer(Timer::from_seconds(FOOD_SPAWN_TIMER, true)))
        .add_system(waste_despawner)
        .add_system(food_dispenser)
        .add_system(codon_executing)
        .add_system(update_cell_genome_executor)
        .add_system(change_cell_genome)
        .add_system(update_cell_genome)
        .run();
}

#[derive(Component)]
struct Particle;
struct FoodSpawnTimer(Timer);

struct WasteDespawnTimer(Timer);

struct CodonExecutionTimer(Timer);

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

struct PhysicsGrid {}

#[derive(PartialEq, Copy, Clone)]
enum CellWallPosition {
    Left,
    Right,
    Top,
    Bottom,
}

struct CellWall {
    wall_health: f32,
}

impl CellWall {
    fn new() -> CellWall {
        CellWall { wall_health: 1.0 }
    }
    fn get_wall_thickness(&self) -> f32 {
        DEFAULT_CELL_WALL_THICKNESS * self.wall_health
    }
    fn get_wall_position_part(&self, cell_wall_position: &CellWallPosition) -> Vec2 {
        match cell_wall_position {
            CellWallPosition::Left => Vec2::new(-0.6 + self.get_wall_thickness() / 2., 0.),
            CellWallPosition::Right => Vec2::new(0.5 - self.get_wall_thickness() / 2., 0.),
            CellWallPosition::Top => Vec2::new(0., 0.5 - self.get_wall_thickness() / 2.),
            CellWallPosition::Bottom => Vec2::new(0., -0.5 + self.get_wall_thickness() / 2.),
        }
    }
}

#[derive(Component)]
struct Wall;

#[derive(Component)]
struct Food;

#[derive(Component)]
struct Waste;

#[derive(Component)]
struct UGO {
    codons: Vec<Codon>,
}

impl UGO {
    fn new(codons: &Vec<Codon>) -> UGO {
        UGO {
            codons: codons.clone(),
        }
    }
}

// PartialEq and copy
#[derive(PartialEq, Copy, Clone)]
#[repr(u8)]
enum CodonType {
    // OG action codons
    None,
    Digest,
    Remove,
    Repair,
    MoveHand,
    Read,
    Write,
    // OG descriptive codons
    Food,
    Waste,
    Wall,
    WeakLoc,
    Inward,
    Outward,
    RGL {
        f_value_a: f32,
        f_value_b: f32,
        i_value_a: i32,
        i_value_b: i32,
    }, // Add value in here
    // Additional descriptive codons
    Energy, // If combined with Read, will read Energy levels into cell memory as hp/max_hp ratio
    LogicIf, // If combined with Write, will execute next codons depending on the last read state (if read num > 0.5)
}

impl CodonType {
    fn get_color(&self) -> Color {
        match self {
            CodonType::None => Color::rgb(0.0, 0.0, 0.0),
            CodonType::Digest => Color::rgb(0.4, 0.0, 0.8),
            CodonType::Remove => Color::rgb(0.7, 0.6, 0.04),
            CodonType::Repair => Color::rgb(0.0, 0.6, 0.0),
            CodonType::MoveHand => Color::rgb(0.8, 0.0, 0.4),
            CodonType::Read => Color::rgb(0.3, 0.3, 1.0),
            CodonType::Write => Color::rgb(0.0, 0.0, 0.85),
            CodonType::Food => Color::rgb(1.0, 0.0, 0.3),
            CodonType::Waste => Color::rgb(0.4, 0.25, 0.0),
            CodonType::Wall => Color::rgb(0.7, 0.3, 0.7),
            CodonType::WeakLoc => Color::rgb(0.3, 0.7, 0.3),
            CodonType::Inward => Color::rgb(0.0, 0.4, 0.4),
            CodonType::Outward => Color::rgb(0.0, 0.8, 0.8),
            CodonType::RGL {
                f_value_a,
                f_value_b,
                i_value_a,
                i_value_b,
            } => Color::rgb(0.5, 0.5, 0.5),
            CodonType::Energy => Color::rgb(0.8, 0.0, 0.1),
            CodonType::LogicIf => Color::rgb(1.0, 0.8, 1.0),
        }
    }
}

#[derive(PartialEq)]
enum HandDirection {
    Inward,
    Outward,
}

#[derive(PartialEq, Clone, Component)]
struct Codon {
    type_: CodonType,
    health: f32,
}

impl Codon {
    fn replace_from_other(&mut self, other: &Codon) {
        self.type_ = other.type_;
        self.health = other.health;
    }
}

#[derive(Component)]
struct CodonEntity {
    codon_index: usize,
}

impl CodonEntity {
    fn new(codon_index: usize) -> CodonEntity {
        CodonEntity { codon_index }
    }
}

#[derive(Component)]
struct Genome {
    codons: Vec<Codon>,
    codon_entities: Vec<Entity>,
    changed_codons: bool,
    update_only: bool,
}

impl Genome {
    // default genome
    fn new(genome: Vec<Codon>, child_builder: &mut Commands) -> Genome {
        let mut codon_entities = Genome::create_codon_entities(&genome, child_builder);
        Genome { codons: genome, codon_entities,changed_codons: false, update_only: false }
    }

    fn default(child_builder: &mut Commands) -> Genome {
        let mut codons = Vec::new();
        // OG first part
        let codon_types = [
            CodonType::MoveHand,
            CodonType::Outward, // - Denotes every two codons
            CodonType::Digest,
            CodonType::Food, //
            CodonType::Remove,
            CodonType::Waste, //
            CodonType::Repair,
            CodonType::Wall, //
            CodonType::Digest,
            CodonType::Food, //
            CodonType::Remove,
            CodonType::Waste, //
            CodonType::Repair,
            CodonType::Wall, //
            CodonType::MoveHand,
            CodonType::Inward, //
            CodonType::MoveHand,
            CodonType::WeakLoc, //
            CodonType::Read,
            CodonType::RGL {
                f_value_a: 0.0,
                f_value_b: 0.0,
                i_value_a: 0,
                i_value_b: 0,
            }, //
            CodonType::Write,
            CodonType::RGL {
                f_value_a: 0.0,
                f_value_b: 0.0,
                i_value_a: 0,
                i_value_b: 0,
            }, //}
               // CodonType::MoveHand, // TEST UGO
               // CodonType::Outward, // TEST UGO
               // CodonType::Write, // TEST UGO
               // CodonType::RGL, // TEST UGO
        ];

        for codon_type in codon_types.iter() {
            match codon_type {
                CodonType::RGL {
                    f_value_a,
                    f_value_b,
                    i_value_a,
                    i_value_b,
                } => {
                    codons.push(Codon {
                        type_: *codon_type,
                        health: 1.0,
                    });
                }
                _ => {
                    codons.push(Codon {
                        type_: *codon_type,
                        health: 1.0,
                    });
                }
            }
        }
        Genome::new(codons, child_builder)
    }

    fn create_codon_entities(codons: &Vec<Codon>, commands: &mut Commands) -> Vec<Entity> {
        let mut codon_entities = spawn_codons(commands, &codons);
        codon_entities
    }

    fn get_weakest_loc(&mut self) -> usize {
        let mut weakest_loc = 0;
        let mut weakest_health = 1.0;
        for (i, codon) in self.codons.iter().enumerate() {
            if codon.health < weakest_health {
                weakest_health = codon.health;
                weakest_loc = i;
            }
        }
        weakest_loc
    }
    fn damage_codon(&mut self, loc: usize, damage: f32) {
        self.codons[loc].health -= damage;
        if self.codons[loc].health <= 0.0 {
            // mutating codon
            println!("Codon at {} has mutated because of low health", loc);
            self.mutate_codon_inplace(loc, true);
        }
        self.update_only = true;
    }
    fn mutate_codon(&mut self, codon: &Codon) -> Codon {
        let mut new_codon = codon.clone();
        // TODO: Implement partial mutation (only some value or type)
        // let mut codon = codons[loc].clone();
        let codon_type = rand::random::<u32>() % (17);
        new_codon.type_ = match codon_type {
            0 => CodonType::None,
            1 => CodonType::Digest,
            2 => CodonType::Remove,
            3 => CodonType::Repair,
            4 => CodonType::MoveHand,
            5 => CodonType::Read,
            6 => CodonType::Write,
            7 => CodonType::Food,
            8 => CodonType::Waste,
            9 => CodonType::Wall,
            10 => CodonType::WeakLoc,
            11 => CodonType::Inward,
            12 => CodonType::Outward,
            13 => CodonType::RGL {
                f_value_a: 0.,
                f_value_b: 0.,
                i_value_a: 0,
                i_value_b: 0,
            },
            14 => CodonType::Energy,
            15 => CodonType::LogicIf,
            _ => CodonType::None,
        };
        new_codon.health = 1.0;
        match new_codon.type_ {
            CodonType::RGL {
                f_value_a,
                f_value_b,
                i_value_a,
                i_value_b,
            } => {
                let mut temp_i_value_a = ((rand::random::<f32>() - 0.5) * 120.0) as i32;
                let mut temp_i_value_b = ((rand::random::<f32>() - 0.5) * 120.0) as i32;
                if temp_i_value_b < temp_i_value_a {
                    let temp = temp_i_value_a;
                    temp_i_value_a = temp_i_value_b;
                    temp_i_value_b = temp;
                }
                new_codon.type_ = CodonType::RGL {
                    f_value_a: rand::random::<f32>(),
                    f_value_b: rand::random::<f32>(),
                    i_value_a: temp_i_value_a,
                    i_value_b: temp_i_value_b,
                };
            }
            _ => {}
        }
        self.update_only = true;

        new_codon
    }

    fn mutate_codon_inplace(&mut self, loc: usize, total_mutation: bool) {
        let new_codon = self.mutate_codon(&self.codons[loc].clone());
        self.codons[loc].replace_from_other(&new_codon);
        println!("Codon at {} has mutated", loc,);
        self.update_only = true;
    }

    fn read_genome(&mut self, start: i32, end: i32) -> Vec<Codon> {
        let mut new_codons = Vec::new();
        for i in start..end + 1 {
            //  random mutate chance

            let mut pos_i = i % (self.codons.len() as i32);
            if pos_i < 0 {
                pos_i = pos_i + (self.codons.len() as i32);
            }
            let pos_u: usize = pos_i as usize;
            let mut codon = self.codons[pos_u].clone();
            if rand::random::<f32>() < CODON_READ_MUTATION_RATE {
                let mutated = self.mutate_codon(&codon);
                codon.replace_from_other(&mutated);
            }
            new_codons.push(codon);
        }
        new_codons
    }
    fn write_genome(&mut self, start: i32, end: i32, new_codons: Vec<Codon>) {
        // TODO: Fix when start and end are > than the genome length
        let mut reading_hand_pos = 0;
        for i in start..end + 1 {
            // println!("Writing genome at {}", i);
            let mut pos_i = i % (self.codons.len() as i32);
            if pos_i < 0 {
                pos_i = pos_i + (self.codons.len() as i32);
            }
            let pos_u: usize = pos_i as usize;
            if reading_hand_pos >= new_codons.len() {
                // TODO: Implement a wrap around for the future
                break;
            }
            let mut new_codon = new_codons[reading_hand_pos].clone();
            new_codon.health = 1.0; // setting health to 1.0 when writing
            self.codons[pos_u].replace_from_other(&new_codon);
            reading_hand_pos += 1;
        }
        self.update_only = true;
    }
    fn insert_genes(&mut self, start: &usize, new_codons: &Vec<Codon>) {
        let mut i = *start;
        for new_codon in new_codons {
            self.codons.insert(i, new_codon.clone());
            i += 1;
        }
        self.changed_codons = true;
    }
}

#[derive(Component)]
struct CellMemory {
    float_level: f32,
    discrete_number_a: usize, // i.e. starting hand position
    disctete_number_b: usize, // i.e. end hand position
    codons: Vec<Codon>,
}

impl CellMemory {
    fn new() -> CellMemory {
        CellMemory {
            float_level: 0.0,
            discrete_number_a: 0,
            disctete_number_b: 0,
            codons: Vec::new(),
        }
    }

    fn read_float_level(&self) -> f32 {
        self.float_level
    }

    fn codons_to_memory(&mut self, codons: Vec<Codon>) {
        self.codons = codons;
    }
}

#[derive(Component)]
struct CodonExecutor {
    current_codon_reader: usize,
    last_executed_codon: usize,
    codon_execution_timer: CodonExecutionTimer,
}

impl CodonExecutor {
    fn new() -> CodonExecutor {
        let codon_execution_rate = CODON_EXECUTION_RATE
            + CODON_EXECUTION_RATE * CODON_EXECUTION_RATE_VARIATION * (rand::random::<f32>() - 0.5);
        CodonExecutor {
            current_codon_reader: 0,
            last_executed_codon: 0,
            codon_execution_timer: CodonExecutionTimer(Timer::from_seconds(
                codon_execution_rate,
                true,
            )),
        }
    }
}

#[derive(Component)]
struct CodonExecutorEntity;

#[derive(Component)]
struct CellEnergy {
    energy: f32,     // 0.0 - 100.0
    max_energy: f32, // maximum amount of energy
    waste: f32, // generated waste that needs to be cleared. If 0.0, then the cell can have max_energy, if > 0.0, then the cell can have max_energy - waste
}

impl CellEnergy {
    fn new() -> CellEnergy {
        CellEnergy {
            energy: DEFAULT_CELL_ENERGY,
            max_energy: DEFAULT_CELL_ENERGY,
            waste: 0.0,
        }
    }
    fn remove_energy(&mut self, amount: f32) -> bool {
        if self.energy >= amount {
            self.energy -= amount;
            self.waste += amount / 1.5;
            true
        } else {
            false
        }
    }
    fn add_energy(&mut self, amount: f32) -> bool {
        if self.energy + self.waste + ENERGY_RESOLUTION + amount >= self.max_energy {
            // If the cell is full of energy + waste
            false
        } else if self.energy + self.waste + amount / 3.0 < self.max_energy {
            // The cell has space for at least a third of the energy
            self.energy += amount;
            true
        } else {
            false
        }
    }
    fn remove_waste(&mut self, amount: f32) -> bool {
        if self.waste >= amount {
            self.waste -= amount;
            true
        } else {
            false
        }
    }
}

#[derive(Component)]
struct CellHand {
    hand_position: usize,
    hand_direction: HandDirection,
}

impl CellHand {
    fn new() -> CellHand {
        CellHand {
            hand_position: 0,
            hand_direction: HandDirection::Inward,
        }
    }
}

#[derive(Component)]
struct Cell {}

impl Cell {
    fn new() -> Cell {
        let codon_execution_rate = CODON_EXECUTION_RATE
            + CODON_EXECUTION_RATE * CODON_EXECUTION_RATE_VARIATION * (rand::random::<f32>() - 0.5);
        Cell {}
    }
}

#[derive(Bundle)]
struct CellBundle {
    #[bundle]
    sprite_bundle: SpriteBundle,
    cell: Cell,
    collider: Collider,
}

impl CellBundle {
    fn new(sprite_bundle: SpriteBundle, cell: Cell, collider: Collider) -> CellBundle {
        CellBundle {
            sprite_bundle,
            cell,
            collider,
        }
    }
}

#[derive(Bundle)]
struct CellWallBundle {
    #[bundle]
    cell_wall: SpriteBundle,
    collider: Collider,
}

impl CellWallBundle {
    fn new(
        cell_wall: &CellWall,
        parent: &ChildBuilder,
        cell_wall_position: &CellWallPosition,
        collider: Collider,
    ) -> CellWallBundle {
        let cell_walls = spawn_cell_wall(&parent, cell_wall, cell_wall_position);
        CellWallBundle {
            cell_wall: cell_walls,
            collider,
        }
    }
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

fn energy_system(mut commands: Commands, mut query: Query<(Entity, &mut Cell, &mut CellEnergy)>) {
    for (entity, mut cell, mut cell_energy) in query.iter_mut() {
        cell_energy.remove_energy(CELL_ENERGY_LOSS_RATE);
        if cell_energy.energy - ENERGY_RESOLUTION <= 0.0 {
            // TODO: Handle despawning in the cell
            commands.entity(entity).despawn_recursive();
        }
    }
}

// for entities with energy, we will change their size depending on their energy
fn energy_color_system(mut query: Query<(&Cell, &CellEnergy, &mut Sprite, &CellType)>) {
    let red = Color::rgb(1.0, 0.0, 0.0);
    for (cell, cell_energy, mut sprite, cell_type) in query.iter_mut() {
        let cell_color = cell_type.to_color();
        let scale = cell_energy.energy / cell_energy.max_energy;
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
                let cell_wall = CellWall::new();
                let cell = Cell::new();
                let parent_cell = commands
                    .spawn_bundle(CellBundle::new(
                        SpriteBundle {
                            sprite: Sprite {
                                color: cell_type.to_color(),
                                ..default()
                            },
                            transform: Transform {
                                translation: cell_position.extend(CELL_Z_LAYER),
                                scale: Vec3::new(CELL_SIZE.x, CELL_SIZE.y, 1.0),
                                ..default()
                            },
                            ..default()
                        },
                        cell,
                        Collider,
                    ))
                    .insert(CellType::Cell)
                    .insert(RigidBody::Fixed)
                    .insert(CellMemory::new())
                    .insert(CodonExecutor::new())
                    .insert(CellHand::new())
                    .insert(CellEnergy::new())
                    .with_children(|parent| {
                        parent
                            .spawn_bundle(CellWallBundle::new(
                                &cell_wall,
                                &parent,
                                &CellWallPosition::Left,
                                Collider,
                            )).insert(Collider)
                            // .insert(RigidBody::Fixed)
                            ;
                    })
                    // .with_chi
                    .id();
                
                let mut genome = Genome::default(&mut commands);
                // let cell_codon_entities = genome.spawn_codons(&mut commands, &mut asset_server);
                let cell_genome_executor_entity =
                    spawn_cell_genome_executor(&mut commands, &genome.codons);
                commands
                    .entity(parent_cell)
                    .push_children(&genome.codon_entities);
                commands
                    .entity(parent_cell)
                    .push_children(&[cell_genome_executor_entity]);
                commands
                    .entity(parent_cell)
                    .insert(genome);
            } else if cell_type == CellType::Wall {
                commands
                    .spawn()
                    .insert_bundle(SpriteBundle {
                        sprite: Sprite {
                            color: cell_type.to_color(),
                            ..default()
                        },
                        transform: Transform {
                            translation: cell_position.extend(CELL_Z_LAYER),
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
    let mut codons = Vec::new();
        // OG first part
    let codon_types = [
            CodonType::MoveHand, // TEST UGO
            CodonType::Inward, // TEST UGO
            CodonType::Read, // TEST UGO
            CodonType::RGL {
            f_value_a: 0.0,
            f_value_b: 0.0,
            i_value_a: -7,
            i_value_b: 7,},
            CodonType::MoveHand, // TEST UGO
            CodonType::Outward, // TEST UGO
            CodonType::Write, // TEST UGO
            CodonType::RGL {
            f_value_a: 0.0,
            f_value_b: 0.0,
            i_value_a: 0,
            i_value_b: 0,
        }, // TEST UGO
    ];

    for codon_type in codon_types.iter() {
        match codon_type {
            CodonType::RGL {
                f_value_a,
                f_value_b,
                i_value_a,
                i_value_b,
            } => {
                codons.push(Codon {
                    type_: *codon_type,
                    health: 1.0,
                });
            }
            _ => {
                codons.push(Codon {
                    type_: *codon_type,
                    health: 1.0,
                });
            }
        }
    }


    spawn_ugo(&mut commands, &Transform {
        translation: Vec3::new(-15.0, -15.0, 1.0),

        scale: Vec3::new(CELL_SIZE.x, CELL_SIZE.y, 1.0),
        ..default()
    }, &codons );
}

fn get_primary_window_size(windows: &Res<Windows>) -> Vec2 {
    let window = windows.get_primary().unwrap();
    let window = Vec2::new(window.width() as f32, window.height() as f32);
    window
}

/// Pan the camera with middle mouse click, zoom with scroll wheel, orbit with right mouse click.
fn pan_orbit_camera(
    windows: Res<Windows>,
    mut ev_motion: EventReader<MouseMotion>,
    mut ev_scroll: EventReader<MouseWheel>,
    input_mouse: Res<Input<MouseButton>>,
    mut query: Query<(&mut Camera, &mut Transform, &mut OrthographicProjection)>,
) {
    // change input mapping for orbit and panning here
    // let orbit_button = MouseButton::Right;
    let pan_button = MouseButton::Right;

    let mut pan = Vec2::ZERO;
    let mut scroll = 0.0;

    if input_mouse.pressed(pan_button) {
        // Pan only if we're not rotating at the moment
        for ev in ev_motion.iter() {
            pan += ev.delta;
        }
    }
    for ev in ev_scroll.iter() {
        scroll += ev.y;
    }

    for (mut pan_orbit, mut transform, mut projection) in query.iter_mut() {
        let mut any = false;
        let mut pan_scale = Vec2::new(1.0, 1.0);
        if pan.length_squared() > 0.0 {
            any = true;
            println!("pan: {:?}", pan);
            // make panning distance independent of resolution,
            // let window = get_primary_window_size(&windows);
            // let window_aspect_ratio = window.x / window.y;
            // pan_scale = Vec2::new(
            //     window_aspect_ratio * 1.0,
            //     1.0,
            // );
        } else if scroll.abs() > 0.0 {
            any = true;
        }

        if any {
            transform.translation.x -= pan.x * pan_scale.x * projection.scale; // x is inverted
            transform.translation.y += pan.y * pan_scale.y * projection.scale;

            projection.scale = projection.scale - scroll * 0.1;
            if projection.scale < 0.1 {
                projection.scale = 0.1;
            }
        }
    }
}

fn spawn_codons(child_builder: &mut Commands, codons: &Vec<Codon>) -> Vec<Entity> {
    // Spawning each codon in circle inside the parent
    let mut codon_angle: f32 = 0.0;
    let codon_width: f32 = CODON_SIZE * 16.0 / codons.len() as f32;
    let mut codon_angle_step: f32 = 2.0 * PI / codons.len() as f32;
    let mut codon_entities = vec![];
    for (idx, codon) in codons.iter().enumerate() {
        let codon_position = Vec2::new(
            CODON_RADIUS * codon_angle.cos(),
            CODON_RADIUS * codon_angle.sin(),
        );
        let codon_position = codon_position.extend(CODON_Z_LAYER);
        let codon_color = codon.type_.get_color();
        let codon_entity = child_builder
            .spawn()
            .insert_bundle(SpriteBundle {
                sprite: Sprite {
                    color: codon_color,
                    ..default()
                },
                transform: Transform {
                    translation: codon_position,
                    scale: Vec3::new(CODON_SIZE, codon_width, 1.0),
                    rotation: Quat::from_rotation_z(codon_angle),
                    ..default()
                },
                ..default()
            })
            .insert(CodonEntity { codon_index: idx })
            .id();
        codon_angle += codon_angle_step;
        codon_entities.push(codon_entity);
    }
    codon_entities
}

fn spawn_cell_genome_executor(commands: &mut Commands, codons: &Vec<Codon>) -> Entity {
    let codon_angle: f32 = 0.0;
    let codon_width: f32 = CODON_SIZE * 16.0 / codons.len() as f32;
    let reader_position = Vec2::new(
        CODON_RADIUS * codon_angle.cos(),
        CODON_RADIUS * codon_angle.sin(),
    );
    let reader_position = reader_position.extend(CODON_READER_Z_LAYER);
    let codon_entity = commands
        .spawn()
        .insert_bundle(SpriteBundle {
            sprite: Sprite {
                color: GENOME_READER_COLOR,
                ..default()
            },
            transform: Transform {
                translation: reader_position,
                scale: Vec3::new(CODON_SIZE * 1.25, codon_width * 1.25, 1.0),
                rotation: Quat::from_rotation_z(codon_angle),
                ..default()
            },
            ..default()
        })
        .insert(CodonExecutorEntity)
        .id();
    codon_entity
}

fn spawn_cell_wall(
    commands: &ChildBuilder,
    cell_wall: &CellWall,
    cell_wall_position: &CellWallPosition,
) -> SpriteBundle {
    // left wall
    let wall_position = cell_wall.get_wall_position_part(cell_wall_position);
    let wall_thickness = cell_wall.get_wall_thickness();
    let scale = match cell_wall_position {
        CellWallPosition::Left => Vec3::new(wall_thickness, 1.0, 1.0),
        CellWallPosition::Right => Vec3::new(wall_thickness, 1.0, 1.0),
        CellWallPosition::Top => Vec3::new(1.0, wall_thickness, 1.0),
        CellWallPosition::Bottom => Vec3::new(1.0, wall_thickness, 1.0),
    };
    SpriteBundle {
        sprite: Sprite {
            color: CELL_WALL_COLOR,
            ..default()
        },
        transform: Transform {
            translation: wall_position.extend(CELL_WALL_Z_LAYER),
            scale: scale,
            ..default()
        },
        ..default()
    }
}

fn spawn_ugo(commands: &mut Commands, cell_pos: &Transform, codons: &Vec<Codon>) -> Entity {
    let ugo_position = {
        let rand_sign = if rand::random() { 1. } else { -1. };
        let rand_wall = rand::random::<bool>();
        let mut ugo_position = Vec2::new(
            cell_pos.translation.x + rand_sign * CELL_SIZE.x / 1.75, // A bit further away from the cell
            cell_pos.translation.y + rand_sign * CELL_SIZE.y / 1.75,
        );
        if rand_wall {
            ugo_position.x += rand::random::<f32>() * CELL_SIZE.x;
        } else {
            ugo_position.y += rand::random::<f32>() * CELL_SIZE.y;
        }
        ugo_position
    };

    let ugo_direction = {
        let left_direction: bool = ugo_position.x < cell_pos.translation.x;
        let up_direction: bool = ugo_position.y < cell_pos.translation.y;
        let direction = Vec2::new(
            // Spawning away from self
            if left_direction { -1. } else { 1. },
            if up_direction { -1. } else { 1. },
        );

        let mut random_speed =
            Vec2::new(random::<f32>() * 25.0 + 25.0, random::<f32>() * 25.0 + 25.0);
        random_speed.x *= direction.x;
        random_speed.y *= direction.y;

        let random_direction_degrees = (random::<f32>() - 0.5) * 90.; // Spawning UGOs 90 degrees away from self
        let random_direction = Vec2::new(
            random_speed.x * random_direction_degrees.cos(),
            random_speed.y * random_direction_degrees.sin(),
        );

        random_direction
    };
    let shape = shapes::Circle {
        radius: 1.0 / 2.0,
        center: Vec2::new(0., 0.),
    };
    let ugo_entity = commands
        .spawn_bundle(GeometryBuilder::build_as(
            &shape,
            DrawMode::Outlined {
                fill_mode: FillMode::color(UGO_COLOR),
                outline_mode: StrokeMode::new(Color::WHITE, 1. / 32.),
            },
            Transform {
                translation: ugo_position.extend(UGO_Z_LAYER),
                scale: Vec3::new(UGO_SIZE, UGO_SIZE, 1.0),
                ..default()
            },
        ))
        .insert(Collider)
        .insert(Velocity(ugo_direction))
        .insert(Particle)
        .insert(UGO::new(codons))
        .id();
    let ugo_codons = spawn_codons(commands, codons);
    commands.entity(ugo_entity).push_children(&ugo_codons);
    ugo_entity
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
    mut particle_query: Query<
        (
            Entity,
            &mut Velocity,
            &Transform,
            Option<&Food>,
            Option<&UGO>,
        ),
        With<Particle>,
    >,
    mut collider_query: Query<
        (
            &Transform,
            // Entity,
            Option<(&mut Cell, &mut CellEnergy)>,
            Option<(&mut Genome, &mut CodonExecutor)>,
        ),
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
    for (particle_entity, mut particle_velocity, particle_transform, maybe_food, maybe_ugo) in
        particle_query.iter_mut()
    {
        // check collision with walls
        for (transform, mut maybe_cell_and_energy, mut maybe_genome_and_executor) in
            collider_query.iter_mut()
        {
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
                    if maybe_cell_and_energy.is_some() {
                        let (cell, energy) = maybe_cell_and_energy.as_mut().unwrap();
                        if energy.add_energy(FOOD_ENERGY) {
                            commands.entity(particle_entity).despawn();
                            break;
                        }
                    }
                }
                if maybe_ugo.is_some() {
                    if maybe_genome_and_executor.is_some() {
                        // insert UGO genes to cell
                        let (genome, executor) = maybe_genome_and_executor.as_mut().unwrap();
                        let ugo = maybe_ugo.unwrap();
                        if genome.codons.len() + ugo.codons.len() > MAX_CELL_CODONS {
                            // UGO has too many genes
                        } else {
                            genome.insert_genes(&executor.current_codon_reader, &ugo.codons);
                            // despawn old cell codons
                            // commands.entity(cell_entity).despawn_descendants();
                            // // Spawning new genes for the UGO
                            // spawn_codons(commands, &genome.codons);
                            // spawn_cell_genome_executor(commands, &genome.codons);
                            commands.entity(particle_entity).despawn_recursive();
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
                // if collide_y || collide_x {
                //     break;
                // }
            }
        }
    }
}

fn food_dispenser(
    mut commands: Commands,
    time: Res<Time>,
    mut timer: ResMut<FoodSpawnTimer>,
    food_query: Query<(&Particle, &Food)>,
    cells_query: Query<&Transform, (With<Collider>, Without<Food>)>,
) {
    if !timer.0.tick(time.delta()).just_finished() {
        // If timer did not finish, do nothing
        return;
    }
    let mut spawned_num: usize = 0;
    // Check if we have enough food
    let center_vec = Vec2::new(0., 0.);
    let food_size_vec = Vec2::new(FOOD_SIZE, FOOD_SIZE);
    while food_query.iter().count() < MAX_FOOD && spawned_num < FOOD_SPAWN_RATE {
        let mut random_food_position: Vec2 = Vec2::new(0., 0.);
        loop {
            let mut found_empty_spot = true;
            random_food_position = Vec2::new(
                random::<f32>() * WORLD_SIZE * WORLD_SCALE - WORLD_SIZE * WORLD_SCALE / 2.0,
                random::<f32>() * WORLD_SIZE * WORLD_SCALE - WORLD_SIZE * WORLD_SCALE / 2.0,
            );
            for cell_transform in cells_query.iter() {
                let collision = collide(
                    random_food_position.extend(FOOD_Z_LAYER),
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
        let random_speed = Vec2::new(random::<f32>() * 25.0 + 25.0, random::<f32>() * 25.0 + 25.0);
        let random_direction_degrees = random::<f32>() * 360.0;
        let random_direction = Vec2::new(
            random_speed.x * random_direction_degrees.cos(),
            random_speed.y * random_direction_degrees.sin(),
        );
        commands
            .spawn_bundle(GeometryBuilder::build_as(
                &shape,
                DrawMode::Outlined {
                    fill_mode: FillMode::color(FOOD_COLOR),
                    outline_mode: StrokeMode::new(Color::BLACK, 1.),
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
        spawned_num += 1;
    }
}

fn codon_executing(
    mut commands: Commands,
    time: Res<Time>,
    mut cell_query: Query<(
        &mut Cell,
        &mut Genome,
        &mut CodonExecutor,
        &mut CellHand,
        &mut CellEnergy,
        &mut CellMemory,
        &mut Children,
        &Transform,
    )>,
    codon_query: Query<&CodonEntity>,
    pool: Res<AsyncComputeTaskPool>,
) {
    for (
        mut cell,
        mut cell_genome,
        mut cell_codon_executor,
        mut cell_hand,
        mut cell_energy,
        mut cell_memory,
        mut children,
        cell_pos,
    ) in cell_query.iter_mut()
    {
        if !cell_codon_executor
            .codon_execution_timer
            .0
            .tick(time.delta())
            .just_finished()
        {
            continue; // If the timer hasn't ticked, skip this cell
        }
        let mut errors = vec![];

        let cell_codons: Vec<&CodonEntity> = children
            .iter()
            .map(|child| codon_query.get(*child))
            .filter_map(|r| r.map_err(|e| errors.push(e)).ok())
            .collect();
        // gets current codon from hand position
        let current_codon_reader = cell_codon_executor.current_codon_reader;
        let cur_codon = cell_genome.codons.get(current_codon_reader).unwrap();
        let prev_codon = cell_genome
            .codons
            .get(cell_codon_executor.last_executed_codon)
            .unwrap();
        // does action depending on current codon
        // println!("Codon executing");
        loop {
            // Looping so we can break if we don't want to nest ifs too deep.
            match cur_codon.type_ {
                CodonType::None => {
                    // println!("Doing nothing");
                    // do nothing
                }
                CodonType::Digest => {
                    // println!("Digesting");
                    // cell.digest();
                }
                CodonType::Remove => {
                    // println!("Removing");
                    // cell.remove();
                }
                CodonType::Repair => {
                    // println!("Repairing");
                    // cell.repair();
                }
                CodonType::MoveHand => {
                    // println!("Moving hand");
                    // cell.move_hand();
                }
                CodonType::Read => {
                    // println!("Reading");
                    // cell.read();
                }
                CodonType::Write => {
                    // println!("Writing");
                    // cell.write();
                }
                CodonType::Food => {
                    // println!("Food");
                    match prev_codon.type_ {
                        CodonType::Digest => {
                            // cell.eat_food();
                        }
                        _ => {
                            // cell.eat_food();
                        }
                    }
                }
                CodonType::Waste => {
                    // println!("Waste");
                    if cell_hand.hand_direction == HandDirection::Inward {
                        // cannot remove waste if hand is in inward direction
                        println!("Cannot remove waste if hand is in inward direction");
                        break;
                    }
                    match prev_codon.type_ {
                        CodonType::Remove => {
                            let removed_waste = cell_energy.remove_waste(FOOD_ENERGY * 0.6); // Frees up some energy potential from waste
                            if removed_waste {
                                // println!("Removed waste");
                                spawn_waste(&mut commands, cell_pos);
                            }
                        }
                        _ => {
                            // println!("Not removed waste??");
                        }
                    }
                }
                CodonType::Wall => {
                    // println!("Wall");
                    // cell.wall();
                }
                CodonType::WeakLoc => {
                    // println!("Weak loc");
                    // cell.weak_loc();
                    if cell_hand.hand_direction == HandDirection::Outward {
                        // cannot weak loc if hand is in Outward direction
                        // println!("Cannot weak loc if hand is in Outward direction");
                        break;
                    }
                    match prev_codon.type_ {
                        CodonType::MoveHand => {
                            let weakest_loc = cell_genome.get_weakest_loc();
                            cell_hand.hand_position = weakest_loc;
                        }
                        _ => {
                            // skip
                        }
                    }
                }
                CodonType::Inward => {
                    // println!("Inward");
                    // cell.inward();
                    match prev_codon.type_ {
                        CodonType::MoveHand => {
                            cell_hand.hand_direction = HandDirection::Inward;
                        }
                        _ => {
                            // cell.repair();
                        }
                    }
                }
                CodonType::Outward => {
                    // println!("Outward");
                    // cell.outward();
                    match prev_codon.type_ {
                        CodonType::MoveHand => {
                            cell_hand.hand_direction = HandDirection::Outward;
                        }
                        _ => {
                            // cell.repair();
                        }
                    }
                }
                CodonType::RGL {
                    f_value_a,
                    f_value_b,
                    i_value_a,
                    i_value_b,
                } => {
                    if cell_hand.hand_direction == HandDirection::Outward {
                        match prev_codon.type_ {
                            CodonType::Write => {
                                spawn_ugo(&mut commands, cell_pos, &cell_memory.codons);
                            }
                            _ => {
                                // cell.repair();
                            }
                        }
                        break;
                    }
                    match prev_codon.type_ {
                        CodonType::Read => {
                            let read_start = i_value_a + cell_hand.hand_position as i32;
                            let read_end = i_value_b + cell_hand.hand_position as i32;
                            let read_codons = cell_genome.read_genome(read_start, read_end);
                            cell_memory.codons = read_codons;
                        }
                        CodonType::Write => {
                            // println!("Writing genome!");
                            let write_start = i_value_a + cell_hand.hand_position as i32;
                            let write_end = i_value_b + cell_hand.hand_position as i32;
                            let write_codons = cell_memory.codons.clone();
                            cell_genome.write_genome(write_start, write_end, write_codons);
                        }
                        _ => {
                            // cell.repair();
                        }
                    }
                    // println!("RGL");
                    // cell.rgl();
                }
                CodonType::Energy => {
                    // println!("Energy");
                    // cell.energy();
                }
                CodonType::LogicIf => {
                    // println!("Logic if");
                    // cell.logicif();
                }
            }
            break;
        }
        // remove health for codon
        // cur_codon.health -= CODON_HEALTH_COST_PER_EXECUTION;
        cell_genome.damage_codon(current_codon_reader, CODON_HEALTH_COST_PER_EXECUTION);

        cell_codon_executor.last_executed_codon = cell_codon_executor.current_codon_reader;
        cell_codon_executor.current_codon_reader += 1;
        if cell_codon_executor.current_codon_reader >= cell_genome.codons.len() {
            cell_codon_executor.current_codon_reader = 0;
        }
        cell_energy.remove_energy(CODON_ENERGY_COST_PER_EXECUTION);
    }
}

fn spawn_waste(commands: &mut Commands, cell_pos: &Transform) {
    let center_vec = Vec2::new(0., 0.);

    let mut random_waste_position = Vec2::new(
        random::<f32>() * CELL_SIZE.x - CELL_SIZE.x / 2.0,
        random::<f32>() * CELL_SIZE.y - CELL_SIZE.y / 2.0,
    );
    random_waste_position.x += cell_pos.translation.x;
    random_waste_position.y += cell_pos.translation.y;
    let shape = shapes::Circle {
        radius: FOOD_SIZE / 2.0,
        center: center_vec,
    };
    let random_direction = Vec2::new(random::<f32>() * 50.0 - 25.0, random::<f32>() * 50.0 - 25.0);

    commands
        .spawn_bundle(GeometryBuilder::build_as(
            &shape,
            DrawMode::Outlined {
                fill_mode: FillMode::color(WASTE_COLOR),
                outline_mode: StrokeMode::new(Color::BLACK, 1.),
            },
            Transform {
                translation: random_waste_position.extend(WASTE_Z_LAYER),
                scale: Vec3::new(1.0, 1.0, 1.0),
                ..default()
            },
        ))
        .insert(Waste)
        .insert(Collider)
        .insert(Velocity(random_direction))
        .insert(Particle);
}

fn waste_despawner(
    mut commands: Commands,
    time: Res<Time>,
    mut timer: ResMut<WasteDespawnTimer>,
    mut query: Query<Entity, With<Waste>>,
) {
    if !timer.0.tick(time.delta()).just_finished() {
        // If timer did not finish, do nothing
        return;
    }
    let mut despawned = 0;
    let total_waste = query.iter().count();
    for entity in query.iter() {
        commands.entity(entity).despawn();
        despawned += 1;
        if despawned > WASTE_DESPAWN_PER_TIMER && total_waste - despawned < MAX_WASTE {
            break;
        }
    }
}

fn update_cell_genome_executor(
    mut query_codon_executor: Query<(&Parent, &mut Transform), With<CodonExecutorEntity>>,
    query_cells: Query<(&CodonExecutor, &Genome)>,
) {
    for (parent, mut executor_transform) in query_codon_executor.iter_mut() {
        let result = query_cells.get(parent.0);
        let (codon_executor, genome) = result.unwrap();
        let codon_reader_pos = codon_executor.current_codon_reader;
        let codon_width: f32 = CODON_SIZE * 16.0 / genome.codons.len() as f32;
        let codon_angle_step: f32 = 2.0 * PI / genome.codons.len() as f32;
        let codon_angle: f32 = codon_angle_step * codon_reader_pos as f32;
        let reader_position = Vec2::new(
            CODON_RADIUS * codon_angle.cos(),
            CODON_RADIUS * codon_angle.sin(),
        );
        let reader_position = reader_position.extend(CODON_READER_Z_LAYER);
        executor_transform.translation = reader_position;
        executor_transform.scale = Vec3::new(CODON_SIZE * 1.25, codon_width * 1.25, 1.0);
        executor_transform.rotation = Quat::from_rotation_z(codon_angle);
    }
}
fn change_cell_genome(
    mut commands: Commands,
    mut query_cells: Query<(&mut Genome, Entity)>,
) {
    for (mut genome, cell_entity) in query_cells.iter_mut() {
        // if parent not found, skip
        if !genome.changed_codons && !genome.update_only {
            continue;
        }
        // popping codon entities and despawning them
        if genome.changed_codons { 
            while genome.codon_entities.len() > 0 {
                let codon_entity = genome.codon_entities.pop().unwrap();
                commands.entity(codon_entity).despawn();
            }
            // spawning codon entities
            genome.codon_entities = spawn_codons(&mut commands, &genome.codons);

            commands.entity(cell_entity).push_children(&genome.codon_entities);
            genome.changed_codons = false;
        }
    }
}

fn update_cell_genome(
    mut query_codons: Query<
        (&Parent, &mut Transform, &mut Sprite, &CodonEntity),
        With<CodonEntity>,
    >,
    query_cells: Query<&Genome>,
) {
    for (parent, mut codon_transform, mut sprite, codon_entity) in query_codons.iter_mut() {
        // if parent not found, skip
        let genome = {
            match query_cells.get(parent.0) {
                Ok(genome) => genome,
                Err(e) => continue,
            }
        };
        let codon_index = codon_entity.codon_index;
        let codon = genome.codons.get(codon_index).unwrap();
        let codon_width = CODON_SIZE * 16.0 / genome.codons.len() as f32 * codon.health;
        codon_transform.scale = Vec3::new(CODON_SIZE, codon_width, 1.0);
        sprite.color = codon.type_.get_color();
    }
}
