use bevy_app::prelude::*;
use bevy_app::App;
use bevy_ecs::prelude::*;
use bevy_tasks::TaskPool;

fn main() {
    println!("Hello, world!");
}

pub fn analyze() {
    App::new()
        .add_plugin(bevy_app::ScheduleRunnerPlugin)
        .run();
}