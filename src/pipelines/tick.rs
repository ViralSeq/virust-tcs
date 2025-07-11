use indicatif::{ProgressBar, ProgressStyle};
use std::io::{self, Read};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use std::thread;
use std::time::Duration;

use crate::helper::tcs_helper::*;

pub fn run_tick() -> () {
    let running = Arc::new(AtomicBool::new(true));

    // Clone Arc for input thread
    let running_input = running.clone();

    // Spawn a thread to listen for any input
    thread::spawn(move || {
        // Wait for any input (blocks until user presses Enter or types anything)
        let _ = io::stdin().read(&mut [0u8]).unwrap();
        running_input.store(false, Ordering::SeqCst);
    });

    // Set up spinner
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template("{spinner} Cat and Mouse: {msg}")
            .unwrap()
            .tick_strings(&CLI_ANIMATION_TICK_STRINGS),
    );

    let mut step = 0;
    // Loop until the input thread sets running to false
    while running.load(Ordering::SeqCst) {
        spinner.set_message(format!("Chase step {}", step));
        thread::sleep(Duration::from_millis(120));
        spinner.inc(1);
        step += 1;
    }

    spinner.finish_with_message("The chase is over! 🐱🐭");
}
