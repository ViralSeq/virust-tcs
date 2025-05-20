use clap::builder::styling::{AnsiColor, Color};
use clap::builder::styling::{Style, Styles};
use clap::{ColorChoice, Parser, Subcommand};

pub const BANNER: &str = "\x1b[0;91m████████  ██████ ███████     ██████  ██ ██████  ███████ ██      ██ ███    ██ ███████\x1b[0m\n\
                      \x1b[0;93m   ██    ██      ██          ██   ██ ██ ██   ██ ██      ██      ██ ████   ██ ██\x1b[0m\n\
                      \x1b[0;92m   ██    ██      ███████     ██████  ██ ██████  █████   ██      ██ ██ ██  ██ █████\x1b[0m\n\
                      \x1b[0;96m   ██    ██           ██     ██      ██ ██      ██      ██      ██ ██  ██ ██ ██\x1b[0m\n\
                      \x1b[0;95m   ██     ██████ ███████     ██      ██ ██      ███████ ███████ ██ ██   ████ ███████\x1b[0m\n";

#[derive(Parser, Debug, Clone)]
#[command(
    name = "TCS pipeline",
    version = env!("CARGO_PKG_VERSION"),
    about = BANNER,
    color = ColorChoice::Always,
    styles = get_styles(),
)]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug, Clone)]
pub enum Commands {
    /// Run the TCS pipeline
    #[command(alias = "r")]
    Run {
        /// Input directory path
        #[arg(short, long)]
        input: String,

        /// param file path
        #[arg(short, long)]
        param: String,

        /// keep original files
        #[arg(long, default_value_t = false)]
        keep_original: bool,
    },

    /// Generate a param file through CLI
    #[command(alias = "g")]
    Generate {},

    /// Run the TCS HIV-1 DR Pipeline,
    DR {
        /// Input directory path
        #[arg(short, long)]
        input: String,

        /// DR version number
        #[arg(short, long, default_value_t = String::from("v1"))]
        version: String,

        /// keep original files
        #[arg(long, default_value_t = false)]
        keep_original: bool,
    },

    /// List param for the DR pipeline, w/o aurguments it will list all available version numbers.
    DrParams {
        /// Print out params for a specific version
        #[arg(short, long)]
        version: Option<String>,
    },

    /// SDRM pipeline followed by HIV-1 DR pipeline
    SDRM {
        /// Input directory path
        #[arg(short, long)]
        input: String,

        /// DR version number
        #[arg(short, long, default_value_t = String::from("v1"))]
        version: String,
    },

    /// Aggregate log files and reorganize the directory structure after TCS or DR pipeline
    Log {
        /// Input directory path
        #[arg(short, long)]
        input: String,
    },
}

pub fn get_styles() -> Styles {
    Styles::styled()
        .usage(
            Style::new()
                .bold()
                .underline()
                .fg_color(Some(Color::Ansi(AnsiColor::Yellow))),
        )
        .header(
            Style::new()
                .bold()
                .underline()
                .fg_color(Some(Color::Ansi(AnsiColor::Yellow))),
        )
        .literal(Style::new().fg_color(Some(Color::Ansi(AnsiColor::Green))))
        .invalid(
            Style::new()
                .bold()
                .fg_color(Some(Color::Ansi(AnsiColor::Red))),
        )
        .error(
            Style::new()
                .bold()
                .fg_color(Some(Color::Ansi(AnsiColor::Red))),
        )
        .valid(
            Style::new()
                .bold()
                .underline()
                .fg_color(Some(Color::Ansi(AnsiColor::Green))),
        )
        .placeholder(Style::new().fg_color(Some(Color::Ansi(AnsiColor::White))))
}
