
use clap::builder::styling::{AnsiColor, Color};
use clap::builder::styling::{Style, Styles};
use clap::{ColorChoice, Parser};

#[derive(Parser, Debug)]
#[command(
    name = "TCS Pipeline",
    version = env!("CARGO_PKG_VERSION"),
    about = "\x1b[1;91mCLI For TCS Pipeline (Rust Version)\x1b[0m",
    color = ColorChoice::Always,
    styles = get_styles(),
)]
pub struct Args {
    /// Query sequence
    #[arg(short, long, use_value_delimiter = true, value_delimiter = ' ', num_args = 1..)]
    pub query: Vec<String>,

    /// Reference genome, either HXB2 or SIVmm239
    #[arg(short, long, default_value = "HXB2")]
    pub reference: String,

    /// Type of query, either nt or aa
    #[arg(short, long, default_value = "nt")]
    pub type_query: String,

    /// algorithm for locator, 1 is accurate but slower, 2 is fast but less accurate, suitable for smaller query sequences
    #[arg(short, long, default_value_t = 1)]
    pub algorithm: u8,
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
