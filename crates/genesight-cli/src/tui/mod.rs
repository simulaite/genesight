mod app;
mod event;
mod ui;

use std::io;
use std::sync::mpsc;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};
use crossterm::{
    event::{DisableMouseCapture, EnableMouseCapture},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use genesight_core::models::Report;
use ratatui::{backend::CrosstermBackend, Terminal};

use app::App;

pub use app::AnalysisProgress;

/// Run the interactive TUI dashboard for the given report (legacy path).
///
/// This takes ownership of the terminal (raw mode + alternate screen),
/// runs the event loop, and restores the terminal on exit -- including
/// on panics.
#[allow(dead_code)]
pub fn run(report: Report) -> Result<()> {
    // Install a panic hook that restores the terminal before printing the panic.
    let original_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = restore_terminal();
        original_hook(info);
    }));

    enable_raw_mode().context("failed to enable raw mode")?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)
        .context("failed to enter alternate screen")?;

    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend).context("failed to create terminal")?;

    let mut app = App::new(report);

    let result = run_dashboard_loop(&mut terminal, &mut app);

    // Always restore terminal, even on error
    restore_terminal()?;

    result
}

/// Run the interactive TUI with a loading screen, receiving analysis progress
/// from a background thread via an `mpsc` channel.
///
/// The TUI starts immediately in loading phase, displaying progress messages
/// as they arrive. Once `AnalysisProgress::Complete` is received, the app
/// transitions to the interactive dashboard.
pub fn run_with_analysis(rx: mpsc::Receiver<AnalysisProgress>) -> Result<()> {
    // Install a panic hook that restores the terminal before printing the panic.
    let original_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = restore_terminal();
        original_hook(info);
    }));

    enable_raw_mode().context("failed to enable raw mode")?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)
        .context("failed to enter alternate screen")?;

    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend).context("failed to create terminal")?;

    let mut app = App::new_loading();

    let result = run_analysis_loop(&mut terminal, &mut app, &rx);

    // Always restore terminal, even on error
    restore_terminal()?;

    result
}

/// Event loop for the analysis-loading mode. Handles both loading and dashboard phases.
fn run_analysis_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    app: &mut App,
    rx: &mpsc::Receiver<AnalysisProgress>,
) -> Result<()> {
    let mut last_spinner_tick = Instant::now();
    let spinner_interval = Duration::from_millis(80);

    loop {
        // 1. Check for progress messages from the analysis thread (non-blocking)
        loop {
            match rx.try_recv() {
                Ok(progress) => app.handle_progress(progress),
                Err(mpsc::TryRecvError::Empty) => break,
                Err(mpsc::TryRecvError::Disconnected) => {
                    // The sender was dropped. If we have no report and no error,
                    // that means the analysis thread panicked or exited unexpectedly.
                    if app.report.is_none() && app.error_message.is_none() {
                        app.error_message =
                            Some("Analysis thread terminated unexpectedly".to_string());
                    }
                    break;
                }
            }
        }

        // 2. Tick the spinner if enough time has passed
        if matches!(app.phase, app::AppPhase::Loading { .. })
            && last_spinner_tick.elapsed() >= spinner_interval
        {
            app.tick_spinner();
            last_spinner_tick = Instant::now();
        }

        // 3. Draw the UI
        terminal.draw(|frame| ui::draw(frame, app))?;

        // 4. Handle terminal events
        if event::handle_event(app)? {
            return Ok(());
        }

        if app.should_quit {
            return Ok(());
        }
    }
}

/// Main event loop for the dashboard (used by the legacy `run` path).
#[allow(dead_code)]
fn run_dashboard_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    app: &mut App,
) -> Result<()> {
    loop {
        terminal.draw(|frame| ui::draw(frame, app))?;

        if event::handle_event(app)? {
            return Ok(());
        }
    }
}

/// Restore terminal to normal mode.
fn restore_terminal() -> Result<()> {
    disable_raw_mode().context("failed to disable raw mode")?;
    execute!(io::stdout(), LeaveAlternateScreen, DisableMouseCapture)
        .context("failed to leave alternate screen")?;
    Ok(())
}
