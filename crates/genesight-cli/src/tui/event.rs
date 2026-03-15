use anyhow::Result;
use crossterm::event::{self, Event, KeyCode, KeyEventKind, KeyModifiers};

use super::app::{App, AppPhase, Panel, ViewMode};
use genesight_core::models::ConfidenceTier;

/// Handle a single terminal event. Returns `true` if the application should quit.
///
/// During the loading phase, only quit keys (q, Ctrl-C) are handled.
/// During the dashboard phase, the full keybinding set is active.
pub fn handle_event(app: &mut App) -> Result<bool> {
    if !event::poll(std::time::Duration::from_millis(50))? {
        return Ok(false);
    }

    let ev = event::read()?;

    // Only handle key press events (not release or repeat)
    let Event::Key(key) = ev else {
        return Ok(false);
    };
    if key.kind != KeyEventKind::Press {
        return Ok(false);
    }

    match &app.phase {
        AppPhase::Loading { .. } => handle_loading_input(app, key.code, key.modifiers),
        AppPhase::Dashboard {
            search_mode,
            view_mode,
            ..
        } => {
            if *search_mode {
                handle_search_input(app, key.code, key.modifiers)
            } else if *view_mode == ViewMode::Summary {
                handle_summary_input(app, key.code, key.modifiers)
            } else {
                handle_dashboard_input(app, key.code, key.modifiers)
            }
        }
    }
}

/// Handle keyboard input during the loading phase. Only quit is supported.
fn handle_loading_input(app: &mut App, code: KeyCode, modifiers: KeyModifiers) -> Result<bool> {
    match code {
        KeyCode::Char('q') => {
            app.should_quit = true;
            Ok(true)
        }
        KeyCode::Char('c') if modifiers.contains(KeyModifiers::CONTROL) => {
            app.should_quit = true;
            Ok(true)
        }
        _ => Ok(false),
    }
}

/// Handle keyboard input during the summary view.
fn handle_summary_input(app: &mut App, code: KeyCode, modifiers: KeyModifiers) -> Result<bool> {
    match code {
        KeyCode::Char('q') => return Ok(true),
        KeyCode::Char('c') if modifiers.contains(KeyModifiers::CONTROL) => return Ok(true),

        // Switch to table view
        KeyCode::Tab => app.toggle_view(),

        // Scroll summary
        KeyCode::Char('j') | KeyCode::Down => app.scroll_summary_down(),
        KeyCode::Char('k') | KeyCode::Up => app.scroll_summary_up(),
        KeyCode::Char('g') | KeyCode::Home => {
            if let AppPhase::Dashboard { summary_scroll, .. } = &mut app.phase {
                *summary_scroll = 0;
            }
        }

        // Tier filters jump to table with filter
        KeyCode::Char('1') => {
            app.toggle_tier(ConfidenceTier::Tier1Reliable);
            app.toggle_view();
        }
        KeyCode::Char('2') => {
            app.toggle_tier(ConfidenceTier::Tier2Probable);
            app.toggle_view();
        }
        KeyCode::Char('3') => {
            app.toggle_tier(ConfidenceTier::Tier3Speculative);
            app.toggle_view();
        }

        // Search jumps to table with search
        KeyCode::Char('/') => {
            app.toggle_view();
            if let AppPhase::Dashboard {
                search_mode,
                search_query,
                ..
            } = &mut app.phase
            {
                *search_mode = true;
                search_query.clear();
            }
        }

        _ => {}
    }

    Ok(false)
}

/// Handle keyboard input during the table view (normal mode).
fn handle_dashboard_input(app: &mut App, code: KeyCode, modifiers: KeyModifiers) -> Result<bool> {
    match code {
        KeyCode::Char('q') => return Ok(true),
        KeyCode::Char('c') if modifiers.contains(KeyModifiers::CONTROL) => return Ok(true),

        // Navigation
        KeyCode::Char('j') | KeyCode::Down => {
            if let AppPhase::Dashboard { active_panel, .. } = &app.phase {
                if *active_panel == Panel::Results {
                    app.next();
                } else {
                    app.scroll_detail_down();
                }
            }
        }
        KeyCode::Char('k') | KeyCode::Up => {
            if let AppPhase::Dashboard { active_panel, .. } = &app.phase {
                if *active_panel == Panel::Results {
                    app.previous();
                } else {
                    app.scroll_detail_up();
                }
            }
        }
        KeyCode::Char('g') | KeyCode::Home => {
            if let AppPhase::Dashboard {
                active_panel,
                detail_scroll,
                ..
            } = &mut app.phase
            {
                if *active_panel == Panel::Results {
                    app.go_to_top();
                } else {
                    *detail_scroll = 0;
                }
            }
        }
        KeyCode::Char('G') | KeyCode::End => {
            if let AppPhase::Dashboard { active_panel, .. } = &app.phase {
                if *active_panel == Panel::Results {
                    app.go_to_bottom();
                }
            }
        }

        // View switching
        KeyCode::Tab => app.toggle_view(),

        // Panel switching (within table view)
        KeyCode::BackTab => {
            if let AppPhase::Dashboard { active_panel, .. } = &mut app.phase {
                *active_panel = match *active_panel {
                    Panel::Results => Panel::Details,
                    Panel::Details => Panel::Results,
                };
            }
        }

        // Tier filters
        KeyCode::Char('1') => app.toggle_tier(ConfidenceTier::Tier1Reliable),
        KeyCode::Char('2') => app.toggle_tier(ConfidenceTier::Tier2Probable),
        KeyCode::Char('3') => app.toggle_tier(ConfidenceTier::Tier3Speculative),
        KeyCode::Char('a') => app.show_all(),

        // Search
        KeyCode::Char('/') => {
            if let AppPhase::Dashboard {
                search_mode,
                search_query,
                ..
            } = &mut app.phase
            {
                *search_mode = true;
                search_query.clear();
            }
        }

        _ => {}
    }

    Ok(false)
}

/// Handle keyboard input while in search mode. Returns `true` if should quit.
fn handle_search_input(app: &mut App, code: KeyCode, modifiers: KeyModifiers) -> Result<bool> {
    match code {
        KeyCode::Esc => {
            if let AppPhase::Dashboard {
                search_mode,
                search_query,
                ..
            } = &mut app.phase
            {
                *search_mode = false;
                search_query.clear();
            }
            app.apply_filter();
        }
        KeyCode::Enter => {
            if let AppPhase::Dashboard { search_mode, .. } = &mut app.phase {
                *search_mode = false;
                // Keep the filter active, just exit input mode
            }
        }
        KeyCode::Backspace => {
            if let AppPhase::Dashboard { search_query, .. } = &mut app.phase {
                search_query.pop();
            }
            app.apply_filter();
        }
        KeyCode::Char('c') if modifiers.contains(KeyModifiers::CONTROL) => {
            return Ok(true);
        }
        KeyCode::Char(c) => {
            if let AppPhase::Dashboard { search_query, .. } = &mut app.phase {
                search_query.push(c);
            }
            app.apply_filter();
        }
        _ => {}
    }
    Ok(false)
}
