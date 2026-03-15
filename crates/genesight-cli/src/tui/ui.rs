use genesight_core::models::{ConfidenceTier, ResultCategory};
use ratatui::{
    layout::{Constraint, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span, Text},
    widgets::{
        Block, Borders, Cell, Paragraph, Row, Scrollbar, ScrollbarOrientation, ScrollbarState,
        Table, TableState, Wrap,
    },
    Frame,
};

use super::app::{App, AppPhase, DashboardView, Panel, ViewMode};

/// Main color constants.
const TIER1_COLOR: Color = Color::Green;
const TIER2_COLOR: Color = Color::Yellow;
const TIER3_COLOR: Color = Color::DarkGray;
const SELECTED_BG: Color = Color::Rgb(30, 60, 90);
const HEADER_BG: Color = Color::Rgb(20, 40, 70);
const BORDER_FOCUSED: Color = Color::Cyan;
const BORDER_UNFOCUSED: Color = Color::DarkGray;

/// Braille spinner frames for the loading animation.
const SPINNER_FRAMES: &[char] = &[
    '\u{280B}', '\u{2819}', '\u{2839}', '\u{2838}', '\u{283C}', '\u{2834}', '\u{2826}', '\u{2827}',
    '\u{2807}', '\u{280F}',
];

/// Color for a confidence tier.
fn tier_color(tier: ConfidenceTier) -> Color {
    match tier {
        ConfidenceTier::Tier1Reliable => TIER1_COLOR,
        ConfidenceTier::Tier2Probable => TIER2_COLOR,
        ConfidenceTier::Tier3Speculative => TIER3_COLOR,
    }
}

/// Short tier label.
fn tier_label(tier: ConfidenceTier) -> &'static str {
    match tier {
        ConfidenceTier::Tier1Reliable => "T1",
        ConfidenceTier::Tier2Probable => "T2",
        ConfidenceTier::Tier3Speculative => "T3",
    }
}

/// Draw the entire TUI -- dispatches to loading or dashboard based on phase.
pub fn draw(frame: &mut Frame, app: &App) {
    match &app.phase {
        AppPhase::Loading { .. } => draw_loading(frame, app),
        AppPhase::Dashboard { .. } => draw_dashboard(frame, app),
    }
}

// =============================================================================
// Loading Screen — DNA Double Helix Animation
// =============================================================================

/// Draw the loading screen: DNA helix spans the full screen as a background,
/// with the info box rendered on top in the center.
fn draw_loading(frame: &mut Frame, app: &App) {
    let area = frame.area();

    let (stages, current_stage, elapsed, spinner_tick) = match &app.phase {
        AppPhase::Loading {
            stages,
            current_stage,
            elapsed,
            spinner_tick,
        } => (stages, current_stage, elapsed, *spinner_tick),
        _ => return,
    };

    let spinner = SPINNER_FRAMES[spinner_tick % SPINNER_FRAMES.len()];
    let elapsed_secs = elapsed.elapsed().as_secs_f64();

    // --- Layer 1: DNA Helix background (full screen) ---
    let helix_lines = generate_dna_helix(area.width, area.height, spinner_tick);
    frame.render_widget(Paragraph::new(Text::from(helix_lines)), area);

    // --- Layer 2: Info box on top (centered) ---
    let box_width: u16 = 50.min(area.width.saturating_sub(4));
    let stage_count = stages.len() + 1;
    let box_height: u16 = ((8 + stage_count) as u16).min(area.height.saturating_sub(2));
    let box_x = area.x + area.width.saturating_sub(box_width) / 2;
    let box_y = area.y + area.height.saturating_sub(box_height) / 2;
    let box_rect = Rect::new(box_x, box_y, box_width, box_height);

    let mut lines: Vec<Line<'static>> = Vec::new();

    // Title
    lines.push(Line::from(Span::styled(
        "GeneSight",
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    )));
    lines.push(Line::from(Span::styled(
        "DNA Analysis",
        Style::default().fg(Color::White),
    )));
    lines.push(Line::raw(""));

    // Completed stages
    for (desc, completed) in stages {
        let icon = if *completed { "\u{2713}" } else { " " };
        let color = if *completed {
            Color::Green
        } else {
            Color::DarkGray
        };
        lines.push(Line::from(vec![
            Span::styled(format!(" {icon} "), Style::default().fg(color)),
            Span::styled(
                desc.clone(),
                Style::default().fg(if *completed {
                    Color::White
                } else {
                    Color::DarkGray
                }),
            ),
        ]));
    }

    // Current stage (with spinner)
    if !current_stage.is_empty() {
        lines.push(Line::from(vec![
            Span::styled(format!(" {spinner} "), Style::default().fg(Color::Cyan)),
            Span::styled(current_stage.clone(), Style::default().fg(Color::Cyan)),
        ]));
    }

    lines.push(Line::raw(""));

    // Elapsed time
    lines.push(Line::from(Span::styled(
        format!(" Elapsed: {elapsed_secs:.1}s"),
        Style::default().fg(Color::DarkGray),
    )));

    // Error message if present
    if let Some(ref err) = app.error_message {
        lines.push(Line::raw(""));
        lines.push(Line::from(Span::styled(
            format!(" Error: {err}"),
            Style::default()
                .fg(Color::Red)
                .add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(Span::styled(
            " Press q to quit",
            Style::default().fg(Color::DarkGray),
        )));
    }

    // Dark background block so text is readable over the helix
    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan))
        .title_alignment(ratatui::layout::Alignment::Center);

    let paragraph = Paragraph::new(Text::from(lines))
        .block(block)
        .style(Style::default().bg(Color::Black))
        .alignment(ratatui::layout::Alignment::Left);

    frame.render_widget(paragraph, box_rect);
}

/// Generate a DNA double helix animation frame.
///
/// Two sinusoidal strands weave around each other with base-pair
/// connections visible between them. The `tick` parameter advances
/// the rotation. Front strand is bright cyan, back strand is dim,
/// and base-pair rungs fade based on viewing angle.
fn generate_dna_helix(width: u16, height: u16, tick: usize) -> Vec<Line<'static>> {
    let w = width as usize;
    let h = height as usize;
    if w < 4 || h < 5 {
        return (0..h).map(|_| Line::raw("")).collect();
    }

    let phase = tick as f64 * 0.15;
    let center = (h - 1) as f64 / 2.0;
    // Helix occupies roughly 1/3 of the screen height — a tight band through the middle
    let amplitude = ((h as f64 / 3.0) - 0.5).max(2.0);
    let freq = std::f64::consts::TAU / 26.0; // one full twist per ~26 columns

    // Character palettes — DNA bases mixed with structural characters
    let strand1_chars: &[u8] = b"aAcCgGtT~=<>{}*+!?";
    let strand2_chars: &[u8] = b"tTgGcCaA[]()#@&%^~";

    // Style table: 0=bg  1=front-strand  2=back-strand  3=bright-rung  4=dim-rung
    let styles: [Style; 5] = [
        Style::default(),
        Style::default().fg(Color::Cyan),
        Style::default().fg(Color::Rgb(0, 80, 105)),
        Style::default().fg(Color::Rgb(0, 60, 80)),
        Style::default().fg(Color::Rgb(0, 35, 50)),
    ];

    // Build character grid: (char, style_index)
    let mut grid: Vec<Vec<(char, usize)>> = vec![vec![(' ', 0); w]; h];

    // Each column computes y-positions via sine and writes to multiple rows,
    // so a simple iterator-based loop would not be clearer.
    #[allow(clippy::needless_range_loop)]
    for x in 0..w {
        let xf = x as f64;
        let angle = xf * freq + phase;
        let sin_v = angle.sin();
        let cos_v = angle.cos();

        // Two strand positions — 180° apart
        let y1 = (center + amplitude * sin_v).round() as isize;
        let y2 = (center - amplitude * sin_v).round() as isize;

        // cos > 0 means strand 1 is "in front" (closer to the viewer)
        let strand1_front = cos_v > 0.0;
        let depth = cos_v.abs(); // 0.0 at crossings, 1.0 at max separation

        // --- Base-pair rungs between strands ---
        // Visible when strands face the viewer (depth > threshold).
        // Denser when more face-on, sparser when edge-on.
        if depth > 0.15 {
            let (top, bot) = if y1 < y2 { (y1, y2) } else { (y2, y1) };
            let gap = bot - top;
            if gap > 2 {
                // Rung density: face-on → every 2 cols, edge-on → every 4-5 cols
                let interval = if depth > 0.6 {
                    2
                } else if depth > 0.35 {
                    3
                } else {
                    5
                };
                if x % interval == 0 {
                    let rung_si = if depth > 0.5 { 3 } else { 4 };
                    for y in (top + 1)..bot {
                        if y >= 0 && (y as usize) < h {
                            let frac = (y - top) as f64 / gap as f64;
                            let c = if (0.35..0.65).contains(&frac) {
                                ':'
                            } else {
                                '.'
                            };
                            grid[y as usize][x] = (c, rung_si);
                        }
                    }
                }
            }
        }

        // --- Strand characters (varied texture) ---
        let ci1 = ((xf * 1.7 + phase * 2.3).abs() as usize) % strand1_chars.len();
        let ci2 = ((xf * 2.1 + phase * 1.9 + 11.0).abs() as usize) % strand2_chars.len();
        let c1 = strand1_chars[ci1] as char;
        let c2 = strand2_chars[ci2] as char;

        // Back strand first, then front (front overwrites at crossings)
        let (front_y, front_c, back_y, back_c) = if strand1_front {
            (y1, c1, y2, c2)
        } else {
            (y2, c2, y1, c1)
        };

        if back_y >= 0 && (back_y as usize) < h {
            grid[back_y as usize][x] = (back_c, 2);
        }
        if front_y >= 0 && (front_y as usize) < h {
            grid[front_y as usize][x] = (front_c, 1);
        }
    }

    // Convert grid → Lines, grouping consecutive same-style chars into one Span
    grid.into_iter()
        .map(|row| {
            let mut spans: Vec<Span<'static>> = Vec::new();
            let mut buf = String::new();
            let mut cur: usize = 0;

            for (ch, si) in row {
                if si != cur {
                    if !buf.is_empty() {
                        spans.push(Span::styled(
                            std::mem::take(&mut buf),
                            styles[cur],
                        ));
                    }
                    cur = si;
                }
                buf.push(ch);
            }
            if !buf.is_empty() {
                spans.push(Span::styled(buf, styles[cur]));
            }

            Line::from(spans)
        })
        .collect()
}

// =============================================================================
// Summary View
// =============================================================================

/// Draw the summary view — a high-level overview of findings.
fn draw_summary(frame: &mut Frame, app: &App, area: Rect) {
    let report = match &app.report {
        Some(r) => r,
        None => return,
    };

    let scroll = app.summary_scroll();
    let (t1, t2, t3) = app.tier_counts();
    let category_counts = app.category_counts();

    let mut lines: Vec<Line<'static>> = Vec::new();

    // --- Key Findings (Tier 1) ---
    lines.push(Line::raw(""));
    lines.push(Line::from(vec![
        Span::styled(
            " Clinically Actionable Findings ",
            Style::default()
                .fg(Color::White)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(
            format!("  ({t1} Tier 1 results)"),
            Style::default().fg(TIER1_COLOR),
        ),
    ]));
    lines.push(Line::from(Span::styled(
        " \u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}",
        Style::default().fg(Color::DarkGray),
    )));

    let tier1_results: Vec<&genesight_core::models::ScoredResult> = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier1Reliable)
        .collect();

    if tier1_results.is_empty() {
        lines.push(Line::from(Span::styled(
            "   No Tier 1 (clinically reliable) findings.",
            Style::default().fg(Color::DarkGray),
        )));
    } else {
        // Group Tier 1 by category
        let monogenic: Vec<_> = tier1_results
            .iter()
            .filter(|r| r.category == ResultCategory::MonogenicDisease)
            .collect();
        let pharma: Vec<_> = tier1_results
            .iter()
            .filter(|r| r.category == ResultCategory::Pharmacogenomics)
            .collect();
        let carrier: Vec<_> = tier1_results
            .iter()
            .filter(|r| r.category == ResultCategory::CarrierStatus)
            .collect();

        if !monogenic.is_empty() {
            lines.push(Line::raw(""));
            lines.push(Line::from(Span::styled(
                format!("   Disease Risk ({} variants)", monogenic.len()),
                Style::default()
                    .fg(Color::Red)
                    .add_modifier(Modifier::BOLD),
            )));
            for r in &monogenic {
                let gene = app.gene_name(r);
                let rsid = r.variant.variant.rsid.as_deref().unwrap_or("\u{2014}");
                let conditions = r
                    .variant
                    .clinvar
                    .as_ref()
                    .map(|cv| {
                        if cv.conditions.is_empty() {
                            "unspecified".to_string()
                        } else {
                            cv.conditions.join(", ")
                        }
                    })
                    .unwrap_or_default();
                let stars_str = r
                    .variant
                    .clinvar
                    .as_ref()
                    .map(|cv| format!(" ({})", stars(cv.review_stars)))
                    .unwrap_or_default();
                lines.push(Line::from(vec![
                    Span::styled("   \u{25cf} ", Style::default().fg(Color::Red)),
                    Span::styled(
                        gene,
                        Style::default()
                            .fg(Color::White)
                            .add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(format!(" ({rsid})"), Style::default().fg(Color::DarkGray)),
                    Span::styled(stars_str, Style::default().fg(Color::Yellow)),
                ]));
                if !conditions.is_empty() {
                    lines.push(Line::from(Span::styled(
                        format!("     {conditions}"),
                        Style::default().fg(Color::Gray),
                    )));
                }
            }
        }

        if !pharma.is_empty() {
            lines.push(Line::raw(""));
            lines.push(Line::from(Span::styled(
                format!("   Drug Interactions ({} variants)", pharma.len()),
                Style::default()
                    .fg(Color::Magenta)
                    .add_modifier(Modifier::BOLD),
            )));
            for r in &pharma {
                let gene = app.gene_name(r);
                let drug = r
                    .variant
                    .pharmacogenomics
                    .as_ref()
                    .map(|p| p.drug.as_str())
                    .unwrap_or("unknown");
                let phenotype = r
                    .variant
                    .pharmacogenomics
                    .as_ref()
                    .and_then(|p| p.phenotype_category.as_deref())
                    .unwrap_or("");
                lines.push(Line::from(vec![
                    Span::styled("   \u{25cf} ", Style::default().fg(Color::Magenta)),
                    Span::styled(
                        gene,
                        Style::default()
                            .fg(Color::White)
                            .add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(
                        format!(" \u{2192} {drug}"),
                        Style::default().fg(Color::White),
                    ),
                ]));
                if !phenotype.is_empty() {
                    lines.push(Line::from(Span::styled(
                        format!("     {phenotype}"),
                        Style::default().fg(Color::Gray),
                    )));
                }
            }
        }

        if !carrier.is_empty() {
            lines.push(Line::raw(""));
            lines.push(Line::from(Span::styled(
                format!("   Carrier Status ({} variants)", carrier.len()),
                Style::default()
                    .fg(Color::Blue)
                    .add_modifier(Modifier::BOLD),
            )));
            // Just show count — carrier status is less urgent
            let benign_count = carrier
                .iter()
                .filter(|r| {
                    r.variant
                        .clinvar
                        .as_ref()
                        .is_some_and(|cv| cv.significance.to_lowercase().contains("benign"))
                })
                .count();
            if benign_count > 0 {
                lines.push(Line::from(Span::styled(
                    format!("     {benign_count} benign classifications"),
                    Style::default().fg(Color::Gray),
                )));
            }
        }
    }

    // --- Category Breakdown ---
    lines.push(Line::raw(""));
    lines.push(Line::raw(""));
    lines.push(Line::from(Span::styled(
        " Category Breakdown",
        Style::default()
            .fg(Color::White)
            .add_modifier(Modifier::BOLD),
    )));
    lines.push(Line::from(Span::styled(
        " \u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}",
        Style::default().fg(Color::DarkGray),
    )));
    lines.push(Line::raw(""));

    let max_count = category_counts
        .iter()
        .map(|(_, c)| *c)
        .max()
        .unwrap_or(1);
    let bar_max_width: usize = 30;

    for (cat, count) in &category_counts {
        let bar_len = (*count as f64 / max_count as f64 * bar_max_width as f64).ceil() as usize;
        let bar = "\u{2588}".repeat(bar_len);
        let cat_color = category_color(*cat);
        let cat_label = format!("{cat}");
        let padding = 22_usize.saturating_sub(cat_label.len());

        lines.push(Line::from(vec![
            Span::styled(format!("   {cat_label}"), Style::default().fg(cat_color)),
            Span::raw(" ".repeat(padding)),
            Span::styled(
                format!("{count:>6}"),
                Style::default()
                    .fg(Color::White)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw("  "),
            Span::styled(bar, Style::default().fg(cat_color)),
        ]));
    }

    // --- Tier Distribution ---
    lines.push(Line::raw(""));
    lines.push(Line::raw(""));
    lines.push(Line::from(Span::styled(
        " Confidence Distribution",
        Style::default()
            .fg(Color::White)
            .add_modifier(Modifier::BOLD),
    )));
    lines.push(Line::from(Span::styled(
        " \u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}",
        Style::default().fg(Color::DarkGray),
    )));
    lines.push(Line::raw(""));

    let total = (t1 + t2 + t3).max(1);
    let pct1 = t1 as f64 / total as f64 * 100.0;
    let pct2 = t2 as f64 / total as f64 * 100.0;
    let pct3 = t3 as f64 / total as f64 * 100.0;

    for (label, count, pct, color, desc) in [
        ("Tier 1", t1, pct1, TIER1_COLOR, "Clinically reliable (>95%)"),
        ("Tier 2", t2, pct2, TIER2_COLOR, "Probable (60-85%)"),
        ("Tier 3", t3, pct3, TIER3_COLOR, "Speculative (50-65%)"),
    ] {
        let bar_len = (pct / 100.0 * bar_max_width as f64).ceil().max(if count > 0 { 1.0 } else { 0.0 }) as usize;
        let bar = "\u{2588}".repeat(bar_len);
        lines.push(Line::from(vec![
            Span::styled(format!("   {label}"), Style::default().fg(color).add_modifier(Modifier::BOLD)),
            Span::raw("   "),
            Span::styled(format!("{count:>6}"), Style::default().fg(Color::White).add_modifier(Modifier::BOLD)),
            Span::styled(format!(" ({pct:.1}%)"), Style::default().fg(Color::DarkGray)),
            Span::raw("  "),
            Span::styled(bar, Style::default().fg(color)),
        ]));
        lines.push(Line::from(Span::styled(
            format!("              {desc}"),
            Style::default().fg(Color::DarkGray),
        )));
    }

    // --- Notable Tier 2 highlights ---
    let tier2_notable: Vec<_> = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier2Probable)
        .take(8)
        .collect();

    if !tier2_notable.is_empty() {
        lines.push(Line::raw(""));
        lines.push(Line::raw(""));
        lines.push(Line::from(vec![
            Span::styled(
                " Notable Tier 2 Findings ",
                Style::default()
                    .fg(Color::White)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(
                format!("  ({t2} total)"),
                Style::default().fg(TIER2_COLOR),
            ),
        ]));
        lines.push(Line::from(Span::styled(
            " \u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}\u{2500}",
            Style::default().fg(Color::DarkGray),
        )));
        lines.push(Line::raw(""));

        for r in &tier2_notable {
            let gene = app.gene_name(r);
            let cat = short_category(r.category);
            let summary = truncate_str(&r.summary, 50);
            lines.push(Line::from(vec![
                Span::styled("   \u{25cb} ", Style::default().fg(TIER2_COLOR)),
                Span::styled(
                    format!("[{cat}]"),
                    Style::default().fg(Color::DarkGray),
                ),
                Span::raw(" "),
                Span::styled(gene, Style::default().fg(Color::White).add_modifier(Modifier::BOLD)),
                Span::styled(format!(" \u{2014} {summary}"), Style::default().fg(Color::Gray)),
            ]));
        }
    }

    // --- Tip ---
    lines.push(Line::raw(""));
    lines.push(Line::raw(""));
    lines.push(Line::from(vec![
        Span::styled(
            "   Tip: ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(
            format!(
                "Most findings ({t3}) are Tier 3 GWAS associations with small effect sizes. ",
            ),
            Style::default().fg(Color::DarkGray),
        ),
    ]));
    lines.push(Line::from(Span::styled(
        "   Press [Tab] to switch to the full table view for detailed exploration.",
        Style::default().fg(Color::DarkGray),
    )));
    lines.push(Line::raw(""));

    let content_height = lines.len();
    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan))
        .title(" Summary ");

    let paragraph = Paragraph::new(Text::from(lines))
        .block(block)
        .wrap(Wrap { trim: false })
        .scroll((scroll, 0));

    frame.render_widget(paragraph, area);

    // Scrollbar
    let visible_height = area.height.saturating_sub(2) as usize;
    if content_height > visible_height {
        let mut scrollbar_state = ScrollbarState::new(content_height).position(scroll as usize);
        frame.render_stateful_widget(
            Scrollbar::new(ScrollbarOrientation::VerticalRight)
                .begin_symbol(None)
                .end_symbol(None),
            area.inner(ratatui::layout::Margin {
                horizontal: 0,
                vertical: 1,
            }),
            &mut scrollbar_state,
        );
    }
}

/// Color for a result category.
fn category_color(cat: ResultCategory) -> Color {
    match cat {
        ResultCategory::MonogenicDisease => Color::Red,
        ResultCategory::CarrierStatus => Color::Blue,
        ResultCategory::Pharmacogenomics => Color::Magenta,
        ResultCategory::PolygenicRiskScore => Color::Yellow,
        ResultCategory::PhysicalTrait => Color::Cyan,
        ResultCategory::ComplexTrait => Color::DarkGray,
        ResultCategory::Ancestry => Color::White,
    }
}

// =============================================================================
// Dashboard
// =============================================================================

/// Draw the full dashboard (header, stats, body, footer).
fn draw_dashboard(frame: &mut Frame, app: &App) {
    let outer = frame.area();

    // Main vertical layout: header, stats, body, footer
    let vertical = Layout::vertical([
        Constraint::Length(1), // header
        Constraint::Length(2), // stats bar
        Constraint::Min(8),    // body (results + details OR summary)
        Constraint::Length(3), // footer
    ])
    .split(outer);

    draw_header(frame, app, vertical[0]);
    draw_stats_bar(frame, app, vertical[1]);

    match app.view_mode() {
        ViewMode::Summary => draw_summary(frame, app, vertical[2]),
        ViewMode::Table => draw_body(frame, app, vertical[2]),
    }

    draw_footer(frame, app, vertical[3]);
}

/// Draw the title header bar with view mode tabs.
fn draw_header(frame: &mut Frame, app: &App, area: Rect) {
    let view = app.view_mode();
    let (summary_style, table_style) = match view {
        ViewMode::Summary => (
            Style::default()
                .fg(Color::Cyan)
                .bg(HEADER_BG)
                .add_modifier(Modifier::BOLD),
            Style::default().fg(Color::DarkGray).bg(HEADER_BG),
        ),
        ViewMode::Table => (
            Style::default().fg(Color::DarkGray).bg(HEADER_BG),
            Style::default()
                .fg(Color::Cyan)
                .bg(HEADER_BG)
                .add_modifier(Modifier::BOLD),
        ),
    };

    let tabs_width: u16 = 55; // approximate width of fixed content
    let padding_width = area.width.saturating_sub(tabs_width) as usize;
    let title = Line::from(vec![
        Span::styled(
            " GeneSight ",
            Style::default()
                .fg(Color::White)
                .bg(HEADER_BG)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(" Summary ", summary_style),
        Span::styled(" | ", Style::default().fg(Color::DarkGray).bg(HEADER_BG)),
        Span::styled(" Table ", table_style),
        Span::styled(
            " ".repeat(padding_width.max(1)),
            Style::default().bg(HEADER_BG),
        ),
        Span::styled(
            "[Tab] Switch  [q] Quit ",
            Style::default().fg(Color::DarkGray).bg(HEADER_BG),
        ),
    ]);
    frame.render_widget(
        Paragraph::new(title).style(Style::default().bg(HEADER_BG)),
        area,
    );
}

/// Draw the statistics summary bar.
fn draw_stats_bar(frame: &mut Frame, app: &App, area: Rect) {
    let report = match &app.report {
        Some(r) => r,
        None => return,
    };
    let (t1, t2, t3) = app.tier_counts();

    let dv = match app.dashboard_state() {
        Some(s) => s,
        None => return,
    };
    let tier_filter = dv.tier_filter;

    let line1 = Line::from(vec![
        Span::raw("  Variants: "),
        Span::styled(
            format!("{}", report.total_variants),
            Style::default()
                .fg(Color::White)
                .add_modifier(Modifier::BOLD),
        ),
        Span::raw("  |  Annotated: "),
        Span::styled(
            format!("{}", report.annotated_variants),
            Style::default()
                .fg(Color::White)
                .add_modifier(Modifier::BOLD),
        ),
        Span::raw("  |  Findings: "),
        Span::styled(
            format!("{}", report.results.len()),
            Style::default()
                .fg(Color::White)
                .add_modifier(Modifier::BOLD),
        ),
    ]);

    let mut line2_spans = vec![
        Span::raw("  "),
        Span::styled("\u{25a0} ", Style::default().fg(TIER1_COLOR)),
        Span::styled(format!("Tier 1: {t1}"), Style::default().fg(TIER1_COLOR)),
        Span::raw("  "),
        Span::styled("\u{25a0} ", Style::default().fg(TIER2_COLOR)),
        Span::styled(format!("Tier 2: {t2}"), Style::default().fg(TIER2_COLOR)),
        Span::raw("  "),
        Span::styled("\u{25a0} ", Style::default().fg(TIER3_COLOR)),
        Span::styled(format!("Tier 3: {t3}"), Style::default().fg(TIER3_COLOR)),
    ];
    if let Some(filter) = tier_filter {
        line2_spans.push(Span::styled(
            format!("    [Filtered: {} only]", tier_label(filter)),
            Style::default()
                .fg(tier_color(filter))
                .add_modifier(Modifier::BOLD),
        ));
    }
    let line2 = Line::from(line2_spans);

    let text = Text::from(vec![line1, line2]);
    frame.render_widget(Paragraph::new(text), area);
}

/// Draw the main body: results table on the left, details panel on the right.
fn draw_body(frame: &mut Frame, app: &App, area: Rect) {
    let horizontal =
        Layout::horizontal([Constraint::Percentage(58), Constraint::Percentage(42)]).split(area);

    draw_results_table(frame, app, horizontal[0]);
    draw_details_panel(frame, app, horizontal[1]);
}

/// Draw the scrollable results table with virtual scrolling.
///
/// Instead of creating Row widgets for all 200K+ filtered results, this only
/// materializes rows for the visible window (typically ~30 rows), making
/// rendering O(visible_rows) instead of O(total_filtered).
fn draw_results_table(frame: &mut Frame, app: &App, area: Rect) {
    let report = match &app.report {
        Some(r) => r,
        None => return,
    };
    let DashboardView {
        filtered_indices,
        selected_index,
        active_panel,
        ..
    } = match app.dashboard_state() {
        Some(s) => s,
        None => return,
    };

    let focused = active_panel == Panel::Results;
    let border_color = if focused {
        BORDER_FOCUSED
    } else {
        BORDER_UNFOCUSED
    };

    let header_cells = ["rsID", "Gene", "Category", "Summary", "Tier"]
        .iter()
        .map(|h| {
            Cell::from(*h).style(
                Style::default()
                    .fg(Color::White)
                    .add_modifier(Modifier::BOLD),
            )
        });
    let header = Row::new(header_cells)
        .style(Style::default().bg(HEADER_BG))
        .height(1);

    // Virtual scrolling: only render visible rows
    // area.height - 3 accounts for top border (1) + header row (1) + bottom border (1)
    let visible_rows = area.height.saturating_sub(3) as usize;
    let total = filtered_indices.len();

    // Calculate the window start to keep the selection roughly centered
    let (window_start, relative_selection) = if total == 0 {
        (0, 0)
    } else {
        let half = visible_rows / 2;
        let start = if selected_index <= half {
            0
        } else if selected_index + half >= total {
            total.saturating_sub(visible_rows)
        } else {
            selected_index.saturating_sub(half)
        };
        let end = (start + visible_rows).min(total);
        let clamped_start = if end == total {
            total.saturating_sub(visible_rows)
        } else {
            start
        };
        (clamped_start, selected_index - clamped_start)
    };

    let window_end = (window_start + visible_rows).min(total);

    let rows: Vec<Row> = filtered_indices[window_start..window_end]
        .iter()
        .map(|&idx| {
            let result = &report.results[idx];
            let tc = tier_color(result.tier);
            let rsid = result.variant.variant.rsid.as_deref().unwrap_or("\u{2014}");
            let gene = app.gene_name(result);
            let category = short_category(result.category);
            let summary = truncate_str(&result.summary, 30);
            let tier = tier_label(result.tier);

            Row::new(vec![
                Cell::from(rsid.to_string()),
                Cell::from(gene),
                Cell::from(category),
                Cell::from(summary),
                Cell::from(tier).style(Style::default().fg(tc).add_modifier(Modifier::BOLD)),
            ])
            .style(Style::default().fg(Color::White))
        })
        .collect();

    let widths = [
        Constraint::Min(12),
        Constraint::Min(8),
        Constraint::Min(10),
        Constraint::Min(14),
        Constraint::Length(4),
    ];

    let title = if focused {
        format!(
            " Results [\u{2191}\u{2193}] ({}/{}) ",
            if total > 0 { selected_index + 1 } else { 0 },
            total
        )
    } else {
        format!(" Results ({total}) ")
    };
    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(border_color))
        .title(title);

    let table = Table::new(rows, widths)
        .header(header)
        .block(block)
        .row_highlight_style(
            Style::default()
                .bg(SELECTED_BG)
                .add_modifier(Modifier::BOLD),
        )
        .highlight_symbol("\u{25ba} ");

    // Selection is relative to the visible window
    let mut state = TableState::default();
    if !filtered_indices.is_empty() {
        state.select(Some(relative_selection));
    }

    frame.render_stateful_widget(table, area, &mut state);

    // Render scrollbar for the results table (tracks overall position)
    if total > visible_rows {
        let mut scrollbar_state = ScrollbarState::new(total).position(selected_index);
        frame.render_stateful_widget(
            Scrollbar::new(ScrollbarOrientation::VerticalRight)
                .begin_symbol(None)
                .end_symbol(None),
            area.inner(ratatui::layout::Margin {
                horizontal: 0,
                vertical: 1,
            }),
            &mut scrollbar_state,
        );
    }
}

/// Draw the detail panel for the selected result.
fn draw_details_panel(frame: &mut Frame, app: &App, area: Rect) {
    let DashboardView {
        detail_scroll,
        active_panel,
        ..
    } = match app.dashboard_state() {
        Some(s) => s,
        None => return,
    };

    let focused = active_panel == Panel::Details;
    let border_color = if focused {
        BORDER_FOCUSED
    } else {
        BORDER_UNFOCUSED
    };

    let title = if focused {
        " Details [\u{2191}\u{2193} scroll] "
    } else {
        " Details "
    };
    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(border_color))
        .title(title);

    let Some(result) = app.selected_result() else {
        let empty = Paragraph::new("No result selected")
            .block(block)
            .style(Style::default().fg(Color::DarkGray));
        frame.render_widget(empty, area);
        return;
    };

    let wrap_width = area.width.saturating_sub(6) as usize;
    let mut lines: Vec<Line<'static>> = Vec::new();

    // Gene
    let gene = app.gene_name(result);
    lines.push(Line::from(vec![
        Span::styled(
            "Gene: ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(gene, Style::default().fg(Color::White)),
    ]));

    // rsID
    let rsid = result
        .variant
        .variant
        .rsid
        .clone()
        .unwrap_or_else(|| "\u{2014}".to_string());
    lines.push(Line::from(vec![
        Span::styled(
            "rsID: ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(rsid, Style::default().fg(Color::White)),
    ]));

    // Genotype
    lines.push(Line::from(vec![
        Span::styled(
            "Genotype: ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(
            format!("{}", result.variant.variant.genotype),
            Style::default().fg(Color::White),
        ),
    ]));

    // Chromosome / Position
    lines.push(Line::from(vec![
        Span::styled(
            "Location: ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(
            format!(
                "chr{}:{}",
                result.variant.variant.chromosome, result.variant.variant.position
            ),
            Style::default().fg(Color::White),
        ),
    ]));

    // Tier
    let tc = tier_color(result.tier);
    lines.push(Line::from(vec![
        Span::styled(
            "Confidence: ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(format!("{}", result.tier), Style::default().fg(tc)),
    ]));

    // Category
    lines.push(Line::from(vec![
        Span::styled(
            "Category: ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(
            format!("{}", result.category),
            Style::default().fg(Color::White),
        ),
    ]));

    lines.push(Line::raw(""));

    // ClinVar details
    if let Some(ref cv) = result.variant.clinvar {
        lines.push(section_header("ClinVar"));
        let classification = format!("{} ({})", cv.significance, stars(cv.review_stars));
        lines.push(detail_line("  Classification", &classification));
        if !cv.conditions.is_empty() {
            lines.push(detail_line("  Conditions", ""));
            for cond in &cv.conditions {
                lines.push(Line::from(Span::styled(
                    format!("    - {cond}"),
                    Style::default().fg(Color::White),
                )));
            }
        }
        lines.push(Line::raw(""));
    }

    // Allele frequency
    if let Some(ref freq) = result.variant.frequency {
        lines.push(section_header("Allele Frequency"));
        let total_str = format_freq(freq.af_total);
        lines.push(detail_line("  Total", &total_str));
        if let Some(af) = freq.af_afr {
            let s = format_freq(af);
            lines.push(detail_line("  AFR", &s));
        }
        if let Some(af) = freq.af_amr {
            let s = format_freq(af);
            lines.push(detail_line("  AMR", &s));
        }
        if let Some(af) = freq.af_eas {
            let s = format_freq(af);
            lines.push(detail_line("  EAS", &s));
        }
        if let Some(af) = freq.af_eur {
            let s = format_freq(af);
            lines.push(detail_line("  EUR", &s));
        }
        if let Some(af) = freq.af_sas {
            let s = format_freq(af);
            lines.push(detail_line("  SAS", &s));
        }
        lines.push(detail_line("  Source", &freq.source));
        lines.push(Line::raw(""));
    }

    // Pharmacogenomics
    if let Some(ref pharma) = result.variant.pharmacogenomics {
        lines.push(section_header("Pharmacogenomics"));
        lines.push(detail_line("  Gene", &pharma.gene));
        lines.push(detail_line("  Drug", &pharma.drug));
        lines.push(detail_line("  Evidence", &pharma.evidence_level));
        if let Some(ref pheno) = pharma.phenotype_category {
            lines.push(detail_line("  Phenotype", pheno));
        }
        if let Some(ref rec) = pharma.clinical_recommendation {
            lines.push(detail_line("  Recommendation", ""));
            for chunk in wrap_text(rec, wrap_width) {
                lines.push(Line::from(Span::styled(
                    format!("    {chunk}"),
                    Style::default().fg(Color::White),
                )));
            }
        }
        lines.push(Line::raw(""));
    }

    // GWAS hits
    if !result.variant.gwas_hits.is_empty() {
        lines.push(section_header("GWAS Associations"));
        for hit in &result.variant.gwas_hits {
            lines.push(detail_line("  Trait", &hit.trait_name));
            let pval = format!("{:.2e}", hit.p_value);
            lines.push(detail_line("  p-value", &pval));
            if let Some(or) = hit.odds_ratio {
                let or_str = format!("{or:.3}");
                lines.push(detail_line("  Odds ratio", &or_str));
            }
            if let Some(ref allele) = hit.risk_allele {
                lines.push(detail_line("  Risk allele", allele));
            }
            if let Some(ref pmid) = hit.pubmed_id {
                lines.push(detail_line("  PubMed", pmid));
            }
            lines.push(Line::raw(""));
        }
    }

    // SNPedia summary
    if let Some(ref snp) = result.variant.snpedia {
        lines.push(section_header("SNPedia"));
        let mag = format!("{:.1}", snp.magnitude);
        lines.push(detail_line("  Magnitude", &mag));
        if let Some(ref repute) = snp.repute {
            lines.push(detail_line("  Repute", repute));
        }
        for chunk in wrap_text(&snp.summary, wrap_width) {
            lines.push(Line::from(Span::styled(
                format!("    {chunk}"),
                Style::default().fg(Color::White),
            )));
        }
        lines.push(Line::raw(""));
    }

    // Summary & Details from the scored result
    lines.push(section_header("Summary"));
    for chunk in wrap_text(&result.summary, wrap_width) {
        lines.push(Line::from(Span::styled(
            format!("  {chunk}"),
            Style::default().fg(Color::White),
        )));
    }
    lines.push(Line::raw(""));
    lines.push(section_header("Details"));
    for chunk in wrap_text(&result.details, wrap_width) {
        lines.push(Line::from(Span::styled(
            format!("  {chunk}"),
            Style::default().fg(Color::White),
        )));
    }

    // Estimate content height for scrollbar before consuming lines
    let content_height = lines.len();

    let text = Text::from(lines);
    let paragraph = Paragraph::new(text)
        .block(block)
        .wrap(Wrap { trim: false })
        .scroll((detail_scroll, 0));

    frame.render_widget(paragraph, area);

    // Detail scrollbar
    let visible_height = area.height.saturating_sub(2) as usize;
    if content_height > visible_height {
        let mut scrollbar_state =
            ScrollbarState::new(content_height).position(detail_scroll as usize);
        frame.render_stateful_widget(
            Scrollbar::new(ScrollbarOrientation::VerticalRight)
                .begin_symbol(None)
                .end_symbol(None),
            area.inner(ratatui::layout::Margin {
                horizontal: 0,
                vertical: 1,
            }),
            &mut scrollbar_state,
        );
    }
}

/// Draw the footer with keybinding hints and optional search bar.
fn draw_footer(frame: &mut Frame, app: &App, area: Rect) {
    let dv = match app.dashboard_state() {
        Some(s) => s,
        None => return,
    };

    if dv.search_mode {
        let search_line = Line::from(vec![
            Span::styled(
                " Search: ",
                Style::default()
                    .fg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(dv.search_query.to_string(), Style::default().fg(Color::White)),
            Span::styled(
                "_",
                Style::default()
                    .fg(Color::White)
                    .add_modifier(Modifier::SLOW_BLINK),
            ),
            Span::styled(
                "  [Enter] confirm  [Esc] cancel",
                Style::default().fg(Color::DarkGray),
            ),
        ]);
        let block = Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(BORDER_FOCUSED))
            .title(" Search ");
        frame.render_widget(Paragraph::new(search_line).block(block), area);
    } else if dv.view_mode == ViewMode::Summary {
        // Summary view footer
        let line1 = Line::from(vec![
            Span::styled(" [Tab]", Style::default().fg(Color::Cyan)),
            Span::raw(" Full Table  "),
            Span::styled("[1]", Style::default().fg(TIER1_COLOR)),
            Span::raw(" Tier 1  "),
            Span::styled("[2]", Style::default().fg(TIER2_COLOR)),
            Span::raw(" Tier 2  "),
            Span::styled("[3]", Style::default().fg(TIER3_COLOR)),
            Span::raw(" Tier 3  "),
            Span::styled("[/]", Style::default().fg(Color::White)),
            Span::raw(" Search"),
        ]);
        let line2 = Line::from(vec![
            Span::styled(" [j/\u{2193}]", Style::default().fg(Color::White)),
            Span::raw(" Scroll Down  "),
            Span::styled("[k/\u{2191}]", Style::default().fg(Color::White)),
            Span::raw(" Scroll Up  "),
            Span::styled("[g]", Style::default().fg(Color::White)),
            Span::raw(" Top  "),
            Span::styled("[q]", Style::default().fg(Color::Red)),
            Span::raw(" Quit"),
        ]);
        let block = Block::default()
            .borders(Borders::TOP)
            .border_style(Style::default().fg(Color::DarkGray));
        frame.render_widget(Paragraph::new(vec![line1, line2]).block(block), area);
    } else {
        // Table view footer
        let filter_indicator = match dv.tier_filter {
            Some(ConfidenceTier::Tier1Reliable) => "[1]* ",
            Some(ConfidenceTier::Tier2Probable) => "[2]* ",
            Some(ConfidenceTier::Tier3Speculative) => "[3]* ",
            None => "",
        };

        let line1 = Line::from(vec![
            Span::styled(" [Tab]", Style::default().fg(Color::Cyan)),
            Span::raw(" Summary  "),
            Span::styled("[1]", Style::default().fg(TIER1_COLOR)),
            Span::raw(" Tier 1  "),
            Span::styled("[2]", Style::default().fg(TIER2_COLOR)),
            Span::raw(" Tier 2  "),
            Span::styled("[3]", Style::default().fg(TIER3_COLOR)),
            Span::raw(" Tier 3  "),
            Span::styled("[a]", Style::default().fg(Color::White)),
            Span::raw(" All  "),
            Span::styled("[/]", Style::default().fg(Color::White)),
            Span::raw(" Search  "),
            Span::styled(filter_indicator, Style::default().fg(Color::Cyan)),
        ]);
        let line2 = Line::from(vec![
            Span::styled(" [j/\u{2193}]", Style::default().fg(Color::White)),
            Span::raw(" Down  "),
            Span::styled("[k/\u{2191}]", Style::default().fg(Color::White)),
            Span::raw(" Up  "),
            Span::styled("[g]", Style::default().fg(Color::White)),
            Span::raw(" Top  "),
            Span::styled("[G]", Style::default().fg(Color::White)),
            Span::raw(" Bottom  "),
            Span::styled("[Shift+Tab]", Style::default().fg(Color::White)),
            Span::raw(" Switch panel  "),
            Span::styled("[q]", Style::default().fg(Color::Red)),
            Span::raw(" Quit"),
        ]);
        let block = Block::default()
            .borders(Borders::TOP)
            .border_style(Style::default().fg(Color::DarkGray));
        frame.render_widget(Paragraph::new(vec![line1, line2]).block(block), area);
    }
}

// --- Helpers -----------------------------------------------------------------

/// Create a styled section header line.
fn section_header(title: &str) -> Line<'static> {
    Line::from(Span::styled(
        format!("--- {title} ---"),
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    ))
}

/// Create a "label: value" detail line. Takes string slices, clones internally
/// to produce an owned `Line<'static>` that can be collected into a Vec.
fn detail_line(label: &str, value: &str) -> Line<'static> {
    Line::from(vec![
        Span::styled(format!("{label}: "), Style::default().fg(Color::Gray)),
        Span::styled(value.to_string(), Style::default().fg(Color::White)),
    ])
}

/// Generate a star string for ClinVar review status.
fn stars(count: u8) -> String {
    let filled = "\u{2605}".repeat(count as usize);
    let empty = "\u{2606}".repeat(4_usize.saturating_sub(count as usize));
    format!("{filled}{empty}")
}

/// Format a frequency value nicely.
fn format_freq(f: f64) -> String {
    if f < 0.001 {
        format!("{f:.4e}")
    } else {
        format!("{f:.4}")
    }
}

/// Truncate a string to a maximum length, adding "..." if needed.
fn truncate_str(s: &str, max: usize) -> String {
    if s.len() <= max {
        s.to_string()
    } else {
        let truncated: String = s.chars().take(max.saturating_sub(3)).collect();
        format!("{truncated}...")
    }
}

/// Short category label for the table.
fn short_category(cat: genesight_core::models::ResultCategory) -> String {
    use genesight_core::models::ResultCategory;
    match cat {
        ResultCategory::MonogenicDisease => "Disease".to_string(),
        ResultCategory::CarrierStatus => "Carrier".to_string(),
        ResultCategory::Pharmacogenomics => "Pharma".to_string(),
        ResultCategory::PolygenicRiskScore => "PRS".to_string(),
        ResultCategory::PhysicalTrait => "Trait".to_string(),
        ResultCategory::ComplexTrait => "Complex".to_string(),
        ResultCategory::Ancestry => "Ancestry".to_string(),
    }
}

/// Simple word-wrap for a string to fit within a given width.
fn wrap_text(s: &str, width: usize) -> Vec<String> {
    if width == 0 {
        return vec![s.to_string()];
    }
    let mut lines = Vec::new();
    for paragraph in s.split('\n') {
        let mut current_line = String::new();
        for word in paragraph.split_whitespace() {
            if current_line.is_empty() {
                current_line = word.to_string();
            } else if current_line.len() + 1 + word.len() <= width {
                current_line.push(' ');
                current_line.push_str(word);
            } else {
                lines.push(current_line);
                current_line = word.to_string();
            }
        }
        lines.push(current_line);
    }
    if lines.is_empty() {
        lines.push(String::new());
    }
    lines
}
