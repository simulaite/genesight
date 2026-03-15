use std::path::Path;

use anyhow::{Context, Result};
use genesight_core::models::Report;
use genesight_core::report::{self, OutputFormat};

/// Export the report as an HTML file.
pub fn export_html(report: &Report, path: &Path) -> Result<()> {
    let html = report::render(report, OutputFormat::Html)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    std::fs::write(path, html.as_bytes())
        .with_context(|| format!("writing HTML to {}", path.display()))?;
    Ok(())
}

/// Export the report as a PDF file using printpdf.
pub fn export_pdf(report: &Report, path: &Path) -> Result<()> {
    use printpdf::*;

    let (doc, page1, layer1) = PdfDocument::new(
        "GeneSight Report",
        Mm(210.0),
        Mm(297.0),
        "Layer 1",
    );

    let font = doc.add_builtin_font(BuiltinFont::Helvetica)
        .map_err(|e| anyhow::anyhow!("font error: {e}"))?;
    let font_bold = doc.add_builtin_font(BuiltinFont::HelveticaBold)
        .map_err(|e| anyhow::anyhow!("font error: {e}"))?;

    let margin_left = Mm(20.0);
    let _page_width = Mm(170.0);
    let mut y = Mm(277.0);
    let line_height = Mm(5.0);

    // Collect all pages/layers we create
    let mut current_page = page1;
    let mut current_layer = layer1;

    macro_rules! get_layer {
        () => {
            doc.get_page(current_page).get_layer(current_layer)
        };
    }

    macro_rules! new_page {
        () => {{
            let (p, l) = doc.add_page(Mm(210.0), Mm(297.0), "Layer 1");
            current_page = p;
            current_layer = l;
            y = Mm(277.0);
        }};
    }

    macro_rules! check_page {
        () => {
            if y.0 < 25.0 {
                new_page!();
            }
        };
    }

    // Title
    {
        let layer = get_layer!();
        layer.use_text("GeneSight Report", 20.0, margin_left, y, &font_bold);
        y -= Mm(8.0);
        layer.use_text(
            &format!(
                "Generated: {}",
                chrono_lite_date()
            ),
            10.0,
            margin_left,
            y,
            &font,
        );
        y -= Mm(10.0);
    }

    // Summary stats
    {
        let layer = get_layer!();
        layer.use_text("Summary", 14.0, margin_left, y, &font_bold);
        y -= Mm(6.0);
        layer.use_text(
            &format!("Total Variants: {}", report.total_variants),
            10.0,
            margin_left,
            y,
            &font,
        );
        y -= line_height;
        layer.use_text(
            &format!("Annotated Variants: {}", report.annotated_variants),
            10.0,
            margin_left,
            y,
            &font,
        );
        y -= line_height;
        layer.use_text(
            &format!("Results: {}", report.results.len()),
            10.0,
            margin_left,
            y,
            &font,
        );
        y -= Mm(8.0);
    }

    // Results
    {
        let layer = get_layer!();
        layer.use_text("Results", 14.0, margin_left, y, &font_bold);
        y -= Mm(7.0);

        // Column headers
        layer.use_text("rsID", 9.0, margin_left, y, &font_bold);
        layer.use_text("Gene", 9.0, Mm(55.0), y, &font_bold);
        layer.use_text("Category", 9.0, Mm(90.0), y, &font_bold);
        layer.use_text("Tier", 9.0, Mm(120.0), y, &font_bold);
        layer.use_text("Summary", 9.0, Mm(138.0), y, &font_bold);
        y -= Mm(5.0);
    }

    for result in &report.results {
        check_page!();
        let layer = get_layer!();

        let rsid = result.variant.variant.rsid.as_deref().unwrap_or("-");
        let gene = crate::state::gene_name(result);
        let cat = crate::state::short_category(result.category);
        let tier = match result.tier {
            genesight_core::models::ConfidenceTier::Tier1Reliable => "T1",
            genesight_core::models::ConfidenceTier::Tier2Probable => "T2",
            genesight_core::models::ConfidenceTier::Tier3Speculative => "T3",
        };
        let summary = truncate_str(&result.summary, 45);

        layer.use_text(rsid, 8.0, margin_left, y, &font);
        layer.use_text(&gene, 8.0, Mm(55.0), y, &font);
        layer.use_text(cat, 8.0, Mm(90.0), y, &font);
        layer.use_text(tier, 8.0, Mm(120.0), y, &font);
        layer.use_text(&summary, 7.0, Mm(138.0), y, &font);
        y -= line_height;
    }

    // Disclaimer
    check_page!();
    y -= Mm(5.0);
    {
        let layer = get_layer!();
        layer.use_text("Disclaimer", 12.0, margin_left, y, &font_bold);
        y -= Mm(5.0);

        // Word-wrap disclaimer text
        for line in word_wrap(&report.disclaimer, 100) {
            check_page!();
            let layer = get_layer!();
            layer.use_text(&line, 8.0, margin_left, y, &font);
            y -= Mm(4.0);
        }
    }

    // Attributions
    check_page!();
    y -= Mm(4.0);
    {
        let layer = get_layer!();
        layer.use_text("Data Sources", 12.0, margin_left, y, &font_bold);
        y -= Mm(5.0);
    }
    for attr in &report.attributions {
        check_page!();
        let layer = get_layer!();
        let truncated = truncate_str(attr, 110);
        layer.use_text(&truncated, 7.0, margin_left, y, &font);
        y -= Mm(4.0);
    }

    doc.save(&mut std::io::BufWriter::new(
        std::fs::File::create(path)
            .with_context(|| format!("creating PDF at {}", path.display()))?,
    ))
    .map_err(|e| anyhow::anyhow!("PDF save error: {e}"))?;

    Ok(())
}

fn truncate_str(s: &str, max: usize) -> String {
    if s.len() <= max {
        s.to_string()
    } else {
        format!("{}...", &s[..max.saturating_sub(3)])
    }
}

fn word_wrap(text: &str, max_width: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut current_line = String::new();

    for word in text.split_whitespace() {
        if current_line.is_empty() {
            current_line = word.to_string();
        } else if current_line.len() + 1 + word.len() > max_width {
            lines.push(current_line);
            current_line = word.to_string();
        } else {
            current_line.push(' ');
            current_line.push_str(word);
        }
    }
    if !current_line.is_empty() {
        lines.push(current_line);
    }
    lines
}

fn chrono_lite_date() -> String {
    use std::time::SystemTime;
    let now = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    // Simple date calculation
    let days = now / 86400;
    let mut y = 1970i64;
    let mut remaining_days = days as i64;
    loop {
        let days_in_year = if is_leap(y) { 366 } else { 365 };
        if remaining_days < days_in_year {
            break;
        }
        remaining_days -= days_in_year;
        y += 1;
    }
    let months = [31, if is_leap(y) { 29 } else { 28 }, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    let mut m = 1;
    for &md in &months {
        if remaining_days < md {
            break;
        }
        remaining_days -= md;
        m += 1;
    }
    format!("{y}-{m:02}-{:02}", remaining_days + 1)
}

fn is_leap(y: i64) -> bool {
    (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
}
