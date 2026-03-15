use egui::{Color32, CornerRadius, FontFamily, FontId, Margin, Stroke, Style, TextStyle, Visuals};

// ── Color palette ──────────────────────────────────────────────────────────

// Backgrounds
pub const BG_PRIMARY: Color32 = Color32::from_rgb(245, 247, 250);
pub const BG_SURFACE: Color32 = Color32::WHITE;
pub const BG_SIDEBAR: Color32 = Color32::from_rgb(248, 249, 252);
pub const BG_CARD: Color32 = Color32::WHITE;

// Text
pub const TEXT_PRIMARY: Color32 = Color32::from_rgb(26, 32, 44);
pub const TEXT_SECONDARY: Color32 = Color32::from_rgb(100, 116, 139);
pub const TEXT_MUTED: Color32 = Color32::from_rgb(148, 163, 184);

// Borders
pub const BORDER: Color32 = Color32::from_rgb(226, 232, 240);
pub const BORDER_LIGHT: Color32 = Color32::from_rgb(241, 245, 249);

// Accent (indigo-based, more modern than bootstrap blue)
pub const ACCENT: Color32 = Color32::from_rgb(79, 70, 229);
pub const ACCENT_LIGHT: Color32 = Color32::from_rgb(238, 242, 255);
pub const ACCENT_HOVER: Color32 = Color32::from_rgb(67, 56, 202);

// Tiers
pub const TIER1: Color32 = Color32::from_rgb(5, 150, 105);
pub const TIER1_BG: Color32 = Color32::from_rgb(209, 250, 229);
pub const TIER2: Color32 = Color32::from_rgb(217, 119, 6);
pub const TIER2_BG: Color32 = Color32::from_rgb(254, 243, 199);
pub const TIER3: Color32 = Color32::from_rgb(100, 116, 139);
pub const TIER3_BG: Color32 = Color32::from_rgb(241, 245, 249);

// Status
pub const SUCCESS: Color32 = Color32::from_rgb(5, 150, 105);
pub const SUCCESS_BG: Color32 = Color32::from_rgb(236, 253, 245);
pub const WARNING: Color32 = Color32::from_rgb(217, 119, 6);
pub const WARNING_BG: Color32 = Color32::from_rgb(255, 251, 235);
pub const DANGER: Color32 = Color32::from_rgb(239, 68, 68);
pub const DANGER_BG: Color32 = Color32::from_rgb(254, 242, 242);

// Categories
pub const CAT_DISEASE: Color32 = Color32::from_rgb(239, 68, 68);
pub const CAT_DISEASE_BG: Color32 = Color32::from_rgb(254, 242, 242);
pub const CAT_CARRIER: Color32 = Color32::from_rgb(59, 130, 246);
pub const CAT_CARRIER_BG: Color32 = Color32::from_rgb(239, 246, 255);
pub const CAT_PHARMA: Color32 = Color32::from_rgb(139, 92, 246);
pub const CAT_PHARMA_BG: Color32 = Color32::from_rgb(245, 243, 255);
pub const CAT_PRS: Color32 = Color32::from_rgb(249, 115, 22);
pub const CAT_PRS_BG: Color32 = Color32::from_rgb(255, 247, 237);
pub const CAT_TRAIT: Color32 = Color32::from_rgb(6, 182, 212);
pub const CAT_TRAIT_BG: Color32 = Color32::from_rgb(236, 254, 255);
pub const CAT_COMPLEX: Color32 = Color32::from_rgb(100, 116, 139);
pub const CAT_COMPLEX_BG: Color32 = Color32::from_rgb(248, 250, 252);

// Selection
pub const SELECTED_ROW: Color32 = Color32::from_rgb(238, 242, 255);

// Shadow (for card borders as egui doesn't do real shadows easily)
pub const SHADOW_BORDER: Color32 = Color32::from_rgb(226, 232, 240);

/// Apply the professional medical theme to egui context.
pub fn apply_theme(ctx: &egui::Context) {
    let mut style = Style {
        visuals: Visuals {
            dark_mode: false,
            panel_fill: BG_PRIMARY,
            window_fill: BG_SURFACE,
            window_corner_radius: CornerRadius::same(12),
            window_stroke: Stroke::new(1.0, BORDER),
            window_shadow: egui::Shadow {
                offset: [0, 2],
                blur: 8,
                spread: 0,
                color: Color32::from_black_alpha(12),
            },
            extreme_bg_color: BG_SURFACE,
            faint_bg_color: BG_SIDEBAR,
            override_text_color: Some(TEXT_PRIMARY),
            hyperlink_color: ACCENT,
            ..Visuals::light()
        },
        ..Default::default()
    };

    // Non-interactive widgets (labels, frames)
    style.visuals.widgets.noninteractive.bg_fill = BG_SURFACE;
    style.visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0, BORDER_LIGHT);
    style.visuals.widgets.noninteractive.corner_radius = CornerRadius::same(6);
    style.visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, TEXT_PRIMARY);

    // Inactive widgets (buttons, checkboxes at rest)
    style.visuals.widgets.inactive.bg_fill = BG_SURFACE;
    style.visuals.widgets.inactive.bg_stroke = Stroke::new(1.0, BORDER);
    style.visuals.widgets.inactive.corner_radius = CornerRadius::same(6);
    style.visuals.widgets.inactive.fg_stroke = Stroke::new(1.0, TEXT_PRIMARY);

    // Hovered widgets
    style.visuals.widgets.hovered.bg_fill = ACCENT_LIGHT;
    style.visuals.widgets.hovered.bg_stroke = Stroke::new(1.0, ACCENT);
    style.visuals.widgets.hovered.corner_radius = CornerRadius::same(6);
    style.visuals.widgets.hovered.fg_stroke = Stroke::new(1.0, ACCENT);

    // Active/pressed widgets
    style.visuals.widgets.active.bg_fill = Color32::from_rgb(224, 231, 255);
    style.visuals.widgets.active.bg_stroke = Stroke::new(1.5, ACCENT);
    style.visuals.widgets.active.corner_radius = CornerRadius::same(6);
    style.visuals.widgets.active.fg_stroke = Stroke::new(1.5, ACCENT_HOVER);

    // Open widgets (dropdown open, etc)
    style.visuals.widgets.open.bg_fill = ACCENT_LIGHT;
    style.visuals.widgets.open.bg_stroke = Stroke::new(1.0, ACCENT);
    style.visuals.widgets.open.corner_radius = CornerRadius::same(6);

    // Selection
    style.visuals.selection.bg_fill = ACCENT_LIGHT;
    style.visuals.selection.stroke = Stroke::new(1.0, ACCENT);

    // Scrollbar
    style.visuals.handle_shape = egui::style::HandleShape::Circle;

    // Spacing
    style.spacing.item_spacing = egui::vec2(8.0, 6.0);
    style.spacing.window_margin = Margin::same(16);
    style.spacing.button_padding = egui::vec2(16.0, 8.0);
    style.spacing.indent = 20.0;
    style.spacing.scroll = egui::style::ScrollStyle {
        bar_width: 6.0,
        bar_inner_margin: 2.0,
        bar_outer_margin: 2.0,
        ..Default::default()
    };

    // Text styles
    style.text_styles = [
        (
            TextStyle::Heading,
            FontId::new(22.0, FontFamily::Proportional),
        ),
        (TextStyle::Body, FontId::new(14.0, FontFamily::Proportional)),
        (
            TextStyle::Monospace,
            FontId::new(13.0, FontFamily::Monospace),
        ),
        (
            TextStyle::Button,
            FontId::new(14.0, FontFamily::Proportional),
        ),
        (
            TextStyle::Small,
            FontId::new(12.0, FontFamily::Proportional),
        ),
    ]
    .into();

    ctx.set_style(style);
}

// ── Helpers ────────────────────────────────────────────────────────────────

/// Color for a confidence tier.
pub fn tier_color(tier: genesight_core::models::ConfidenceTier) -> Color32 {
    use genesight_core::models::ConfidenceTier;
    match tier {
        ConfidenceTier::Tier1Reliable => TIER1,
        ConfidenceTier::Tier2Probable => TIER2,
        ConfidenceTier::Tier3Speculative => TIER3,
    }
}

/// Background color for a confidence tier badge.
pub fn tier_bg(tier: genesight_core::models::ConfidenceTier) -> Color32 {
    use genesight_core::models::ConfidenceTier;
    match tier {
        ConfidenceTier::Tier1Reliable => TIER1_BG,
        ConfidenceTier::Tier2Probable => TIER2_BG,
        ConfidenceTier::Tier3Speculative => TIER3_BG,
    }
}

/// Short label for a tier.
pub fn tier_label(tier: genesight_core::models::ConfidenceTier) -> &'static str {
    use genesight_core::models::ConfidenceTier;
    match tier {
        ConfidenceTier::Tier1Reliable => "Tier 1: Reliable",
        ConfidenceTier::Tier2Probable => "Tier 2: Probable",
        ConfidenceTier::Tier3Speculative => "Tier 3: Speculative",
    }
}

/// Color for a result category.
pub fn category_color(cat: genesight_core::models::ResultCategory) -> Color32 {
    use genesight_core::models::ResultCategory;
    match cat {
        ResultCategory::MonogenicDisease => CAT_DISEASE,
        ResultCategory::CarrierStatus => CAT_CARRIER,
        ResultCategory::Pharmacogenomics => CAT_PHARMA,
        ResultCategory::GwasAssociation => CAT_PRS,
        ResultCategory::PhysicalTrait => CAT_TRAIT,
        ResultCategory::ComplexTrait
        | ResultCategory::Ancestry
        | ResultCategory::ClinVarConflicting => CAT_COMPLEX,
    }
}

/// Background color for a category badge.
pub fn category_bg(cat: genesight_core::models::ResultCategory) -> Color32 {
    use genesight_core::models::ResultCategory;
    match cat {
        ResultCategory::MonogenicDisease => CAT_DISEASE_BG,
        ResultCategory::CarrierStatus => CAT_CARRIER_BG,
        ResultCategory::Pharmacogenomics => CAT_PHARMA_BG,
        ResultCategory::GwasAssociation => CAT_PRS_BG,
        ResultCategory::PhysicalTrait => CAT_TRAIT_BG,
        ResultCategory::ComplexTrait
        | ResultCategory::Ancestry
        | ResultCategory::ClinVarConflicting => CAT_COMPLEX_BG,
    }
}

/// Draw a card frame with subtle elevation.
pub fn card_frame() -> egui::Frame {
    egui::Frame::NONE
        .fill(BG_CARD)
        .stroke(Stroke::new(1.0, SHADOW_BORDER))
        .corner_radius(CornerRadius::same(12))
        .inner_margin(Margin::same(24))
        .shadow(egui::Shadow {
            offset: [0, 1],
            blur: 4,
            spread: 0,
            color: Color32::from_black_alpha(8),
        })
}

/// Draw a subtle section frame.
pub fn section_frame() -> egui::Frame {
    egui::Frame::NONE
        .fill(BG_SIDEBAR)
        .corner_radius(CornerRadius::same(8))
        .inner_margin(Margin::same(12))
}
