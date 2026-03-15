mod app;
mod export;
mod state;
mod theme;
mod views;

fn main() -> eframe::Result<()> {
    tracing_subscriber::fmt::init();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 860.0])
            .with_min_inner_size([900.0, 560.0])
            .with_title("GeneSight"),
        ..Default::default()
    };

    eframe::run_native(
        "GeneSight",
        options,
        Box::new(|cc| Ok(Box::new(app::GeneSightApp::new(cc)))),
    )
}
