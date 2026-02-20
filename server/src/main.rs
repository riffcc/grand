use axum::{
    routing::get,
    Router,
};
use tower_http::services::ServeDir;

#[tokio::main]
async fn main() {
    let app = Router::new()
        .route("/", get(serve_index))
        .route("/api/health", get(health))
        .nest_service("/static", ServeDir::new("../frontend"));

    let listener = tokio::net::TcpListener::bind("0.0.0.0:9000").await.unwrap();
    println!("🌌 GUTOE Server running at http://localhost:9000");
    println!("   Void awaits...");

    axum::serve(listener, app).await.unwrap();
}

async fn serve_index() -> axum::response::Html<String> {
    let html = std::fs::read_to_string("/mnt/riffcastle/castle/garage/grand-2026/frontend/index.html").unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        "<html><body><h1>Void</h1></body></html>".to_string()
    });
    axum::response::Html(html)
}

async fn health() -> &'static str {
    r#"{"status": "void", "coherence": 1.0, "veracity": 1.0}"#
}
