# Finding 026: Kerr Metrics Baseline (m87star, 640x360)

Date: 2026-02-24
Command:
`cargo run -p gutoe-gpu --bin bh_render -- kerr_metrics m87star 640x360`

Output CSV: `/tmp/bh_renders/kerr_metrics_m87star.csv`

Spin sweep summary:
- a*=0.00: cx=319.60, cy=185.32, D_sh=186.70, ring_w=1.00, A_LR=-0.0119, CP=-4.33°
- a*=0.30: cx=320.03, cy=185.32, D_sh=187.07, ring_w=1.00, A_LR=-0.0115, CP=-3.36°
- a*=0.60: cx=320.46, cy=185.28, D_sh=187.43, ring_w=1.00, A_LR=-0.0110, CP=-0.98°
- a*=0.90: cx=320.90, cy=185.18, D_sh=188.29, ring_w=1.00, A_LR=-0.0096, CP=-7.67°

Interpretation:
- Increasing spin shifts centroid rightward (cx increases), as expected from frame dragging.
- Shadow diameter drift is small but monotonic in this setup.
- This dataset serves as a baseline target for live `bh_viewer` Kerr parity sweeps.
