#!/usr/bin/env python3
"""
Fig 1C - Seven-Dimension Evaluation Framework + Reproducibility

Usage (in Docker):
    docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
        python3 /benchmark/scripts/draw_fig1c_evaluation_framework.py
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from pathlib import Path
import numpy as np
import os, sys

print(f"Python: {sys.version}")
print(f"matplotlib: {matplotlib.__version__}")

# ── Paths ──
if Path('/benchmark/icons').exists():
    ICONS_DIR = Path('/benchmark/icons')
    FIGURES_DIR = Path('/benchmark/figures/main')
else:
    ICONS_DIR = Path('benchmark/icons')
    FIGURES_DIR = Path('benchmark/figures/main')

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
print(f"ICONS_DIR: {ICONS_DIR} (exists={ICONS_DIR.exists()})")
print(f"Icons found: {sorted(p.name for p in ICONS_DIR.glob('*.png'))}")

# ── Colors ──
CC_PINK = '#E91E63'
CC_BG = '#FCE4EC'
TXT1 = '#212121'
TXT2 = '#616161'
CARD_BG = '#F5F5F5'
CARD_BORDER = '#BDBDBD'
BADGE_BG = '#FAFAFA'
BADGE_BD = '#9E9E9E'
ACCENT_GOLD = '#FFC107'
DIM_AXIS = '#E0E0E0'
AXIS_COLOR = '#9E9E9E'


def load_icon(filename, target_h=20):
    """Load icon PNG; return None if missing."""
    path = ICONS_DIR / filename
    if not path.exists():
        print(f"  WARN: icon missing: {path}")
        return None
    try:
        img = mpimg.imread(str(path))
        zoom = target_h / img.shape[0]
        return OffsetImage(img, zoom=zoom)
    except Exception as e:
        print(f"  ERROR loading {path}: {e}")
        return None



# ── 7 Evaluation Dimensions ──
DIMENSIONS = [
    {'label': 'Conversion\nRobustness',    'icon': 'icon_shield.png',    'emoji': '🛡️', 'angle_idx': 0},
    {'label': 'Runtime\nPerformance',      'icon': 'icon_clock.png',     'emoji': '⏱️', 'angle_idx': 1},
    {'label': 'Memory\nEfficiency',        'icon': 'icon_chip.png',      'emoji': '💾', 'angle_idx': 2},
    {'label': 'Expression\nFidelity',      'icon': 'icon_checkmark.png', 'emoji': '✓',  'angle_idx': 3},
    {'label': 'Type\nIntegrity',           'icon': 'icon_lock.png',      'emoji': '🔒', 'angle_idx': 4},
    {'label': 'V5\nCompatibility',         'icon': 'icon_layers.png',    'emoji': '📚', 'angle_idx': 5},
    {'label': 'Large-scale\nScalability',  'icon': 'icon_scaleup.png',   'emoji': '📈', 'angle_idx': 6},
]

# ── Reproducibility items ──
REPRO_ITEMS = [
    {'icon': 'icon_docker.png',    'emoji': '🐳', 'text': 'Version-locked environment'},
    {'icon': 'icon_replicate.png', 'emoji': '3×', 'text': '3× replicate runs per test'},
    {'icon': 'icon_github.png',    'emoji': '🔗', 'text': 'All code, data, notebooks public'},
    {'icon': 'icon_grid.png',      'emoji': '📊', 'text': '65 datasets × 5 tools × 7 dim\n= 2,275 test points'},
]


def draw():
    fig_w, fig_h = 13.0, 7.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)
    ax.axis('off')
    fig.patch.set_facecolor('white')

    # ── Title ──
    ax.text(fig_w / 2, fig_h - 0.3,
            'Seven-Dimension Evaluation Framework',
            fontsize=15, fontweight='bold', ha='center', va='center',
            color=TXT1, fontfamily='sans-serif')

    # ════════════════════════════════════════════
    # LEFT: Radar chart skeleton (7 axes)
    # ════════════════════════════════════════════
    cx, cy = 3.8, 3.4   # center of radar
    R_outer = 2.4        # outer radius for axis lines
    R_label = 2.95       # label position
    R_icon = 2.55        # icon position
    N = 7
    angles = [np.pi/2 + 2*np.pi*i/N for i in range(N)]  # start from top

    # Draw concentric rings (3 levels)
    for r_frac in [0.33, 0.66, 1.0]:
        r = R_outer * r_frac
        ring_angles = np.linspace(0, 2*np.pi, 200)
        xs = cx + r * np.cos(ring_angles + np.pi/2)
        ys = cy + r * np.sin(ring_angles + np.pi/2)
        ax.plot(xs, ys, color=DIM_AXIS, linewidth=0.6, alpha=0.5)

    # Draw axis lines and labels
    for dim in DIMENSIONS:
        i = dim['angle_idx']
        a = angles[i]
        dx, dy = np.cos(a), np.sin(a)

        # Axis line from center to outer
        ax.plot([cx, cx + R_outer*dx], [cy, cy + R_outer*dy],
                color=AXIS_COLOR, linewidth=1.0, alpha=0.7)

        # Axis endpoint dot
        ax.plot(cx + R_outer*dx, cy + R_outer*dy,
                'o', color=AXIS_COLOR, markersize=4, zorder=5)

        # Icon at axis tip
        im = load_icon(dim['icon'], target_h=22)
        if im is not None:
            ab = AnnotationBbox(im, (cx + R_icon*dx, cy + R_icon*dy),
                                frameon=False, zorder=6)
            ax.add_artist(ab)
        else:
            # Fallback: emoji/text
            ax.text(cx + R_icon*dx, cy + R_icon*dy, dim['emoji'],
                    fontsize=14, ha='center', va='center', zorder=6)

        # Label
        lx = cx + R_label * dx
        ly = cy + R_label * dy
        ha = 'center'
        if dx > 0.3:
            ha = 'left'
        elif dx < -0.3:
            ha = 'right'

        ax.text(lx, ly, dim['label'],
                fontsize=8.5, ha=ha, va='center',
                color=TXT1, fontfamily='sans-serif',
                linespacing=1.1)

    # Draw a sample filled polygon (CrossCell = perfect heptagon)
    poly_r = R_outer * 0.85
    poly_xs = [cx + poly_r * np.cos(a) for a in angles]
    poly_ys = [cy + poly_r * np.sin(a) for a in angles]
    poly_xs.append(poly_xs[0])
    poly_ys.append(poly_ys[0])
    ax.fill(poly_xs, poly_ys, color=CC_PINK, alpha=0.12, zorder=2)
    ax.plot(poly_xs, poly_ys, color=CC_PINK, linewidth=1.5, alpha=0.6, zorder=3)

    # Center label
    ax.text(cx, cy - 0.15, 'CrossCell', fontsize=9, fontweight='bold',
            ha='center', va='center', color=CC_PINK, alpha=0.8)

    # ════════════════════════════════════════════
    # RIGHT: Reproducibility box
    # ════════════════════════════════════════════
    box_x = 8.2
    box_y = 1.2
    box_w = 3.5
    box_h = 4.4

    # Box background
    repro_box = mpatches.FancyBboxPatch(
        (box_x, box_y), box_w, box_h,
        boxstyle='round,pad=0.15',
        facecolor='#FAFAFA', edgecolor=CARD_BORDER, linewidth=1.2)
    ax.add_patch(repro_box)

    # Box title
    ax.text(box_x + box_w/2, box_y + box_h - 0.3,
            'Reproducibility', fontsize=11, fontweight='bold',
            ha='center', va='center', color=TXT1, fontfamily='sans-serif')

    # Divider line under title
    ax.plot([box_x + 0.3, box_x + box_w - 0.3],
            [box_y + box_h - 0.55, box_y + box_h - 0.55],
            color=CARD_BORDER, linewidth=0.8)

    # Reproducibility items
    item_start_y = box_y + box_h - 0.85
    item_spacing = 0.9

    for idx, item in enumerate(REPRO_ITEMS):
        iy = item_start_y - idx * item_spacing
        icon_x = box_x + 0.45
        text_x = box_x + 0.9

        # Icon
        im = load_icon(item['icon'], target_h=18)
        if im is not None:
            ab = AnnotationBbox(im, (icon_x, iy),
                                frameon=False, zorder=6)
            ax.add_artist(ab)
        else:
            ax.text(icon_x, iy, item['emoji'],
                    fontsize=13, ha='center', va='center')

        # Text
        ax.text(text_x, iy, item['text'],
                fontsize=8.5, ha='left', va='center',
                color=TXT2, fontfamily='sans-serif',
                linespacing=1.2)

    # ── Save ──
    for fmt in ['pdf', 'png']:
        out = FIGURES_DIR / f'fig1c_evaluation_framework.{fmt}'
        print(f"Saving {out} ...")
        fig.savefig(str(out), dpi=300, bbox_inches='tight', facecolor='white')
        if out.exists():
            print(f"  OK: {out} ({out.stat().st_size} bytes)")
        else:
            print(f"  FAILED: {out} not created")
    plt.close(fig)
    print('Done!')


if __name__ == '__main__':
    draw()
