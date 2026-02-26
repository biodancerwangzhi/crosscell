#!/usr/bin/env python3
"""
Figure 1 Composite: A (dataset landscape) + B (tool landscape) + C (evaluation framework)

Layout: Left column = A (top) + B (bottom), Right column = C (full height)

Usage (in Docker):
    docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
        python3 /benchmark/scripts/draw_fig1_composite.py
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from pathlib import Path
import numpy as np
import json, os, sys

print(f"Python: {sys.version}")
print(f"matplotlib: {matplotlib.__version__}")

# ── Paths ──
if Path('/benchmark/icons').exists():
    ICONS_DIR = Path('/benchmark/icons')
    FIGURES_DIR = Path('/benchmark/figures/main')
    CONFIG_DIR = Path('/benchmark/config')
else:
    ICONS_DIR = Path('benchmark/icons')
    FIGURES_DIR = Path('benchmark/figures/main')
    CONFIG_DIR = Path('benchmark/config')

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
print(f"ICONS_DIR: {ICONS_DIR}")
print(f"Icons: {sorted(p.name for p in ICONS_DIR.glob('*.png'))}")

# ── Load test cases for panel A ──
TESTCASES_FILE = CONFIG_DIR / 'benchmark_testcases.json'
with open(TESTCASES_FILE) as f:
    testcases = json.load(f)
print(f"Loaded {testcases['summary']['total_test_cases']} test cases")


# ══════════════════════════════════════════════════════════════
# Shared constants
# ══════════════════════════════════════════════════════════════
CC_PINK = '#E91E63'
CC_BG   = '#FCE4EC'
TXT1    = '#212121'
TXT2    = '#616161'
CARD_BG = '#F5F5F5'
CARD_BD = '#BDBDBD'
BADGE_BG = '#FAFAFA'
BADGE_BD = '#9E9E9E'
DIM_AXIS = '#E0E0E0'
AXIS_CLR = '#9E9E9E'

SCALE_COLORS = {'small': '#4CAF50', 'medium': '#2196F3', 'large': '#E91E63'}


def load_icon(filename, target_h=20):
    path = ICONS_DIR / filename
    if not path.exists():
        return None
    try:
        img = mpimg.imread(str(path))
        return OffsetImage(img, zoom=target_h / img.shape[0])
    except:
        return None


# ══════════════════════════════════════════════════════════════
# Prepare dataset table for Panel A
# ══════════════════════════════════════════════════════════════
rows = []
for tc in testcases['test_cases'].get('rds_to_h5ad', []):
    sv = tc.get('seurat_version', 'V5')
    rows.append({'test_id': tc['test_id'], 'source': f'SeuratData {sv}', 'format': 'RDS',
                 'n_cells': tc['n_cells'], 'n_genes': tc['n_genes'],
                 'spatial': tc.get('spatial', False)})
for tc in testcases['test_cases'].get('h5ad_to_rds', []):
    src = tc.get('source', 'unknown')
    source_map = {'scanpy': 'scanpy', 'scvelo': 'scvelo', 'cellxgene': 'CELLxGENE', 'squidpy': 'squidpy'}
    rows.append({'test_id': tc['test_id'], 'source': source_map.get(src, src), 'format': 'H5AD',
                 'n_cells': tc['n_cells'], 'n_genes': tc['n_genes'],
                 'spatial': tc.get('spatial', False)})

cellxgene_large = [
    ('cellxgene_skin_bcc_10k', 9841, 26886), ('cellxgene_tabula_liver_22k', 22214, 60606),
    ('cellxgene_kidney_atacseq_37k', 37747, 19276), ('cellxgene_brain_multiome_102k', 101924, 35451),
    ('cellxgene_pancreas_122k', 121916, 32356), ('cellxgene_brain_dlpfc_172k', 172120, 37490),
    ('cellxgene_gut_428k', 428469, 32383), ('cellxgene_heart_486k', 486134, 32383),
    ('cellxgene_hlca_core_585k', 584944, 27402), ('cellxgene_combat_pbmc_836k', 836148, 36306),
    ('cellxgene_eqtl_autoimmune_1.2M', 1248980, 35528),
]
existing_ids = {r['test_id'] for r in rows}
for tid, nc, ng in cellxgene_large:
    if tid not in existing_ids:
        rows.append({'test_id': tid, 'source': 'CELLxGENE', 'format': 'H5AD',
                     'n_cells': nc, 'n_genes': ng, 'spatial': False})

for r in rows:
    r['scale'] = 'small' if r['n_cells'] < 10000 else ('medium' if r['n_cells'] < 100000 else 'large')

print(f"Dataset rows: {len(rows)}")


# ══════════════════════════════════════════════════════════════
# DRAW
# ══════════════════════════════════════════════════════════════
def draw():
    # ── Global figure: left col (A top, B bottom), right col (C full) ──
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 2, figure=fig,
                  width_ratios=[1, 1.15],
                  height_ratios=[1.1, 1],
                  hspace=0.12, wspace=0.08)

    # Consistent font sizes
    TITLE_FS = 11
    LABEL_FS = 9
    TICK_FS = 8
    SMALL_FS = 7.5

    # ════════════════════════════════════════
    # Panel A: Dataset Landscape (top-left)
    # ════════════════════════════════════════
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_title('A', fontsize=14, fontweight='bold', loc='left', pad=6)

    source_order = ['SeuratData V4', 'SeuratData V5', 'scanpy', 'scvelo', 'CELLxGENE', 'squidpy']
    source_x = {s: i for i, s in enumerate(source_order)}
    np.random.seed(42)

    for r in rows:
        x_base = source_x.get(r['source'], 0)
        x = x_base + np.random.uniform(-0.25, 0.25)
        size = np.clip(r['n_genes'] / 500, 15, 200)
        color = SCALE_COLORS.get(r['scale'], '#999')
        marker = 's' if r['spatial'] else 'o'
        ax_a.scatter(x, r['n_cells'], s=size, c=color, marker=marker,
                     alpha=0.7, edgecolor='white', linewidth=0.5, zorder=3)

    ax_a.axvline(x=1.6, color='#BDBDBD', linestyle=':', linewidth=1, alpha=0.5)
    ax_a.set_yscale('log')
    ax_a.set_xticks(range(len(source_order)))
    ax_a.set_xticklabels(source_order, rotation=0, ha='center', fontsize=TICK_FS)
    ax_a.set_ylabel('Number of Cells (log scale)', fontsize=LABEL_FS)
    ax_a.tick_params(axis='y', labelsize=TICK_FS)
    ax_a.text(0.5, 1.08, 'Dataset Landscape: 65 Benchmark Datasets',
              transform=ax_a.transAxes, fontsize=TITLE_FS, fontweight='bold',
              ha='center', va='center', color=TXT1, fontfamily='sans-serif')

    # Legends
    tier_patches = [mpatches.Patch(color=SCALE_COLORS[s], label=l)
                    for s, l in [('small', 'Small (<10k)'), ('medium', 'Medium (10k-100k)'), ('large', 'Large (≥100k)')]]
    leg1 = ax_a.legend(handles=tier_patches, loc='upper left', fontsize=SMALL_FS,
                       title='Cell Count', title_fontsize=TICK_FS, bbox_to_anchor=(0.0, 1.0))
    ax_a.add_artist(leg1)

    spatial_h = [Line2D([0], [0], marker='o', color='gray', label='Transcriptomic', markersize=6, linestyle=''),
                 Line2D([0], [0], marker='s', color='gray', label='Spatial', markersize=6, linestyle='')]
    leg2 = ax_a.legend(handles=spatial_h, loc='upper left', fontsize=SMALL_FS,
                       title='Data Type', title_fontsize=TICK_FS, bbox_to_anchor=(0.28, 1.0))
    ax_a.add_artist(leg2)

    # Gene Count legend (bubble size)
    size_handles = [
        Line2D([0], [0], marker='o', color='#BDBDBD', markersize=np.sqrt(np.clip(500/500, 15, 200)),
               label='<1k genes', linestyle=''),
        Line2D([0], [0], marker='o', color='#BDBDBD', markersize=np.sqrt(np.clip(10000/500, 15, 200)),
               label='1k–20k genes', linestyle=''),
        Line2D([0], [0], marker='o', color='#BDBDBD', markersize=np.sqrt(np.clip(40000/500, 15, 200)),
               label='>20k genes', linestyle=''),
    ]
    leg3 = ax_a.legend(handles=size_handles, loc='upper left', fontsize=SMALL_FS,
                       title='Gene Count', title_fontsize=TICK_FS, bbox_to_anchor=(0.0, 0.85))
    ax_a.add_artist(leg3)

    n_rds = sum(1 for r in rows if r['format'] == 'RDS')
    n_h5ad = sum(1 for r in rows if r['format'] == 'H5AD')
    n_cx = sum(1 for r in rows if r['source'] == 'CELLxGENE')
    ax_a.text(0.28, 0.85, f'{len(rows)} datasets\n{n_rds} RDS + {n_h5ad} H5AD\nCELLxGENE: {n_cx}',
              transform=ax_a.transAxes, ha='left', va='top', fontsize=SMALL_FS,
              bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#ccc', linewidth=0.8))

    ax_a.grid(axis='y', alpha=0.3)
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)


    # ════════════════════════════════════════
    # Panel B: Tool Landscape (bottom-left)
    # ════════════════════════════════════════
    ax_b = fig.add_subplot(gs[1, 0])
    ax_b.set_xlim(0, 10)
    ax_b.set_ylim(0, 5.25)
    ax_b.axis('off')
    ax_b.set_title('B', fontsize=14, fontweight='bold', loc='left', pad=6)

    TOOLS_1B = [
        {'name': 'Zellkonverter',   'desc': 'R→Python bridge (basilisk) — needs R + Python',
         'icons': ['r_logo.png', 'python_logo.png'], 'direction': '↔', 'deps': 'Runtime Deps: R, Python', 'hl': False},
        {'name': 'anndataR',        'desc': 'Native R (hdf5r) — needs R only',
         'icons': ['r_logo.png', 'file_logo.png'],   'direction': '←', 'deps': 'Runtime Deps: R', 'hl': False},
        {'name': 'convert2anndata', 'desc': 'R→Python bridge (reticulate) — needs R + Python',
         'icons': ['r_logo.png', 'python_logo.png'], 'direction': '→', 'deps': 'Runtime Deps: R, Python', 'hl': False},
        {'name': 'easySCF',         'desc': 'Custom HDF5 intermediate — needs R + Python',
         'icons': ['hdf5_logo.png', 'file_logo.png'],'direction': '↔', 'deps': 'Runtime Deps: R, Python', 'hl': False},
        {'name': 'CrossCell',       'desc': 'Native Rust (no runtime) — standalone binary',
         'icons': ['rust_logo.png', 'file_logo.png'],'direction': '↔', 'deps': 'Runtime Deps: None', 'hl': True},
    ]

    W_b = 10.0
    ROW_H = 0.85
    N_tools = len(TOOLS_1B)
    TITLE_H_B = 0.45
    y_start_b = 5.25 - TITLE_H_B - 0.1

    ax_b.text(W_b / 2, 5.25 - 0.15, 'Tool Comparison: Technical Approaches',
              fontsize=TITLE_FS, fontweight='bold', ha='center', va='center',
              color=TXT1, fontfamily='sans-serif')

    outer_b = mpatches.FancyBboxPatch((0.15, 0.1), W_b - 0.3, 5.25 - TITLE_H_B - 0.05,
                                       boxstyle='round,pad=0.08', facecolor='white',
                                       edgecolor=CARD_BD, linewidth=1.0)
    ax_b.add_patch(outer_b)

    CL, CW = 0.3, W_b - 0.6
    for i, t in enumerate(TOOLS_1B):
        yt = y_start_b - i * ROW_H
        yc = yt - ROW_H / 2
        cy = yt - ROW_H + 0.08

        bg, bc, bw = (CC_BG, CC_PINK, 1.8) if t['hl'] else (CARD_BG, CARD_BD, 0.8)
        card = mpatches.FancyBboxPatch((CL, cy), CW, ROW_H - 0.14,
                                        boxstyle='round,pad=0.06', facecolor=bg,
                                        edgecolor=bc, linewidth=bw)
        ax_b.add_patch(card)

        nx = CL + 0.2
        ax_b.text(nx, yc + 0.12, t['name'], fontsize=LABEL_FS, fontweight='bold',
                  ha='left', va='center', color=TXT1, fontfamily='sans-serif')

        for j, ic_file in enumerate(t['icons']):
            im = load_icon(ic_file, target_h=12)
            if im is not None:
                ab = AnnotationBbox(im, (nx + j * 0.35, yc - 0.1),
                                    frameon=False, box_alignment=(0, 0.5))
                ax_b.add_artist(ab)

        dx = nx + len(t['icons']) * 0.35 + 0.15
        ax_b.text(dx, yc - 0.1, t['desc'], fontsize=SMALL_FS, ha='left', va='center',
                  color=TXT2, fontfamily='sans-serif')

        # Direction badge (shifted left)
        ax_r = CL + CW - 3.6
        db = mpatches.FancyBboxPatch((ax_r - 0.15, yc - 0.12), 0.3, 0.24,
                                      boxstyle='round,pad=0.04', facecolor=BADGE_BG,
                                      edgecolor=BADGE_BD, linewidth=0.6)
        ax_b.add_patch(db)
        ax_b.text(ax_r, yc, t['direction'], fontsize=12, ha='center', va='center',
                  color=TXT1, fontfamily='sans-serif')

        # Deps badge (shifted left, wider)
        rx = CL + CW - 1.8
        dtxt = t['deps']
        bw2 = max(len(dtxt) * 0.08, 2.2)
        rb = mpatches.FancyBboxPatch((rx - bw2/2, yc - 0.11), bw2, 0.22,
                                      boxstyle='round,pad=0.04', facecolor=BADGE_BG,
                                      edgecolor=BADGE_BD, linewidth=0.6)
        ax_b.add_patch(rb)
        ax_b.text(rx, yc, dtxt, fontsize=SMALL_FS, ha='center', va='center',
                  color=TXT2, fontfamily='sans-serif')


    # ════════════════════════════════════════
    # Panel C: Evaluation Framework (right, full height)
    # ════════════════════════════════════════
    ax_c = fig.add_subplot(gs[:, 1])
    ax_c.set_xlim(0, 10)
    ax_c.set_ylim(0, 10)
    ax_c.axis('off')
    ax_c.set_title('C', fontsize=14, fontweight='bold', loc='left', pad=6)

    ax_c.text(5.0, 9.7, 'Seven-Dimension Evaluation Framework',
              fontsize=TITLE_FS, fontweight='bold', ha='center', va='center',
              color=TXT1, fontfamily='sans-serif')

    # ── Radar chart (compact height) ──
    cx_c, cy_c = 5.0, 6.2
    R_outer = 2.2
    R_label = 2.8
    R_icon = 2.45
    N_dim = 7
    angles_c = [np.pi/2 + 2*np.pi*i/N_dim for i in range(N_dim)]

    DIMENSIONS = [
        {'label': 'Conversion\nRobustness',   'icon': 'icon_shield.png',    'emoji': '🛡️', 'idx': 0},
        {'label': 'Runtime\nPerformance',     'icon': 'icon_clock.png',     'emoji': '⏱️', 'idx': 1},
        {'label': 'Memory\nEfficiency',       'icon': 'icon_chip.png',      'emoji': '💾', 'idx': 2},
        {'label': 'Expression\nFidelity',     'icon': 'icon_checkmark.png', 'emoji': '✓',  'idx': 3},
        {'label': 'Type\nIntegrity',          'icon': 'icon_lock.png',      'emoji': '🔒', 'idx': 4},
        {'label': 'V5\nCompatibility',        'icon': 'icon_layers.png',    'emoji': '📚', 'idx': 5},
        {'label': 'Large-scale\nScalability', 'icon': 'icon_scaleup.png',   'emoji': '📈', 'idx': 6},
    ]

    # Concentric rings
    for r_frac in [0.33, 0.66, 1.0]:
        r = R_outer * r_frac
        ring_a = np.linspace(0, 2*np.pi, 200)
        ax_c.plot(cx_c + r * np.cos(ring_a + np.pi/2),
                  cy_c + r * np.sin(ring_a + np.pi/2),
                  color=DIM_AXIS, linewidth=0.6, alpha=0.5)

    # Axes, icons, labels
    for dim in DIMENSIONS:
        a = angles_c[dim['idx']]
        ddx, ddy = np.cos(a), np.sin(a)

        ax_c.plot([cx_c, cx_c + R_outer*ddx], [cy_c, cy_c + R_outer*ddy],
                  color=AXIS_CLR, linewidth=0.8, alpha=0.7)
        ax_c.plot(cx_c + R_outer*ddx, cy_c + R_outer*ddy,
                  'o', color=AXIS_CLR, markersize=3.5, zorder=5)

        im = load_icon(dim['icon'], target_h=18)
        if im is not None:
            ab = AnnotationBbox(im, (cx_c + R_icon*ddx, cy_c + R_icon*ddy),
                                frameon=False, zorder=6)
            ax_c.add_artist(ab)
        else:
            ax_c.text(cx_c + R_icon*ddx, cy_c + R_icon*ddy, dim['emoji'],
                      fontsize=12, ha='center', va='center', zorder=6)

        lx = cx_c + R_label * ddx
        ly = cy_c + R_label * ddy
        ha_l = 'center'
        if ddx > 0.3: ha_l = 'left'
        elif ddx < -0.3: ha_l = 'right'
        ax_c.text(lx, ly, dim['label'], fontsize=TICK_FS, ha=ha_l, va='center',
                  color=TXT1, fontfamily='sans-serif', linespacing=1.1)

    # CrossCell polygon
    poly_r = R_outer * 0.85
    pxs = [cx_c + poly_r * np.cos(a) for a in angles_c]
    pys = [cy_c + poly_r * np.sin(a) for a in angles_c]
    pxs.append(pxs[0]); pys.append(pys[0])
    ax_c.fill(pxs, pys, color=CC_PINK, alpha=0.12, zorder=2)
    ax_c.plot(pxs, pys, color=CC_PINK, linewidth=1.5, alpha=0.6, zorder=3)
    ax_c.text(cx_c, cy_c - 0.15, 'CrossCell', fontsize=LABEL_FS, fontweight='bold',
              ha='center', va='center', color=CC_PINK, alpha=0.8)

    # ── Reproducibility box (below radar, vertical 1-per-row) ──
    bx, by, bw, bh = 1.0, 0.2, 8.0, 2.2
    repro = mpatches.FancyBboxPatch((bx, by), bw, bh, boxstyle='round,pad=0.12',
                                     facecolor='#FAFAFA', edgecolor=CARD_BD, linewidth=1.0)
    ax_c.add_patch(repro)
    ax_c.text(bx + bw/2, by + bh - 0.2, 'Reproducibility',
              fontsize=LABEL_FS, fontweight='bold', ha='center', va='center',
              color=TXT1, fontfamily='sans-serif')
    ax_c.plot([bx + 0.3, bx + bw - 0.3], [by + bh - 0.38, by + bh - 0.38],
              color=CARD_BD, linewidth=0.6)

    REPRO = [
        {'icon': 'icon_docker.png',    'emoji': '🐳', 'text': 'Version-locked environment'},
        {'icon': 'icon_replicate.png', 'emoji': '3×', 'text': '3× replicate runs per test'},
        {'icon': 'icon_github.png',    'emoji': '🔗', 'text': 'All code, data, notebooks public'},
        {'icon': 'icon_grid.png',      'emoji': '📊', 'text': '65 datasets × 5 tools × 7 dim = 2,275 test points'},
    ]

    # Vertical layout: one item per row, tighter spacing
    item_spacing = 0.42
    item_start_y = by + bh - 0.58
    for idx, item in enumerate(REPRO):
        iy = item_start_y - idx * item_spacing
        icon_x = bx + 0.5

        im = load_icon(item['icon'], target_h=14)
        if im is not None:
            ab = AnnotationBbox(im, (icon_x, iy), frameon=False, zorder=6)
            ax_c.add_artist(ab)
        else:
            ax_c.text(icon_x, iy, item['emoji'], fontsize=10, ha='center', va='center')

        ax_c.text(icon_x + 0.35, iy, item['text'], fontsize=SMALL_FS, ha='left', va='center',
                  color=TXT2, fontfamily='sans-serif')

    # ── Save ──
    for fmt in ['pdf', 'png']:
        out = FIGURES_DIR / f'fig1_composite.{fmt}'
        print(f"Saving {out} ...")
        fig.savefig(str(out), dpi=300, bbox_inches='tight', facecolor='white')
        if out.exists():
            print(f"  OK: {out} ({out.stat().st_size} bytes)")
    plt.close(fig)
    print('Done!')


if __name__ == '__main__':
    draw()
