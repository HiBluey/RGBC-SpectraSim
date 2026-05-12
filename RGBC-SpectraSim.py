import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons

# ==========================================
# 0. 字体支持
# ==========================================
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'PingFang SC', 'Heiti TC', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# ==========================================
# 1. 配置物理参数
# ==========================================
WAVELENGTHS = np.arange(380, 781, 1)

PEAKS = {'R': 630, 'G': 530, 'C': 495, 'B': 460}
SIGMAS = {'R': 8, 'G': 10, 'C': 8, 'B': 8}

XY = {
    'R': np.array([0.68, 0.32]),
    'G': np.array([0.21, 0.71]),
    'C': np.array([0.12, 0.45]),
    'B': np.array([0.15, 0.06])
}

CIE_LOCUS = np.array([
    [0.1741, 0.0050], [0.1730, 0.0048], [0.1708, 0.0055], [0.1666, 0.0083],
    [0.1584, 0.0163], [0.1450, 0.0336], [0.1259, 0.0673], [0.1011, 0.1232],
    [0.0718, 0.2057], [0.0401, 0.3168], [0.0136, 0.4497], [0.0039, 0.5878],
    [0.0069, 0.7046], [0.0227, 0.7766], [0.0543, 0.8143], [0.1096, 0.8256],
    [0.1804, 0.8023], [0.2612, 0.7320], [0.3475, 0.6436], [0.4316, 0.5574],
    [0.5113, 0.4786], [0.5843, 0.4109], [0.6457, 0.3541], [0.6908, 0.3090],
    [0.7212, 0.2787], [0.7347, 0.2653]
])

# ==========================================
# 2. 数学与空间几何算法
# ==========================================
def gaussian(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def wavelength_to_rgb(wl):
    if 380 <= wl <= 440:
        r = -(wl - 440) / (440 - 380); g = 0.0; b = 1.0
    elif 440 < wl <= 490:
        r = 0.0; g = (wl - 440) / (490 - 440); b = 1.0
    elif 490 < wl <= 510:
        r = 0.0; g = 1.0; b = -(wl - 510) / (510 - 490)
    elif 510 < wl <= 580:
        r = (wl - 510) / (580 - 510); g = 1.0; b = 0.0
    elif 580 < wl <= 645:
        r = 1.0; g = -(wl - 645) / (645 - 580); b = 0.0
    elif 645 < wl <= 780:
        r = 1.0; g = 0.0; b = 0.0
    else:
        r = g = b = 0.0

    if 380 <= wl <= 420: factor = 0.3 + 0.7 * (wl - 380) / (420 - 380)
    elif 420 < wl <= 700: factor = 1.0
    elif 700 < wl <= 780: factor = 0.3 + 0.7 * (780 - wl) / (780 - 700)
    else: factor = 0.0
    return (r * factor, g * factor, b * factor)

def barycentric_coords(p, a, b, c):
    x, y = p
    x1, y1 = a; x2, y2 = b; x3, y3 = c
    det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    if abs(det) < 1e-8: return 0, 0, 0
    w1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det
    w2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det
    return w1, w2, 1.0 - w1 - w2

def point_in_gamut(p):
    wr1, wg1, wc1 = barycentric_coords(p, XY['R'], XY['G'], XY['C'])
    if wr1 >= -1e-4 and wg1 >= -1e-4 and wc1 >= -1e-4: return True
    wr2, wb2, wc2 = barycentric_coords(p, XY['R'], XY['B'], XY['C'])
    if wr2 >= -1e-4 and wb2 >= -1e-4 and wc2 >= -1e-4: return True
    return False

def restrict_to_gamut(p):
    if point_in_gamut(p): return p
    edges = [(XY['R'], XY['G']), (XY['G'], XY['C']), (XY['C'], XY['B']), (XY['B'], XY['R'])]
    closest_p, min_dist = None, float('inf')
    for a, b in edges:
        ap = p - a; ab = b - a
        t = max(0.0, min(1.0, np.dot(ap, ab) / np.dot(ab, ab)))
        cp = a + t * ab
        dist = np.linalg.norm(p - cp)
        if dist < min_dist:
            min_dist = dist; closest_p = cp
    return closest_p

def calc_rc_split(target_xy):
    w_r, w_g, w_c = barycentric_coords(target_xy, XY['R'], XY['G'], XY['C'])
    if w_g >= -0.001:
        return np.clip([w_r, w_g, w_c, 0.0], 0, None)
    w_r, w_b, w_c = barycentric_coords(target_xy, XY['R'], XY['B'], XY['C'])
    return np.clip([w_r, 0.0, w_c, w_b], 0, None)

def calc_gb_split(target_xy):
    w_r, w_g, w_b = barycentric_coords(target_xy, XY['R'], XY['G'], XY['B'])
    if w_r >= -0.001:
        return np.clip([w_r, w_g, 0.0, w_b], 0, None)
    w_c, w_g, w_b = barycentric_coords(target_xy, XY['C'], XY['G'], XY['B'])
    return np.clip([0.0, w_g, w_c, w_b], 0, None)

def calculate_weights(target_xy, mode):
    w_rc = calc_rc_split(target_xy)
    w_gb = calc_gb_split(target_xy)
    
    if mode == '革新派':
        weights = w_rc
    elif mode == '保守派':
        weights = w_gb
    else: 
        weights = (w_rc + w_gb) / 2.0
        
    return weights[0], weights[1], weights[2], weights[3]

# ==========================================
# 3. UI 初始化与渲染
# ==========================================
plt.style.use('dark_background')
fig = plt.figure(figsize=(13, 7))
fig.canvas.manager.set_window_title("RGBC Gamut Simulator")

# --- 左侧：色度图 ---
ax_xy = fig.add_axes([0.05, 0.25, 0.4, 0.65])
ax_xy.set_xlim(0, 0.8)
ax_xy.set_ylim(0, 0.9)
ax_xy.set_aspect('equal')
ax_xy.set_title("CIE 1931 xy Gamut", fontsize=14)
ax_xy.set_xlabel("x")
ax_xy.set_ylabel("y")
ax_xy.grid(True, alpha=0.2)

ax_xy.plot(CIE_LOCUS[:, 0], CIE_LOCUS[:, 1], '-', color='#888888', linewidth=1.5, zorder=1)
ax_xy.plot([CIE_LOCUS[-1, 0], CIE_LOCUS[0, 0]], [CIE_LOCUS[-1, 1], CIE_LOCUS[0, 1]], '--', color='#888888', linewidth=1.5, zorder=1)

poly_x = [XY['R'][0], XY['G'][0], XY['C'][0], XY['B'][0], XY['R'][0]]
poly_y = [XY['R'][1], XY['G'][1], XY['C'][1], XY['B'][1], XY['R'][1]]
ax_xy.plot(poly_x, poly_y, 'w-', alpha=0.8, zorder=2)

colors = {'R': '#FF3333', 'G': '#33FF33', 'C': '#33FFFF', 'B': '#3333FF'}
for key, pos in XY.items():
    ax_xy.plot(pos[0], pos[1], marker='o', color=colors[key], markersize=8, zorder=3)
    ax_xy.text(pos[0]+0.02, pos[1]+0.02, key, color=colors[key], fontsize=12, fontweight='bold')

initial_target = restrict_to_gamut(np.array([0.3127, 0.3290])) # D65
target_pt, = ax_xy.plot(initial_target[0], initial_target[1], 'w+', markersize=15, markeredgewidth=2, zorder=4)

# --- 右侧：光谱图 ---
ax_spd = fig.add_axes([0.55, 0.25, 0.4, 0.65])
ax_spd.set_xlim(380, 780)
ax_spd.set_ylim(0, 1.1)
ax_spd.set_title("Spectral Power Distribution", fontsize=14)
ax_spd.set_xlabel("Wavelength (nm)")
ax_spd.set_ylabel("Relative Intensity")
ax_spd.grid(True, alpha=0.2)

gradient_colors = [wavelength_to_rgb(w) for w in WAVELENGTHS]
fill_lines = ax_spd.vlines(WAVELENGTHS, 0, np.zeros_like(WAVELENGTHS), colors=gradient_colors, alpha=0.6, zorder=1)

line_r, = ax_spd.plot(WAVELENGTHS, np.zeros_like(WAVELENGTHS), color=colors['R'], label='Red', zorder=2, alpha=0.7)
line_g, = ax_spd.plot(WAVELENGTHS, np.zeros_like(WAVELENGTHS), color=colors['G'], label='Green', zorder=2, alpha=0.7)
line_c, = ax_spd.plot(WAVELENGTHS, np.zeros_like(WAVELENGTHS), color=colors['C'], label='Cyan', zorder=2, alpha=0.7)
line_b, = ax_spd.plot(WAVELENGTHS, np.zeros_like(WAVELENGTHS), color=colors['B'], label='Blue', zorder=2, alpha=0.7)
line_total, = ax_spd.plot(WAVELENGTHS, np.zeros_like(WAVELENGTHS), 'w-', linewidth=2, label='Total SPD', zorder=3)

ax_spd.legend(loc='upper right')
text_weights = ax_spd.text(0.05, 0.95, "", transform=ax_spd.transAxes, va='top', fontsize=11, family='monospace', zorder=4)

# --- 底部：UI 控件 ---
ax_radio = fig.add_axes([0.05, 0.05, 0.45, 0.12])
radio = RadioButtons(ax_radio, ('革新派', '保守派', '环保少女How dare you派'), activecolor='white')

for label in radio.labels:
    label.set_fontsize(12)

# ==========================================
# 4. 交互逻辑
# ==========================================
state = {
    'target': initial_target,
    'mode': '革新派',
    'is_dragging': False
}

def update_plot():
    target_pt.set_data([state['target'][0]], [state['target'][1]])
    
    r, g, c, b = calculate_weights(state['target'], state['mode'])
    
    spd_r = r * gaussian(WAVELENGTHS, PEAKS['R'], SIGMAS['R'])
    spd_g = g * gaussian(WAVELENGTHS, PEAKS['G'], SIGMAS['G'])
    spd_c = c * gaussian(WAVELENGTHS, PEAKS['C'], SIGMAS['C'])
    spd_b = b * gaussian(WAVELENGTHS, PEAKS['B'], SIGMAS['B'])
    spd_total = spd_r + spd_g + spd_c + spd_b
    
    line_r.set_ydata(spd_r)
    line_g.set_ydata(spd_g)
    line_c.set_ydata(spd_c)
    line_b.set_ydata(spd_b)
    line_total.set_ydata(spd_total)
    
    segs = [np.array([[w, 0], [w, h]]) for w, h in zip(WAVELENGTHS, spd_total)]
    fill_lines.set_segments(segs)
    
    text_weights.set_text(
        f"Point: x={state['target'][0]:.3f}, y={state['target'][1]:.3f}\n"
        f"----------------------\n"
        f"Drive Levels (0-1):\n"
        f"Red   : {r:5.3f}\n"
        f"Green : {g:5.3f}\n"
        f"Cyan  : {c:5.3f}\n"
        f"Blue  : {b:5.3f}"
    )
    fig.canvas.draw_idle()

def handle_mouse(event):
    if event.inaxes == ax_xy:
        state['target'] = restrict_to_gamut(np.array([event.xdata, event.ydata]))
        update_plot()

def on_radio_click(label):
    state['mode'] = label
    update_plot()

radio.on_clicked(on_radio_click)
fig.canvas.mpl_connect('button_press_event', lambda e: state.update({'is_dragging': True}) or handle_mouse(e))
fig.canvas.mpl_connect('motion_notify_event', lambda e: handle_mouse(e) if state['is_dragging'] else None)
fig.canvas.mpl_connect('button_release_event', lambda e: state.update({'is_dragging': False}))

update_plot()
plt.show()
