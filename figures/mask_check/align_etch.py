"""Anode trim that removes the island leakage AND widens the misalignment tolerance (organic and cathode masks unchanged).
Keeps: a 3.9 mm neck (x 8.5-12.4) along the whole length, the two pixel strips trimmed to x 7.5-8.97 / 11.93-13.4 between
y 5.0 and 18.0, and the outer cathode-contact ITO (x < 3.9; C-lead contacts x 3.96-5.5 at y < 3.3 and y > 19.7).
Tolerance maps: anode shifted by (dx, dy) relative to the (holder-referenced) organic and cathode masks, raster 0.025 mm."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch
from matplotlib.colors import ListedColormap
plt.rcParams['font.family'] = 'WenQuanYi Zen Hei'; plt.rcParams['axes.unicode_minus'] = False

W, H, CX = 21.0, 23.0, 10.45
C_AN, C_ISL, C_ORG, C_CAT, C_PAD = '#7fb3d5', '#aed6f1', '#f5b7b1', '#34495e', '#566573'
R = lambda x0, y0, x1, y1: [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
mir = lambda p: [(2 * CX - x, y) for x, y in p]
both = lambda ps: ps + [mir(p) for p in ps]
ORG = (5.93, 4.83, 14.98, 18.17)
STRIP = [(7.95, 0), (12.95, 0), (12.95, 2.5), (13.94, 3.46), (13.94, 19.54), (12.95, 20.5), (12.95, 23), (7.95, 23),
         (7.95, 20.5), (6.96, 19.54), (6.96, 3.46), (7.95, 2.5)]
SLOT = R(8.97, 4.97, 11.93, 18.03)
ISL_OLD = both([R(3.96, 0, 6.46, 8.0), R(1.2, 8.99, 6.46, 10.95), R(1.2, 12.05, 6.46, 14.01), R(3.96, 15.0, 6.46, 23)])
NEW_CENTRAL = [R(8.5, 0, 12.4, 23), R(7.5, 5.0, 13.4, 18.0)]
ISL_NEW = both([R(3.96, 0, 5.5, 3.3), R(1.2, 8.99, 3.9, 10.95), R(1.2, 12.05, 3.9, 14.01), R(3.96, 19.7, 5.5, 23)])
WIN = both([R(3.9, 3.3, 7.5, 19.7), R(5.5, 0, 8.5, 5.0), R(5.5, 18.0, 8.5, 23)])        # etch windows (union)
TL = [(3.96, 0), (6.91, 0), (6.91, 2.5), (5.97, 2.5), (5.97, 6.0), (9.98, 6.0), (9.98, 8.0), (3.96, 8.0)]
CAT = both([TL, R(1.2, 8.99, 9.98, 10.95), R(1.2, 12.05, 9.98, 14.01), [(x, H - y) for x, y in TL]])
PADS = [R(8.0, 0, 12.9, 2.5), R(8.0, 20.5, 12.9, 23)]

# ---------------- raster tolerance analysis ----------------
D = 0.025; PADPX = 100
xs = np.arange(-PADPX, int(W / D) + PADPX) * D + D / 2; ys = np.arange(-PADPX, int(H / D) + PADPX) * D + D / 2
XX, YY = np.meshgrid(xs, ys); PTS = np.c_[XX.ravel(), YY.ravel()]
def ras(polys, holes=()):
    m = np.zeros(PTS.shape[0], bool)
    for p in polys: m |= Path(p).contains_points(PTS)
    for p in holes: m &= ~Path(p).contains_points(PTS)
    return m.reshape(XX.shape)
CATM, PADM = ras(CAT), ras(PADS); ORGM = ras([R(*ORG)])
CAT_OUT = CATM & ~ORGM
def shift(m, k, l):
    o = np.zeros_like(m); h, w = m.shape
    o[max(0, l):h + min(0, l), max(0, k):w + min(0, k)] = m[max(0, -l):h - max(0, l), max(0, -k):w - max(0, k)]
    return o
def tolmap(central, islands, step=0.05, rng=2.0):
    AC, IS = ras(central, [SLOT]), ras(islands)
    s = np.arange(-rng, rng + 1e-9, step); out = np.zeros((len(s), len(s)), int)
    for i, dy in enumerate(s):
        for j, dx in enumerate(s):
            k, l = int(round(dx / D)), int(round(dy / D))
            a, isl = shift(AC, k, l), shift(IS, k, l)
            short = (CAT_OUT & a).any() or (PADM & isl).any()
            leak = (isl & ORGM).any()
            out[i, j] = 2 if short else (1 if leak else 0)
    return s, out
s, T_old = tolmap([STRIP], ISL_OLD)
_, T_new = tolmap(NEW_CENTRAL, ISL_NEW)
def radius(T, lvl=1):                       # largest shift radius free of state >= lvl in every direction
    rr = np.hypot(*np.meshgrid(s, s)); bad = rr[T >= lvl]
    return bad.min() if bad.size else np.inf
def axis_lim(T):
    i0 = len(s) // 2; row, col = T[i0], T[:, i0]
    f = lambda v, sgn: min([abs(x) for x, t in zip(s, v) if t > 0 and np.sign(x) == sgn] or [np.inf])
    return f(row, -1), f(row, 1), f(col, -1), f(col, 1)
print('short-only radius old %.2f new %.2f' % (radius(T_old, 2), radius(T_new, 2)))
print('old: centre state', T_old[len(s)//2, len(s)//2], 'safe radius %.2f' % radius(T_old), 'axis limits (-x,+x,-y,+y)', axis_lim(T_old))
print('new: centre state', T_new[len(s)//2, len(s)//2], 'safe radius %.2f' % radius(T_new), 'axis limits', axis_lim(T_new))

# ---------------- figure ----------------
def frame(ax, title, extra=5.2):
    sub = Rectangle((0, 0), W, H, fc='white', ec='k', lw=1); ax.add_patch(sub)
    ax.set_xlim(-0.6, W + 0.6); ax.set_ylim(H + extra, -0.6); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, fontsize=11, loc='left', fontweight='bold'); return sub
def poly(ax, ps, **kw):
    for p in ps: ax.add_patch(Polygon(p, **kw))
def central_new(ax):
    poly(ax, NEW_CENTRAL, fc=C_AN, ec='#2471a3', lw=0.9); poly(ax, [SLOT], fc='white', ec='#2471a3', lw=0.7)
def overlay(ax):
    ax.add_patch(Rectangle(ORG[:2], ORG[2] - ORG[0], ORG[3] - ORG[1], fc=C_ORG, ec='#c0392b', lw=1.2, alpha=0.5))
    poly(ax, CAT, fc=C_CAT, ec='none', alpha=0.72); poly(ax, PADS, fc=C_PAD, ec='none', alpha=0.72, hatch='///')
def box(ax, x, y, t, col, size=8.4):
    ax.text(x, y, t, fontsize=size, color=col, fontweight='bold', va='top', bbox=dict(fc='white', ec=col, lw=0.7, pad=2.5))
def dimh(ax, x0, x1, y, t, col='#1e8449', dy=-0.35):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='<->', mutation_scale=8, color=col, lw=1.2, zorder=9))
    ax.text((x0 + x1) / 2, y + dy, t, ha='center', fontsize=7.8, color=col, fontweight='bold', zorder=9)

fig = plt.figure(figsize=(18.5, 12.6))
gs = fig.add_gridspec(2, 3, wspace=0.06, hspace=0.12)

# (a) where to etch, on the present anode
ax = fig.add_subplot(gs[0, 0]); frame(ax, '(a) 지울 곳 (빨간 빗금) — 지금 anode 위')
poly(ax, ISL_OLD, fc=C_ISL, ec='#5dade2', lw=0.7); poly(ax, [STRIP], fc=C_AN, ec='#2471a3', lw=0.9); poly(ax, [SLOT], fc='white', ec='#2471a3', lw=0.7)
poly(ax, WIN, fc='#e74c3c', ec='#922b21', lw=0.9, alpha=0.38, hatch='xx', zorder=6)
for x, t in ((3.9, '3.9'), (5.5, '5.5'), (7.5, '7.5'), (8.5, '8.5'), (12.4, '12.4'), (13.4, '13.4'), (15.4, '15.4'), (17.0, '17.0')):
    ax.plot([x, x], [-0.5, H + 0.4], color='#922b21', lw=0.5, ls=':', zorder=7)
    ax.text(x, H + 0.95, t, ha='center', fontsize=7.4, color='#922b21', rotation=90, va='top')
for y, t in ((3.3, '3.3'), (5.0, '5.0'), (18.0, '18.0'), (19.7, '19.7')):
    ax.plot([-0.5, W + 0.5], [y, y], color='#922b21', lw=0.5, ls=':', zorder=7); ax.text(W + 0.3, y, t, fontsize=7.4, color='#922b21', va='center')
box(ax, 0, H + 3.0, '창(좌): x 3.9–7.5 × y 3.3–19.7  +  x 5.5–8.5 × y 0–5.0 / 18.0–23\n창(우): 좌우 대칭 (x 13.4–17.0, x 12.4–15.4)   단위 mm, 21×23 패턴 기준', '#922b21', 8.2)

# (b) tape layout
ax = fig.add_subplot(gs[0, 1]); frame(ax, '(b) 테이프로 가릴 곳 (주황) = 남길 곳')
poly(ax, ISL_OLD, fc=C_ISL, ec='#5dade2', lw=0.7); poly(ax, [STRIP], fc=C_AN, ec='#2471a3', lw=0.9); poly(ax, [SLOT], fc='white', ec='#2471a3', lw=0.7)
TAPE = [R(8.5, -0.3, 12.4, 23.3), R(7.5, 5.0, 13.4, 18.0)] + both([R(-0.3, -0.3, 3.9, 23.3), R(3.9, -0.3, 5.5, 3.3), R(3.9, 19.7, 5.5, 23.3)])
poly(ax, TAPE, fc='#f5cba7', ec='#ca6f1e', lw=0.8, alpha=0.8, zorder=7)
ax.text(10.45, 2.0, '①', ha='center', fontsize=12, color='#784212', fontweight='bold', zorder=9)
ax.text(10.45, 11.5, '②', ha='center', fontsize=12, color='#784212', fontweight='bold', zorder=9)
ax.text(1.9, 11.5, '③', ha='center', fontsize=12, color='#784212', fontweight='bold', zorder=9)
ax.text(2 * CX - 1.9, 11.5, '③', ha='center', fontsize=12, color='#784212', fontweight='bold', zorder=9)
ax.text(4.7, 1.6, '④', ha='center', fontsize=10, color='#784212', fontweight='bold', zorder=9)
box(ax, 0, H + 1.0, '① 목: 폭 3.9 mm 띠, 가운데 슬롯(유리 창) 중심에 맞춰 위→아래 전체\n'
    '② 픽셀 띠: 폭 5.9 mm × 길이 13 mm, 슬롯 위아래 끝(y 5.0 / 18.0)에 맞춤\n'
    '   = 슬롯 양옆으로 1.47 mm 씩 (기존 띠 바깥 0.54 mm 깎임)\n'
    '③ 바깥 컨택: 기판 좌/우 끝 ~ x 3.9 (섬의 바깥 2.7 mm)\n'
    '④ C자 리드 컨택: 섬 바깥 1.54 mm 폭 × 위/아래 끝에서 3.3 mm', '#784212', 8.2)

# (c) result + overlay
ax = fig.add_subplot(gs[0, 2]); frame(ax, '(c) 결과 anode + 유기층·cathode (그대로)')
poly(ax, ISL_NEW, fc=C_ISL, ec='#5dade2', lw=0.7); central_new(ax); overlay(ax)
dimh(ax, 5.93, 7.5, 12.5 - 1.0 + 0.0, '1.57')
dimh(ax, 5.97, 8.5, 3.6, '2.5')
dimh(ax, 6.91, 8.5, 1.3, '1.59')
dimh(ax, 3.9, 5.93, 15.0 - 0.45, '2.0')
for (x, y) in ((8.24, 7.0), (12.67, 7.0), (8.24, 10.0), (12.67, 10.0), (8.24, 13.0), (12.67, 13.0), (8.24, 16.0), (12.67, 16.0)):
    ax.add_patch(Rectangle((x - 0.735, y - 0.98), 1.47, 1.96, fc='none', ec='#f1c40f', lw=1.7, zorder=8))
box(ax, 0, H + 1.0, '유기층 아래 cathode 전위 ITO 없음 → 섬 누설 소멸\n'
    '초록 = 금속↔ITO 최소 간격 (모두 ≥ 1.5 mm, 기존 1.0 mm / 위로 밀릴 때 0.05 mm)\n'
    '픽셀(노랑) 1.47 × 1.96 mm ≈ 2.9 mm² (면적은 현미경 실측)\n'
    'anode 컨택(맨 위/아래 가운데)은 목(3.9 mm) 위 금속 패드로 그대로', '#1e8449', 8.2)

# (d)(e) tolerance maps
cm = ListedColormap(['#abebc6', '#d7bde2', '#f1948a'])
for k, (T, ttl) in enumerate(((T_old, '(d) 지금 anode: anode 위치 오차 허용 지도'), (T_new, '(e) 식각 후 anode: 허용 지도'))):
    ax = fig.add_subplot(gs[1, k])
    ax.imshow(T, origin='upper', extent=[s[0], s[-1], s[-1], s[0]], cmap=cm, vmin=0, vmax=2, interpolation='nearest')
    for r in (0.5, 1.0, 1.5):
        ax.add_patch(plt.Circle((0, 0), r, fc='none', ec='k', lw=0.6, ls='--'))
        ax.text(r * 0.72, -r * 0.72, '%.1f' % r, fontsize=7.5)
    ax.axhline(0, color='k', lw=0.4); ax.axvline(0, color='k', lw=0.4)
    ax.set_xlabel('anode x 오차 (mm, + = 오른쪽)', fontsize=9); ax.set_ylabel('anode y 오차 (mm, + = 아래쪽)', fontsize=9)
    ax.set_title(ttl, fontsize=11, loc='left', fontweight='bold'); ax.set_aspect('equal')
    ax.text(0.02, 0.02, '단락 없는 반경 %.2f mm / 누설까지 없는 반경 %.2f mm' % (radius(T, 2), radius(T)), transform=ax.transAxes, fontsize=9, fontweight='bold',
            bbox=dict(fc='white', ec='k', lw=0.6))
ax = fig.add_subplot(gs[1, 2]); ax.axis('off')
ax.set_title('(f) 읽는 법과 작업 순서', fontsize=11, loc='left', fontweight='bold')
for i, (c, t) in enumerate((('#abebc6', '초록: 단락·누설 없음'), ('#d7bde2', '보라: 섬이 유기층 아래로 들어감 → 측면 누설'), ('#f1948a', '빨강: cathode 금속이 맨 anode ITO에 닿음 → 단락'))):
    ax.add_patch(Rectangle((0.0, 0.93 - i * 0.07), 0.06, 0.045, fc=c, ec='k', lw=0.5, transform=ax.transAxes))
    ax.text(0.08, 0.95 - i * 0.07, t, transform=ax.transAxes, fontsize=9.5, va='center')
ax.text(0.0, 0.66, '작업 순서\n'
        '1. 뒷면(MLA/원뿔) 보호 필름\n'
        '2. 실체현미경 아래서 ①→②→③④ 순서로 Kapton\n'
        '   (IZO는 비스듬히 비춰야 보임; 슬롯 끝이 y 기준점)\n'
        '3. 노출 창만 묽은 HCl(1:5–1:10) 또는 옥살산, 짧게\n'
        '   (면봉이면 폼 스왑 + 최소량, 테이프 경계가 둑 역할)\n'
        '4. DI 헹굼 → 테이프 제거 → 아세톤/IPA/DI 초음파\n'
        '5. 멀티미터: 목↔섬 개방, 창 안 개방, 섬 바깥부 도통,\n'
        '   목 저항(위 컨택↔아래 컨택) 기록\n'
        '6. 증착 직전 UV-오존, 첫 배치는 미식각 기판 대조군\n\n'
        '테이프 허용 오차 ≈ ±0.25 mm (섬 쪽 잔여 = 고립 조각, 무해;\n'
        '띠 쪽 과식각 = 픽셀 폭 감소, 면적 실측으로 보정)\n'
        '이 식각 후에도 심으로 0.3 mm 이내 정렬하면 여유가 가장 큼', transform=ax.transAxes, fontsize=9.3, va='top', color='#2c3e50')

h = [Rectangle((0, 0), 1, 1, fc=C_AN, ec='#2471a3'), Rectangle((0, 0), 1, 1, fc=C_ISL, ec='#5dade2'),
     Rectangle((0, 0), 1, 1, fc='#e74c3c', alpha=0.38, hatch='xx', ec='#922b21'), Rectangle((0, 0), 1, 1, fc='#f5cba7', ec='#ca6f1e'),
     Rectangle((0, 0), 1, 1, fc=C_ORG, ec='#c0392b', alpha=0.5), Rectangle((0, 0), 1, 1, fc=C_CAT, alpha=0.72)]
fig.legend(h, ['anode ITO (가운데)', 'cathode 컨택용 ITO', '식각 창', 'Kapton (남길 곳)', '유기층 (그대로)', 'cathode (그대로)'],
           loc='upper center', ncol=6, fontsize=9.5, frameon=False, bbox_to_anchor=(0.5, 0.935))
fig.savefig('anode_align_etch.png', dpi=150, bbox_inches='tight')
print('saved')
