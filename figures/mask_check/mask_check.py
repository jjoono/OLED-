"""Schematic of the 8-pixel anode/organic/cathode mask set (dimensions read from the CAD screenshots, ~0.1 mm):
nominal overlay, the two misalignment failure modes, the ITO-island leakage path, and the fixes."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch
plt.rcParams['font.family'] = 'WenQuanYi Zen Hei'
plt.rcParams['axes.unicode_minus'] = False

W, H, CX = 21.0, 23.0, 10.45
C_AN, C_ISL, C_ORG, C_CAT, C_PAD = '#7fb3d5', '#aed6f1', '#f5b7b1', '#34495e', '#566573'
mx = lambda xs: [2 * CX - x for x in xs]
my = lambda ys: [H - y for y in ys]


def anode_polys(dx=0.0, dy=0.0):
    """central anode (two strips joined at the ends) as outer polygon + slot; cathode-lead ITO islands"""
    xo = [7.95, 12.95, 12.95, 13.94, 13.94, 12.95, 12.95, 7.95, 7.95, 6.96, 6.96, 7.95]
    yo = [0, 0, 2.5, 3.46, 19.54, 20.5, 23, 23, 20.5, 19.54, 3.46, 2.5]
    outer = [(x + dx, y + dy) for x, y in zip(xo, yo)]
    slot = [(x + dx, y + dy) for x, y in [(8.97, 4.97), (11.93, 4.97), (11.93, 18.03), (8.97, 18.03)]]
    isl = []
    left = [[(3.96, 0), (6.46, 0), (6.46, 8.0), (3.96, 8.0)],                     # under the top-left C lead
            [(1.2, 8.99), (6.46, 8.99), (6.46, 10.95), (1.2, 10.95)],             # middle stubs (lead runs to the edge)
            [(1.2, 12.05), (6.46, 12.05), (6.46, 14.01), (1.2, 14.01)]]
    left.append([(x, H - y) for x, y in left[0]])
    for p in left:
        isl.append([(x + dx, y + dy) for x, y in p])
        isl.append([(2 * CX - x + dx, y + dy) for x, y in p])
    return outer, slot, isl


def cathode_polys():
    tl = [(3.96, 0), (6.91, 0), (6.91, 2.5), (5.97, 2.5), (5.97, 6.0), (9.98, 6.0), (9.98, 8.0), (3.96, 8.0)]
    m1 = [(1.2, 8.99), (9.98, 8.99), (9.98, 10.95), (1.2, 10.95)]
    m2 = [(1.2, 12.05), (9.98, 12.05), (9.98, 14.01), (1.2, 14.01)]
    bl = [(x, H - y) for x, y in tl]
    out = []
    for p in (tl, m1, m2, bl):
        out += [p, [(2 * CX - x, y) for x, y in p]]
    pads = [[(8.0, 0), (12.9, 0), (12.9, 2.5), (8.0, 2.5)], [(8.0, 20.5), (12.9, 20.5), (12.9, 23), (8.0, 23)]]
    return out, pads


ORG = (5.93, 4.83, 14.98, 18.17)


def draw(ax, dx=0.0, dy=0.0, org=ORG, title='', marks=(), note='', ncol='#c0392b'):
    sub = Rectangle((0, 0), W, H, fc='white', ec='k', lw=1); ax.add_patch(sub)
    n0 = len(ax.patches)
    outer, slot, isl = anode_polys(dx, dy)
    for p in isl: ax.add_patch(Polygon(p, fc=C_ISL, ec='#5dade2', lw=0.6))
    ax.add_patch(Polygon(outer, fc=C_AN, ec='#2471a3', lw=0.8))
    ax.add_patch(Polygon(slot, fc='white', ec='#2471a3', lw=0.6))
    x0, y0, x1, y1 = org
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=C_ORG, ec='#c0392b', lw=1.2, alpha=0.55))
    cat, pads = cathode_polys()
    for p in cat: ax.add_patch(Polygon(p, fc=C_CAT, ec='none', alpha=0.78))
    for p in pads: ax.add_patch(Polygon(p, fc=C_PAD, ec='none', alpha=0.78, hatch='///'))
    for p in ax.patches[n0:]: p.set_clip_path(sub)
    if note: ax.text(0, H + 1.0, note, fontsize=8.6, color=ncol, fontweight='bold', va='top')
    for (x, y, s) in marks:
        ax.plot(x, y, marker='*', ms=17, mfc='#e74c3c', mec='k', mew=0.8, zorder=9)
        if s: ax.annotate(s[0], (x, y), s[1], fontsize=8.5, color='#c0392b', fontweight='bold', zorder=10,
                          bbox=dict(fc='white', ec='#c0392b', lw=0.6, pad=1.5),
                          arrowprops=dict(arrowstyle='-', color='#c0392b', lw=0.9))
    ax.set_xlim(-0.6, W + 0.6); ax.set_ylim(H + 5.2, -0.6); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, fontsize=10.5, loc='left', fontweight='bold')


def dim(ax, x0, x1, y, txt, col='#1e8449', vertical=False, off=0.0):
    if vertical:
        ax.add_patch(FancyArrowPatch((y, x0), (y, x1), arrowstyle='<->', mutation_scale=8, color=col, lw=1.1, zorder=8))
        ax.text(y + 0.25 + off, (x0 + x1) / 2, txt, color=col, fontsize=8, va='center', fontweight='bold', zorder=8)
    else:
        ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='<->', mutation_scale=8, color=col, lw=1.1, zorder=8))
        ax.text((x0 + x1) / 2, y - 0.35 + off, txt, color=col, fontsize=8, ha='center', fontweight='bold', zorder=8)


fig = plt.figure(figsize=(16, 12.2))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.10, wspace=0.08)

# (a) nominal
ax = fig.add_subplot(gs[0, 0])
draw(ax, title='(a) 설계 위치 (정렬 OK)', ncol='#1e8449',
     note='노란 사각형 = 픽셀 8개 (약 2×2 mm)\n초록 = anode가 움직일 수 있는 여유 (약 1 mm)\n보라 = ITO 섬 ↔ anode 틈 0.5 mm (정렬과 무관)')
for (x, y) in ((7.97, 7.0), (12.93, 7.0), (7.97, 10.0), (12.93, 10.0), (7.97, 13.0), (12.93, 13.0), (7.97, 16.0), (12.93, 16.0)):
    ax.add_patch(Rectangle((x - 1.0, y - 1.0), 2.0, 2.0, fc='none', ec='#f1c40f', lw=1.8, zorder=7))
dim(ax, 5.93, 6.96, 11.5, '')
ax.annotate('유기층 끝 ↔ anode 띠\n1.0 mm', (6.45, 11.5), (-0.4, 15.6), fontsize=8, color='#1e8449', fontweight='bold',
            bbox=dict(fc='white', ec='#1e8449', lw=0.6, pad=1.5), arrowprops=dict(arrowstyle='-', color='#1e8449', lw=0.9), zorder=10)
dim(ax, 5.97, 6.96, 3.4, '')
ax.annotate('C자 리드 ↔ anode 끝\n1.0 mm (유기층 밖)', (6.45, 3.4), (13.8, 2.6), fontsize=8, color='#1e8449', fontweight='bold',
            bbox=dict(fc='white', ec='#1e8449', lw=0.6, pad=1.5), arrowprops=dict(arrowstyle='-', color='#1e8449', lw=0.9), zorder=10)
dim(ax, 4.83, 6.0, 10.45, '', vertical=True)
ax.annotate('유기층 끝 ↔ 손가락\n1.2 mm', (10.45, 5.4), (14.6, 5.5), fontsize=8, color='#1e8449', fontweight='bold',
            bbox=dict(fc='white', ec='#1e8449', lw=0.6, pad=1.5), arrowprops=dict(arrowstyle='-', color='#1e8449', lw=0.9), zorder=10)
dim(ax, 6.46, 6.96, 10.0, '', col='#8e44ad')
ax.annotate('0.5 mm', (6.7, 10.0), (0.0, 7.6), fontsize=8, color='#8e44ad', fontweight='bold',
            bbox=dict(fc='white', ec='#8e44ad', lw=0.6, pad=1.5), arrowprops=dict(arrowstyle='-', color='#8e44ad', lw=0.9), zorder=10)

# (b) x shift
ax = fig.add_subplot(gs[0, 1])
draw(ax, dx=-1.3, title='(b) anode가 왼쪽으로 1.3 mm 밀림 → 단락',
     marks=[(5.8, 7.0, ('손가락이 유기층 밖에서\nanode ITO에 바로 닿음', (11.5, 3.9))), (5.8, 10.0, None), (5.8, 13.0, None), (5.8, 16.0, None),
            (6.78, 1.3, ('C자 리드 ↔ anode 끝 겹침', (9.3, -0.1))), (6.78, 21.7, None)],
     note='유기층 왼쪽 끝(5.93 mm)보다 anode 띠(5.66 mm)가 밖으로 나옴\n→ 왼쪽 픽셀 4개 단락 + 위/아래 C자 리드도 anode 끝에 겹침\n오른쪽으로 밀리면 오른쪽 4개가 같은 식으로 단락')

# (c) y shift
ax = fig.add_subplot(gs[0, 2])
draw(ax, dy=-1.0, title='(c) anode가 위로 1.0 mm 밀림 → 모서리 단락',
     marks=[(6.93, 2.5, ('간격 ≈ 0.05 mm', (8.6, 4.3))), (13.97, 2.5, None)],
     note='anode 띠 끝(모따기)과 C자 리드는 서로 피하도록 설계되어 있음\n위로 1 mm 밀리면 넓은 부분끼리 같은 높이 → 간격 ≈ 0.05 mm\n이 곳은 유기층(4.8 mm) 위쪽 = 맨 ITO → 조금만 옆으로 어긋나도 단락\n(대각선으로 밀리면 가장 위험, 아래로 밀리면 아래 모서리)')
ax.add_patch(Rectangle((4.4, 0.05), 12.1, 4.6, fc='none', ec='#c0392b', lw=1.1, ls='--', zorder=8))

# (d) cross-section of the island leakage path
ax = fig.add_subplot(gs[1, 0])
ax.set_title('(d) 정렬과 무관: ITO 섬 ↔ anode 측면 누설 (단면)', fontsize=10.5, loc='left', fontweight='bold')
ax.add_patch(Rectangle((0, 0), 10, 0.6, fc='#d5d8dc', ec='k', lw=0.6)); ax.text(5, 0.22, '유리 기판', ha='center', fontsize=8)
ax.add_patch(Rectangle((0, 0.6), 4.2, 0.35, fc=C_ISL, ec='#2471a3')); ax.text(2.4, 0.7, 'ITO 섬 (cathode 전위)', ha='center', fontsize=7.6)
ax.add_patch(Rectangle((5.2, 0.6), 4.8, 0.35, fc=C_AN, ec='#2471a3')); ax.text(7.6, 0.7, 'anode ITO 띠', ha='center', fontsize=7.6)
ax.annotate('0.5 mm 틈', (4.7, 0.78), (4.7, -0.55), ha='center', fontsize=8, color='#8e44ad', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=0.9))
ax.add_patch(Rectangle((0.8, 0.95), 9.2, 0.25, fc='#f9e79f', ec='#b7950b')); ax.text(9.9, 1.0, 'HATCN', ha='right', fontsize=7.6)
ax.add_patch(Rectangle((0.8, 1.2), 9.2, 0.6, fc='#fdebd0', ec='#ca6f1e')); ax.text(9.9, 1.55, 'TAPC …', ha='right', fontsize=7.6)
ax.add_patch(Rectangle((0.8, 1.8), 9.2, 0.6, fc=C_ORG, ec='#c0392b', alpha=0.7)); ax.text(9.9, 2.05, 'EML / ETL', ha='right', fontsize=7.6)
ax.add_patch(Rectangle((0, 2.4), 10, 0.35, fc=C_CAT, ec='none', alpha=0.85)); ax.text(5, 2.5, 'cathode Ag', ha='center', color='w', fontsize=8)
ax.add_patch(Rectangle((0, 0.95), 0.8, 1.45, fc=C_CAT, ec='none', alpha=0.85))
ax.text(0.4, 1.55, '유기층\n밖', ha='center', color='w', fontsize=6.8)
ax.annotate('', (1.4, 1.27), (8.8, 1.27), arrowprops=dict(arrowstyle='->', color='#c0392b', lw=2.4))
ax.text(5.0, 1.42, '정공이 HATCN/TAPC 계면을 따라 옆으로 이동', ha='center', fontsize=7.8, color='#c0392b', fontweight='bold')
ax.text(0.0, -1.1, '유기층 밖에서 cathode 금속이 ITO 섬에 붙어 있으므로\nanode → 계면 → ITO 섬 = EML을 거치지 않는 누설 경로\n어림: 4 V에서 수십 nA – 1 µA (저휘도 픽셀 전류 ~0.1 µA와 비슷)',
        fontsize=8.5, color='#2c3e50', va='top')
ax.text(0.0, -2.75, '→ cathode 올리기 전에 anode 컨택 ↔ ITO 섬 I–V(±4 V) 측정\n   1 nA 미만이면 안전', fontsize=8.8, color='#1e8449', fontweight='bold', va='top')
ax.set_xlim(-0.3, 10.3); ax.set_ylim(-3.7, 3.0); ax.axis('off')

# (e) fixes: enlarged organic + shims
ax = fig.add_subplot(gs[1, 1])
draw(ax, org=(5.4, 3.0, 15.5, 20.0), title='(e) 대책 ①: 유기층 영역 확대 (예시)', ncol='#1e8449',
     note='점선 = 현재(x 5.9–15.0, y 4.8–18.2), 실선 = 확대안(x 5.4–15.5, y 3–20)\n위/아래 C자 리드 근처까지 덮음 → 어긋나도 단락 대신 작은 기생 다이오드\n좌우 여유 1.0 → 1.5 mm. anode 컨택 패드(맨 위/아래)는 계속 노출')
ax.add_patch(Rectangle((ORG[0], ORG[1]), ORG[2] - ORG[0], ORG[3] - ORG[1], fc='none', ec='#c0392b', lw=1.2, ls=':', zorder=8))

ax = fig.add_subplot(gs[1, 2])
ax.set_title('(f) 대책 ②: 23 mm 기판을 심(shim)으로 고정', fontsize=10.5, loc='left', fontweight='bold')
ax.add_patch(Rectangle((0, 0), 24, 24, fc='#eaecee', ec='k', lw=1.2))
ax.add_patch(Rectangle((0.6, 0.4), 23, 23, fc='white', ec='#2471a3', lw=1.2)); ax.text(12.1, 11.0, '23 × 23 mm 기판', ha='center', fontsize=9.5, color='#2471a3')
ax.add_patch(Rectangle((0, 0), 0.6, 24, fc='#f39c12', ec='k', lw=0.5)); ax.add_patch(Rectangle((0, 0), 24, 0.4, fc='#f39c12', ec='k', lw=0.5))
ax.annotate('심 (두 변)', (0.3, 6.0), (3.0, 5.0), fontsize=8.5, color='#a04000', fontweight='bold', arrowprops=dict(arrowstyle='->', color='#a04000'))
ax.text(3.0, 7.6, '두께 = 설계 위치 - 실측 위치\n(예: 0.4–0.6 mm, Kapton / Al 포일 / 얇은 유리)', fontsize=8.3, color='#a04000', va='top')
ax.text(12.1, 14.5, '두 변을 기준면에 밀착\n→ 위치와 회전을 동시에 고정, 목표 오차 ≤ 0.3 mm\n넣기 전 라이트 테이블에서 마스크와 겹쳐 확인', ha='center', fontsize=8.6, color='#1e8449', fontweight='bold')
ax.text(0.0, 25.2, '홀더 자리 24 × 24 mm (회색)', fontsize=8.3, color='#566573', va='top')
ax.text(0.0, 26.8, '정렬이 불확실하면: 모서리 C자 리드 4개를 Kapton으로 막고 가운데 픽셀 4개만 사용\n측정 전 각 픽셀 -2~-3 V로 단락·누설 픽셀 선별', fontsize=8.6, color='#2c3e50', va='top')
ax.set_xlim(-0.6, 24.6); ax.set_ylim(29.8, -0.6); ax.set_aspect('equal'); ax.axis('off')

handles = [Rectangle((0, 0), 1, 1, fc=C_AN, ec='#2471a3'), Rectangle((0, 0), 1, 1, fc=C_ISL, ec='#5dade2'),
           Rectangle((0, 0), 1, 1, fc=C_ORG, ec='#c0392b', alpha=0.55), Rectangle((0, 0), 1, 1, fc=C_CAT, alpha=0.78),
           Rectangle((0, 0), 1, 1, fc=C_PAD, alpha=0.78, hatch='///')]
fig.legend(handles, ['anode ITO (공통, 가운데 띠 2줄)', 'cathode 리드용 ITO 섬', '유기층 영역', 'cathode 금속', 'anode 컨택 금속 패드 (cathode 마스크)'],
           loc='upper center', ncol=5, fontsize=9.5, frameon=False, bbox_to_anchor=(0.5, 0.955))
fig.text(0.5, 0.085, '치수는 CAD 화면에서 읽은 값 (±0.1 mm). 위는 기판 위쪽, 단위 mm (21 × 23 패턴 영역 기준).', ha='center', fontsize=8.5, color='#566573')
fig.savefig('mask_alignment_check.png', dpi=160, bbox_inches='tight')
print('ok')
