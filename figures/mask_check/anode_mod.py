"""Recommended anode modification with the organic and cathode masks unchanged: etch the cathode-lead ITO islands
only where they lie within ~1.5 mm of the organic area, keep the central anode strip and the outer contact parts."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch, Circle
plt.rcParams['font.family'] = 'WenQuanYi Zen Hei'; plt.rcParams['axes.unicode_minus'] = False

W, H, CX = 21.0, 23.0, 10.45
C_AN, C_ISL, C_ORG, C_CAT, C_PAD = '#7fb3d5', '#aed6f1', '#f5b7b1', '#34495e', '#566573'
ORG = (5.93, 4.83, 14.98, 18.17)
BAND_L = (3.9, 3.3, 6.7, 19.7)                       # etch window, left (x0, y0, x1, y1); right = mirror
BAND_R = (2 * CX - 6.7, 3.3, 2 * CX - 3.9, 19.7)
R = lambda x0, y0, x1, y1: [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
mirror = lambda p: [(2 * CX - x, y) for x, y in p]

STRIP = [(7.95, 0), (12.95, 0), (12.95, 2.5), (13.94, 3.46), (13.94, 19.54), (12.95, 20.5), (12.95, 23), (7.95, 23),
         (7.95, 20.5), (6.96, 19.54), (6.96, 3.46), (7.95, 2.5)]
SLOT = R(8.97, 4.97, 11.93, 18.03)
ISL_OLD = [R(3.96, 0, 6.46, 8.0), R(1.2, 8.99, 6.46, 10.95), R(1.2, 12.05, 6.46, 14.01), R(3.96, 15.0, 6.46, 23)]
ISL_NEW = [R(3.96, 0, 6.46, 3.3), R(1.2, 8.99, 3.9, 10.95), R(1.2, 12.05, 3.9, 14.01), R(3.96, 19.7, 6.46, 23)]
ISL_GONE = [R(3.96, 3.3, 6.46, 8.0), R(3.9, 8.99, 6.46, 10.95), R(3.9, 12.05, 6.46, 14.01), R(3.96, 15.0, 6.46, 19.7)]
TL = [(3.96, 0), (6.91, 0), (6.91, 2.5), (5.97, 2.5), (5.97, 6.0), (9.98, 6.0), (9.98, 8.0), (3.96, 8.0)]
CAT = [TL, R(1.2, 8.99, 9.98, 10.95), R(1.2, 12.05, 9.98, 14.01), [(x, H - y) for x, y in TL]]
PADS = [R(8.0, 0, 12.9, 2.5), R(8.0, 20.5, 12.9, 23)]


def frame(ax, title, extra=4.8):
    sub = Rectangle((0, 0), W, H, fc='white', ec='k', lw=1); ax.add_patch(sub)
    ax.set_xlim(-0.6, W + 0.6); ax.set_ylim(H + extra, -0.6); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, fontsize=11, loc='left', fontweight='bold')
    return sub


def anode(ax, islands, gone=(), gone_style='hatch'):
    for p in islands:
        for q in (p, mirror(p)): ax.add_patch(Polygon(q, fc=C_ISL, ec='#5dade2', lw=0.7))
    for p in gone:
        for q in (p, mirror(p)):
            ax.add_patch(Polygon(q, fc='white', ec='#c0392b', lw=0.9, ls='--', hatch='xxx' if gone_style == 'hatch' else None))
    ax.add_patch(Polygon(STRIP, fc=C_AN, ec='#2471a3', lw=0.9)); ax.add_patch(Polygon(SLOT, fc='white', ec='#2471a3', lw=0.7))


def overlay(ax):
    x0, y0, x1, y1 = ORG
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=C_ORG, ec='#c0392b', lw=1.2, alpha=0.5))
    for p in CAT:
        for q in (p, mirror(p)): ax.add_patch(Polygon(q, fc=C_CAT, ec='none', alpha=0.72))
    for p in PADS: ax.add_patch(Polygon(p, fc=C_PAD, ec='none', alpha=0.72, hatch='///'))


def box(ax, x, y, s, col, size=8.6):
    ax.text(x, y, s, fontsize=size, color=col, fontweight='bold', va='top', bbox=dict(fc='white', ec=col, lw=0.7, pad=2.5))


def dimh(ax, x0, x1, y, col='#1e8449'):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='<->', mutation_scale=8, color=col, lw=1.2, zorder=9))


fig = plt.figure(figsize=(17.5, 6.9))
gs = fig.add_gridspec(1, 4, wspace=0.07)

# (a) now
ax = fig.add_subplot(gs[0, 0]); frame(ax, '(a) 지금 anode (스퍼터 완료)')
anode(ax, ISL_OLD); overlay(ax)
for y in (7.0, 10.0, 13.0, 16.0):
    for x in (6.2, 2 * CX - 6.2):
        ax.add_patch(Circle((x, y), 0.55, fc='none', ec='#8e44ad', lw=2.0, zorder=9))
box(ax, 0, H + 0.9, '보라 원 = 유기층 아래로 들어간 ITO 섬 끝\n(cathode 전위) → anode 띠에서 0.5 mm\n→ HATCN/TAPC 계면 따라 누설', '#8e44ad')

# (b) recommended anode alone
ax = fig.add_subplot(gs[0, 1]); frame(ax, '(b) 추천 anode: 섬의 안쪽 부분만 식각')
anode(ax, ISL_NEW, ISL_GONE)
for b in (BAND_L, BAND_R):
    ax.add_patch(Rectangle((b[0], b[1]), b[2] - b[0], b[3] - b[1], fc='none', ec='#e67e22', lw=1.8, ls='-', zorder=8))
ax.text(5.3, 11.5, '식각\n창', ha='center', va='center', fontsize=9, color='#c0392b', fontweight='bold', zorder=9,
        bbox=dict(fc='white', ec='none', alpha=0.85, pad=1))
ax.text(2 * CX - 5.3, 11.5, '식각\n창', ha='center', va='center', fontsize=9, color='#c0392b', fontweight='bold', zorder=9,
        bbox=dict(fc='white', ec='none', alpha=0.85, pad=1))
box(ax, 0, H + 0.9, '남김: 가운데 anode 띠 + 섬의 바깥 컨택부(파랑)\n제거(빨간 빗금): 주황 창 안의 섬\n창: x 3.9–6.7 / 14.2–17.0, y 3.3–19.7 mm', '#a04000')

# (c) overlay with unchanged organic + cathode
ax = fig.add_subplot(gs[0, 2]); frame(ax, '(c) 추천 anode + 유기층·cathode (그대로)')
anode(ax, ISL_NEW); overlay(ax)
dimh(ax, 3.9, 5.93, 11.5); ax.text(4.9, 11.15, '2.0', ha='center', fontsize=8, color='#1e8449', fontweight='bold', zorder=9)
ax.add_patch(FancyArrowPatch((5.2, 3.3), (5.2, 4.83), arrowstyle='<->', mutation_scale=8, color='#1e8449', lw=1.2, zorder=9))
ax.text(3.0, 4.35, '1.5', fontsize=8, color='#1e8449', fontweight='bold', zorder=9)
box(ax, 0, H + 0.9, '유기층 아래에 cathode 전위 ITO 없음 → 누설 경로 소멸\n남은 섬 ↔ 유기층: x 2.0 mm, y 1.5 mm (정렬 1 mm 여유 OK)\ncathode 컨택은 기판 끝에서 IZO 위 Ag (그대로 튼튼)', '#1e8449')

# (d) how to etch
ax = fig.add_subplot(gs[0, 3]); frame(ax, '(d) 식각 방법 (테이프 배치)')
anode(ax, ISL_OLD)
tape = [R(-0.2, -0.2, 3.9, 23.2), R(6.7, -0.2, 2 * CX - 6.7, 23.2), R(2 * CX - 3.9, -0.2, 21.2, 23.2),
        R(3.9, -0.2, 6.7, 3.3), R(3.9, 19.7, 6.7, 23.2), R(2 * CX - 6.7, -0.2, 2 * CX - 3.9, 3.3), R(2 * CX - 6.7, 19.7, 2 * CX - 3.9, 23.2)]
for p in tape: ax.add_patch(Polygon(p, fc='#f5cba7', ec='#ca6f1e', lw=0.6, alpha=0.75, zorder=7))
ax.text(10.45, 11.5, 'Kapton\n(가운데 띠 보호)', ha='center', va='center', fontsize=8.5, color='#784212', fontweight='bold', zorder=9)
ax.text(1.95, 17.5, 'Kapton', ha='center', fontsize=8, color='#784212', fontweight='bold', rotation=90, zorder=9)
ax.text(2 * CX - 1.95, 17.5, 'Kapton', ha='center', fontsize=8, color='#784212', fontweight='bold', rotation=90, zorder=9)
box(ax, 0, H + 0.9, '① 뒷면 MLA/원뿔 면 전체 보호 필름\n② 앞면 Kapton: 노출 = 세로 창 2개뿐\n③ 묽은 옥살산(~3.4%) 또는 HCl 1:10, 짧게 반복\n④ 멀티미터: 띠 ↔ 섬 개방, 섬 바깥부 도통\n⑤ 아세톤/IPA/DI → 증착 직전 UV-오존', '#784212', size=8.3)

h = [Rectangle((0, 0), 1, 1, fc=C_AN, ec='#2471a3'), Rectangle((0, 0), 1, 1, fc=C_ISL, ec='#5dade2'),
     Rectangle((0, 0), 1, 1, fc='white', ec='#c0392b', ls='--', hatch='xxx'), Rectangle((0, 0), 1, 1, fc=C_ORG, ec='#c0392b', alpha=0.5),
     Rectangle((0, 0), 1, 1, fc=C_CAT, alpha=0.72), Rectangle((0, 0), 1, 1, fc='#f5cba7', ec='#ca6f1e')]
fig.legend(h, ['anode ITO (가운데 띠)', 'cathode 컨택용 ITO (남김)', '식각으로 제거', '유기층 (그대로)', 'cathode (그대로)', 'Kapton 테이프'],
           loc='upper center', ncol=6, fontsize=9.5, frameon=False, bbox_to_anchor=(0.5, 0.93))
fig.text(0.5, 0.115, '테이프 허용 오차 ±0.25 mm: 섬 쪽으로 남으면 고립 조각(무해), 띠 쪽으로 깎이면 픽셀 폭만 약간 감소(면적 실측).   '
         '기판 정렬(심)은 여전히 필요 — 손가락·C자 리드가 anode 띠에 닿는 단락은 섬과 무관.   치수: CAD 화면 판독(±0.1 mm), 단위 mm',
         ha='center', fontsize=8.6, color='#566573')
fig.savefig('anode_modification.png', dpi=160, bbox_inches='tight')
print('ok')
