// Two slides: screening descriptors (1) and per-material values (2).
//   NODE_PATH=<dir with pptxgenjs> node slides/make_seed_slides.js
// Data: runs/slide_descriptors.csv (scripts/152_slide_descriptors.py).
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");

const ROOT = path.join(__dirname, "..");
const OUT = path.join(__dirname, "seed_screening_2p.pptx");

// ---------- data ----------
const lines = fs.readFileSync(path.join(ROOT, "runs", "slide_descriptors.csv"), "utf8")
  .trim().split("\n");
const head = lines[0].split(",");
const rows = lines.slice(1).map((l) => {
  const v = l.split(",");
  const r = {};
  head.forEach((h, i) => (r[h] = i >= 2 && i <= 6 ? parseFloat(v[i]) : v[i]));
  return r;
});
const DISPLAY = { Cs2CO3: "Cs₂CO₃" };
rows.forEach((r) => (r.label = DISPLAY[r.name] || r.name));
const organics = rows.filter((r) => !r.name.startsWith("MoOx")); // + Cs2CO3 (inorganic EIL)
const moox = rows.find((r) => r.name.startsWith("MoOx"));
const klass = (r) =>
  r.name.startsWith("MoOx") ? "ref" : r.role.startsWith("HIL") ? "HIL"
    : (r.role === "host" || r.role === "HTL") ? "HOST" : "ETL";

// ---------- style ----------
const FONT = "Malgun Gothic";
const INK = "1F2937", SUB = "52514E", MUTED = "8A8F98";
const C = {
  elec: "1C5CAB", disp: "86B6EF",                 // one-hue ordinal pair (E_b parts)
  HIL: "EB6834", ETL: "1BAF7A", HOST: "6250D6",   // categorical slots, validated
  fail: "D03B3B",                                 // status: critical
  card: "F3F5F8", sub: "CDE2FB", ag: "C3C7CE", agLine: "6B7280", navy: "1E2761",
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.333 x 7.5 in
pres.theme = { headFontFace: FONT, bodyFontFace: FONT };

function title(slide, t, sub) {
  slide.addText(t, { x: 0.5, y: 0.28, w: 12.3, h: 0.6, fontFace: FONT, fontSize: 28,
    bold: true, color: INK, margin: 0 });
  slide.addText(sub, { x: 0.5, y: 0.88, w: 12.3, h: 0.4, fontFace: FONT, fontSize: 13,
    color: SUB, margin: 0 });
}

// ---------- slide 1: descriptors ----------
const s1 = pres.addSlide();
s1.background = { color: "FFFFFF" };
title(s1, "Thin Ag를 위한 seed layer 스크리닝 지표",
  "seed layer 특성에 영향을 줄 수 있는 지표를 가정 → 어떤 지표가 크게 작용하는지 계산으로 확인 → 최적 물질 선정");

const CW = 2.95, CH = 4.85, CY = 1.5, GAP = 0.2;
const cardX = (i) => 0.5 + i * (CW + GAP);

function substrate(slide, x, y, w, label) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h: 0.34, rectRadius: 0.06,
    fill: { color: C.sub }, line: { color: C.sub } });
  slide.addText(label, { x, y, w, h: 0.34, fontFace: FONT, fontSize: 9, color: SUB,
    align: "center", valign: "middle", margin: 0 });
}
function ag(slide, x, y, opts = {}) {
  const d = 0.46;
  slide.addShape(pres.shapes.OVAL, { x, y, w: d, h: d,
    fill: opts.ghost ? { color: "FFFFFF", transparency: 100 } : { color: C.ag },
    line: { color: C.agLine, width: 1, dashType: opts.ghost ? "dash" : "solid" } });
  slide.addText(opts.label || [{ text: "Ag" }], { x, y, w: d, h: d, fontFace: FONT,
    fontSize: 10, bold: true, color: opts.ghost ? MUTED : INK, align: "center",
    valign: "middle", margin: 0 });
}
function arrow(slide, x, y, w, h, opts = {}) {
  slide.addShape(pres.shapes.LINE, { x, y, w, h, flipV: !!opts.flipV,
    line: { color: opts.color || INK, width: 1.5, endArrowType: "triangle",
      beginArrowType: opts.both ? "triangle" : undefined } });
}
function note(slide, text, x, y, w, opts = {}) {
  slide.addText(text, { x, y, w, h: 0.3, fontFace: FONT, fontSize: opts.size || 10,
    bold: !!opts.bold, color: opts.color || INK, align: opts.align || "left",
    valign: "middle", margin: 0 });
}

const cards = [
  {
    name: "① E_b  흡착 에너지",
    def: "Ag 원자 1개가 분자의 가장 강한 자리에 붙는 에너지",
    why: "높을수록 Ag 원자가 재증발·확산하지 않고 붙어 핵이 촘촘히 생김 → 더 얇은 두께에서 연결",
    tag: "채택 · 순위 지표", tagColor: C.navy,
    draw: (x, y) => {
      substrate(s1, x + 0.25, y + 1.05, 2.45, "seed 분자");
      ag(s1, x + 1.2, y + 0.15);
      arrow(s1, x + 1.9, y + 0.3, 0, 0.7, { both: true });
      note(s1, "E_b", x + 2.0, y + 0.5, 0.6, { bold: true });
    },
  },
  {
    name: "② E_d  확산 장벽",
    def: "Ag가 표면에서 옆 자리로 옮겨 가는 데 드는 에너지",
    why: "낮으면 Ag가 멀리 이동해 큰 섬으로 뭉침 → 핵 밀도 감소",
    tag: "보조 · 판별력 낮음", tagColor: "9CA3AF",
    draw: (x, y) => {
      substrate(s1, x + 0.25, y + 1.05, 2.45, "seed 분자");
      ag(s1, x + 0.55, y + 0.5);
      ag(s1, x + 1.95, y + 0.5, { ghost: true });
      arrow(s1, x + 1.07, y + 0.73, 0.82, 0);
      note(s1, "E_d", x + 1.2, y + 0.18, 0.6, { bold: true, align: "center" });
    },
  },
  {
    name: "③ Ag 이온화  q(Ag₁)",
    def: "첫 Ag 원자가 기판에 넘겨준 전하",
    why: "결합 방식 구분: 크면 전하이동형 이온 결합, 0 근처면 배위 결합·물리흡착",
    tag: "채택 · 결합 성격 해석", tagColor: C.navy,
    draw: (x, y) => {
      substrate(s1, x + 0.25, y + 1.05, 2.45, "seed 분자");
      ag(s1, x + 1.0, y + 0.15, { label: [{ text: "Ag" }, { text: "δ+", options: { superscript: true } }] });
      arrow(s1, x + 1.23, y + 0.62, 0, 0.42, { color: C.HIL });
      note(s1, "e⁻", x + 1.35, y + 0.68, 0.4, { bold: true, color: C.HIL });
    },
  },
  {
    name: "④ 전하 되돌림  Δq",
    def: "Δq(Ag₁→Ag₂): 두 번째 Ag가 왔을 때 Ag 전체 전하의 변화",
    why: "되돌려 주면 Ag가 금속이 되어 섬끼리 합쳐짐. 계속 빼앗으면 Ag⁺로 고정 → 연결돼도 막이 닫히지 않음",
    tag: "채택 · 통과/탈락 기준", tagColor: C.navy,
    draw: (x, y) => {
      substrate(s1, x + 0.25, y + 1.05, 2.45, "seed 분자");
      ag(s1, x + 0.8, y + 0.15);
      ag(s1, x + 1.24, y + 0.15);
      arrow(s1, x + 1.25, y + 1.04, 0, -0.42, { color: C.ETL, flipV: true });
      note(s1, "e⁻ 되돌림", x + 1.4, y + 0.68, 1.1, { bold: true, color: C.ETL });
      note(s1, "Ag₂ 금속", x + 1.8, y + 0.22, 0.9, { size: 9, color: SUB });
    },
  },
];

cards.forEach((c, i) => {
  const x = cardX(i);
  s1.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: CY, w: CW, h: CH, rectRadius: 0.12,
    fill: { color: C.card }, line: { color: C.card } });
  c.draw(x, CY + 0.2);
  s1.addText(c.name, { x: x + 0.2, y: CY + 1.75, w: CW - 0.4, h: 0.45, fontFace: FONT,
    fontSize: 15, bold: true, color: INK, margin: 0 });
  s1.addText(c.def, { x: x + 0.2, y: CY + 2.22, w: CW - 0.4, h: 0.6, fontFace: FONT,
    fontSize: 11, color: SUB, margin: 0, valign: "top" });
  s1.addText([{ text: "왜 영향을 주나", options: { bold: true, breakLine: true } },
    { text: c.why }], { x: x + 0.2, y: CY + 2.85, w: CW - 0.4, h: 1.3, fontFace: FONT,
    fontSize: 11.5, color: INK, margin: 0, valign: "top", paraSpaceAfter: 3 });
  s1.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: x + 0.2, y: CY + CH - 0.6, w: CW - 0.4,
    h: 0.4, rectRadius: 0.2, fill: { color: c.tagColor }, line: { color: c.tagColor } });
  s1.addText(c.tag, { x: x + 0.2, y: CY + CH - 0.6, w: CW - 0.4, h: 0.4, fontFace: FONT,
    fontSize: 11, bold: true, color: "FFFFFF", align: "center", valign: "middle", margin: 0 });
});

// flow: gate -> rank -> result
const FY = 6.6, steps = [
  { t: "④ Δq로 거르기", s: "+0.3 e 이상 탈락", fill: "E5E7EB", ink: INK },
  { t: "① E_b로 순위", s: "통과한 물질끼리", fill: "E5E7EB", ink: INK },
  { t: "HATCN 1위 · MoOx 탈락", s: "→ SEM 비교", fill: C.navy, ink: "FFFFFF" },
];
steps.forEach((st, i) => {
  const w = i === 2 ? 4.3 : 3.6, x = 0.5 + [0, 4.05, 8.1][i];
  s1.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: FY, w, h: 0.6, rectRadius: 0.1,
    fill: { color: st.fill }, line: { color: st.fill } });
  s1.addText([{ text: st.t, options: { bold: true } }, { text: "   " + st.s,
    options: { fontSize: 10.5 } }], { x, y: FY, w, h: 0.6, fontFace: FONT, fontSize: 13,
    color: st.ink, align: "center", valign: "middle", margin: 0 });
  if (i < 2) arrow(s1, x + w + 0.07, FY + 0.3, 0.3, 0, { color: MUTED });
});
s1.addNotes("① E_b가 순위 변수, ④ Δq는 통과/탈락 게이트. ② E_d는 비정질 유기막에서 단일 분자 경로가 " +
  "대표성이 없고 계산값이 불안정해 보조로만 사용. 전자친화도(EA)·LUMO 깊이도 검토했으나 " +
  "F4TCNQ(EA 4.02 eV) < HATCN(3.38 eV) 결합 역전, MoOx는 더 깊지만 다공성 → 기각.");

// ---------- slide 2: values ----------
const s2 = pres.addSlide();
s2.background = { color: "FFFFFF" };
title(s2, "물질별 지표 계산값 (내림차순)",
  "유기물 19종 + Cs₂CO₃ · E_b, q(Ag₁): PBE0-D3/def2-SVP · Δq: GFN2-xTB · MoOx는 참고값");

const GX = [0.35, 4.67, 8.99], GW = 4.05, GY = 1.35, GH = 5.05;
const LAY = { x: 0.33, y: 0.1, w: 0.62, h: 0.76 }; // inner plot area, fraction of chart
const fmt2 = (v) => (v >= 0 ? v.toFixed(2) : v.toFixed(2));
const xAt = (i, v, lo, hi) => GX[i] + GW * (LAY.x + LAY.w * (v - lo) / (hi - lo));

function baseOpts(i, lo, hi, unit, chartTitle) {
  return {
    x: GX[i], y: GY, w: GW, h: GH, barDir: "bar", barGrouping: "stacked",
    catAxisLabelPos: "low",
    catAxisLabelFontSize: 9, catAxisLabelFontFace: FONT, catAxisLabelColor: INK,
    catAxisLineShow: false, catGridLine: { style: "none" },
    valAxisMinVal: lo, valAxisMaxVal: hi, valAxisMajorUnit: unit,
    valAxisLabelFontSize: 8.5, valAxisLabelFontFace: FONT, valAxisLabelColor: MUTED,
    valAxisLabelFormatCode: "0.0", valGridLine: { color: "E5E7EB", size: 0.5 },
    valAxisLineShow: false,
    showTitle: true, title: chartTitle, titleFontFace: FONT, titleFontSize: 12.5,
    titleColor: INK, barGapWidthPct: 45,
    showLegend: true, legendPos: "b", legendFontSize: 9, legendFontFace: FONT,
    legendColor: SUB, layout: LAY,
  };
}
function vline(i, v, lo, hi, color, label) {
  const x = xAt(i, v, lo, hi);
  s2.addShape(pres.shapes.LINE, { x, y: GY + GH * LAY.y, w: 0, h: GH * LAY.h,
    line: { color, width: 1.25, dashType: "dash" } });
  // label at the foot of the line, inside the plot, clear of the chart title
  s2.addText(label, { x: x + 0.05, y: GY + GH * (LAY.y + LAY.h) - 0.26, w: 1.2, h: 0.2,
    fontFace: FONT, fontSize: 8.5, bold: true, color, align: "left", margin: 0 });
}
function takeaway(i, text) {
  s2.addText(text, { x: GX[i] + 0.1, y: GY + GH + 0.08, w: GW - 0.2, h: 0.55, fontFace: FONT,
    fontSize: 10.5, color: INK, margin: 0, valign: "top" });
}

// chart 1: E_b, electronic + dispersion, sorted by total
{
  // ascending: PowerPoint draws the first category at the bottom
  const d = [...organics].sort((a, b) => a.E_b - b.E_b);
  const labels = d.map((r) => `${r.label}  ${r.E_b.toFixed(2)}`);
  s2.addChart(pres.charts.BAR, [
    { name: "화학 결합(전자적)", labels, values: d.map((r) => Math.max(0, r.E_elec)) },
    { name: "분산력", labels, values: d.map((r) => r.E_disp) },
  ], Object.assign(baseOpts(0, 0, 1.8, 0.3, "① E_b  흡착 에너지 (eV)"),
    { chartColors: [C.elec, C.disp] }));
  vline(0, 0.5, 0, 1.8, C.fail, "기준 0.5 eV");
  takeaway(0, "0.5 eV 통과: HIL 3종 + phen 킬레이트 2종 + Cs₂CO₃. PO-T2T는 대부분 분산력, BCP는 화학 결합");
}

// charts 2-3: one series per class, zeros hidden by the label format
function classSeries(d, key, withRef) {
  const labels = d.map((r) => r.label);
  const groups = [["HIL", "HIL"], ["ETL", "ETL · EIL"], ["HOST", "호스트 · HTL"]];
  if (withRef) groups.push(["ref", "MoOx (참고)"]);
  return groups.map(([k, nm]) => ({ name: nm, labels,
    values: d.map((r) => (klass(r) === k ? r[key] : 0)) }));
}
// one series per class, fully overlapped (not stacked) so labels can sit
// outside the bar end; the empty third section of the format hides zeros
const lblOpts = { barGrouping: "clustered", barOverlapPct: 100, showValue: true,
  dataLabelFormatCode: "+0.00;-0.00;", dataLabelFontSize: 8, dataLabelFontFace: FONT,
  dataLabelColor: INK, dataLabelPosition: "outEnd" };

{
  const d = [...organics].sort((a, b) => a.q_Ag1 - b.q_Ag1);
  s2.addChart(pres.charts.BAR, classSeries(d, "q_Ag1", false),
    Object.assign(baseOpts(1, -0.4, 0.8, 0.2, "③ Ag 이온화  q(Ag₁) (e)"), lblOpts,
      { chartColors: [C.HIL, C.ETL, C.HOST], valAxisLabelFormatCode: "+0.0;-0.0;0" }));
  takeaway(1, "Ag를 이온화하는 것은 HIL(받개)뿐 (+0.55~0.64). ETL은 배위(≈ −0.1), HTL·호스트는 물리흡착(≈ 0)");
}
{
  const d = [moox, ...organics].sort((a, b) => a.dq_Ag1_2 - b.dq_Ag1_2);
  s2.addChart(pres.charts.BAR, classSeries(d, "dq_Ag1_2", true),
    Object.assign(baseOpts(2, -0.8, 0.6, 0.2, "④ 전하 되돌림  Δq(Ag₁→Ag₂) (e)"), lblOpts,
      { chartColors: [C.HIL, C.ETL, C.HOST, C.fail], valAxisLabelFormatCode: "+0.0;-0.0;0" }));
  vline(2, 0.3, -0.8, 0.6, C.fail, "탈락 +0.3 e");
  takeaway(2, "유기물 19종과 Cs₂CO₃ 모두 통과(≤ 0: Ag₂ 금속). MoOx만 +0.46으로 탈락");
}

s2.addText("E_d(확산 장벽)는 비정질 유기막에서 계산 대표성이 낮아 그래프에서 제외. CBP·TCTA·TPD·TAPC는 mCP·NPB와 같은 π 계열(0.25–0.28 eV). " +
  "MoOx: E_b·q는 Mo₃O₉ 클러스터(PBE0, 2.16 eV / +0.72 e), Δq는 MoO₃(010) 표면 모델(GFN2). Ag가 수소에 막히는 자리는 PBE-D3로 위치를 다시 잡고 PBE0로 계산.",
  { x: 0.4, y: 7.0, w: 12.5, h: 0.4, fontFace: FONT, fontSize: 8.5, color: MUTED, margin: 0,
    valign: "top" });
s2.addNotes("각 그래프는 한 가지 계산법만 사용. E_b 막대는 화학 결합(진한 파랑)과 분산력(연한 파랑)을 쌓은 것. " +
  "Δq는 순위가 아니라 통과/탈락 기준이므로 유기물끼리의 차이는 강조하지 않음.");

pres.writeFile({ fileName: OUT }).then(() => console.log("wrote", OUT));
