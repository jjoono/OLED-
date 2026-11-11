// Two slides: design pipeline + method finding (1), candidates vs formulations (2).
//   NODE_PATH=<dir with pptxgenjs> node slides/make_design_slides.js
// Data: runs/design/*.json (scripts 160-167), notes/DESIGN_RESULTS.md.
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");

const D = path.join(__dirname, "..", "runs", "design");
const OUT = path.join(__dirname, "seed_design_2p.pptx");
const J = (f) => JSON.parse(fs.readFileSync(path.join(D, f), "utf8"));
const JL = (f) => Object.fromEntries(fs.readFileSync(path.join(D, f), "utf8").trim().split("\n")
  .map((l) => JSON.parse(l)).map((r) => [r.name, r]));

// ---------- data ----------
const anchors = J("pbe0_anchor.json");
const sites = JL("sites.jsonl");
const props = JL("props.jsonl");
const td = J("tddft.json");
const form = J("formulations.json").rows;
const ANCH = ["HAT-CN-CN-F2benzo", "HATCN", "HAT-CN-CN-CN2benzo",
  "HAT-CN2benzo-CN2benzo-CN2benzo", "HAT-F4benzo-F4benzo-F4benzo"];
const SHORT = {
  "HAT-CN-CN-F2benzo": "CN-CN-F₂benzo", HATCN: "HATCN", "HAT-CN-CN-CN2benzo": "CN-CN-(CN)₂benzo",
  "HAT-CN2benzo-CN2benzo-CN2benzo": "(CN)₂benzo ×3", "HAT-F4benzo-F4benzo-F4benzo": "F₄benzo ×3",
};
const A = ANCH.map((n) => ({
  n, pbe0: anchors[n].result.E_b, gfn2: sites[n].E_b_gfn2, ea: props[n].EA,
}));
const pearson = (x, y) => {
  const mx = x.reduce((a, b) => a + b) / x.length, my = y.reduce((a, b) => a + b) / y.length;
  let sxy = 0, sxx = 0, syy = 0;
  x.forEach((xi, i) => { sxy += (xi - mx) * (y[i] - my); sxx += (xi - mx) ** 2; syy += (y[i] - my) ** 2; });
  return sxy / Math.sqrt(sxx * syy);
};
const rG = pearson(A.map((a) => a.gfn2), A.map((a) => a.pbe0));
const rE = pearson(A.map((a) => a.ea), A.map((a) => a.pbe0));

// ---------- style ----------
const FONT = "Malgun Gothic";
const INK = "1F2937", SUB = "52514E", MUTED = "8A8F98";
const C = {
  accent: "1C5CAB", light: "CDE2FB", card: "F3F5F8", navy: "1E2761",
  good: "1BAF7A", warn: "E59A1C", bad: "D03B3B", grey: "B8BEC8",
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.333 x 7.5 in
pres.theme = { headFontFace: FONT, bodyFontFace: FONT };

function title(slide, t, sub) {
  slide.addText(t, { x: 0.5, y: 0.28, w: 12.3, h: 0.6, fontFace: FONT, fontSize: 26,
    bold: true, color: INK, margin: 0 });
  slide.addText(sub, { x: 0.5, y: 0.88, w: 12.3, h: 0.4, fontFace: FONT, fontSize: 13,
    color: SUB, margin: 0 });
}

function banner(slide, y, text) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.5, y, w: 12.33, h: 0.62, rectRadius: 0.08,
    fill: { color: C.navy }, line: { color: C.navy } });
  slide.addText(text, { x: 0.75, y, w: 11.9, h: 0.62, fontFace: FONT, fontSize: 14,
    color: "FFFFFF", valign: "middle", margin: 0 });
}

// ---------- slide 1: pipeline + method finding ----------
const s1 = pres.addSlide();
s1.background = { color: "FFFFFF" };
title(s1, "AI 기반 seed 분자 설계: 파이프라인과 핵심 발견",
  "역할 고정 (n-ETL / seed / Ag 의 CGL p-side) → 라이브러리 → 저가 계산 → 계열별 PBE0 보정 → 파레토 → TD-DFT");

// funnel (left)
const STAGES = [
  ["701", "후보 라이브러리", "HAT · dppz · Bphen 유도체, 합성 경로 표시"],
  ["225", "저가 물성 (xTB)", "EA · IP · HOMO–LUMO 갭 · 결정화 위험"],
  ["83", "Ag 사이트 탐색 (GFN2)", "결합 사이트 · E_b · Ag₂ 전하 변화 Δq"],
  ["5", "PBE0 앵커 → 계열 보정", "불확실한 영역만 고정밀 계산"],
  ["17→3", "조건 통과 → TD-DFT", "중성 + 음이온 흡수"],
];
const FX = 0.5, FY = 1.55, FW = 5.6, FH = 0.78, FG = 0.14;
STAGES.forEach(([num, head, desc], i) => {
  const inset = i * 0.22;
  const y = FY + i * (FH + FG);
  s1.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: FX + inset, y, w: FW - 2 * inset, h: FH,
    rectRadius: 0.08, fill: { color: i === 3 ? C.accent : C.card },
    line: { color: i === 3 ? C.accent : "DDE1E7", width: 0.75 } });
  const ink = i === 3 ? "FFFFFF" : INK;
  s1.addText(num, { x: FX + inset + 0.15, y, w: 1.15, h: FH, fontFace: FONT, fontSize: 22,
    bold: true, color: i === 3 ? "FFFFFF" : C.accent, valign: "middle", margin: 0 });
  s1.addText([
    { text: head, options: { bold: true, fontSize: 13, color: ink, breakLine: true } },
    { text: desc, options: { fontSize: 10.5, color: i === 3 ? "E6EEF9" : SUB } },
  ], { x: FX + inset + 1.3, y, w: FW - 2 * inset - 1.4, h: FH, fontFace: FONT,
    valign: "middle", margin: 0 });
});

// two scatter charts (right)
function scatter(slide, x, y, xs, ys, xTitle, chartTitle, r, color, xmin, xmax) {
  slide.addChart(pres.charts.SCATTER, [
    { name: "X", values: xs },
    { name: "PBE0 E_b", values: ys },
  ], {
    x, y, w: 3.08, h: 3.0,
    showTitle: true, title: chartTitle, titleFontFace: FONT, titleFontSize: 12, titleColor: INK,
    lineSize: 0, lineDataSymbol: "circle", lineDataSymbolSize: 10,
    chartColors: [color], showLegend: false,
    valAxisMinVal: 1.0, valAxisMaxVal: 1.5, valAxisMajorUnit: 0.1,
    valAxisTitle: "PBE0 E_b (eV)", showValAxisTitle: true, valAxisTitleFontSize: 10,
    valAxisTitleColor: SUB, valAxisLabelFontSize: 9, valAxisLabelColor: MUTED,
    valAxisLabelFormatCode: "0.0",
    catAxisTitle: xTitle, showCatAxisTitle: true, catAxisTitleFontSize: 10, catAxisTitleColor: SUB,
    catAxisLabelFontSize: 9, catAxisLabelColor: MUTED, catAxisMinVal: xmin, catAxisMaxVal: xmax,
    catAxisLabelFormatCode: "0.0",
    valGridLine: { color: "E5E7EB", size: 0.5 }, catGridLine: { style: "none" },
  });
  slide.addText(`r = ${r >= 0 ? "+" : "−"}${Math.abs(r).toFixed(2)}`, { x: x + 1.98, y: y + 0.45,
    w: 1.0, h: 0.32, fontFace: FONT, fontSize: 14, bold: true, color, margin: 0, align: "right" });
}
const SX = 6.55, SY = 1.5;
scatter(s1, SX, SY, A.map((a) => a.gfn2), A.map((a) => a.pbe0), "GFN2 E_b (eV)",
  "GFN2 결합 에너지 → 순서 맞음", rG, C.accent, 2.9, 3.6);
scatter(s1, SX + 3.2, SY, A.map((a) => a.ea), A.map((a) => a.pbe0), "전자친화도 EA (eV)",
  "EA → 예측 못 함", rE, C.bad, 3.2, 3.9);
s1.addText("같은 HAT 계열 PBE0 앵커 5개 (같은 계산 방식). EA가 가장 낮은 F₄benzo×3가 Ag를 가장 강하게 붙잡음 → Ag 결합은 분자 전체 EA가 아니라 붙는 자리의 구조가 결정",
  { x: SX, y: SY + 3.08, w: 6.28, h: 0.62, fontFace: FONT, fontSize: 10.5, color: SUB, margin: 0 });

// callout: cross-family failure
s1.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: SX, y: 5.25, w: 6.28, h: 1.0, rectRadius: 0.08,
  fill: { color: "FDECEC" }, line: { color: "F5C2C2", width: 0.75 } });
s1.addText([
  { text: "계열이 다르면 GFN2 순위가 뒤집힘", options: { bold: true, color: C.bad, breakLine: true } },
  { text: "GFN2: F6-TCNNQ 3.52 > HATCN 3.00 > F4TCNQ 2.86 eV   ↔   PBE0: HATCN 1.63 > F4TCNQ 1.21 ≈ F6-TCNNQ 1.20 eV",
    options: { color: INK } },
], { x: SX + 0.2, y: 5.25, w: 5.9, h: 1.0, fontFace: FONT, fontSize: 11, valign: "middle", margin: 0 });

banner(s1, 6.45, "도구 설계 원칙: 계열마다 PBE0 앵커 ~5개로 전용 보정식 (교차검증 오차 0.18 eV) — 전체 공통 보정은 강한 억셉터에서 실패");

// ---------- slide 2: candidates vs formulations ----------
const s2 = pres.addSlide();
s2.background = { color: "FFFFFF" };
title(s2, "결과: 분자 하나로는 HATCN을 넘기 어렵고, 결정화는 조성으로 푼다",
  "Ag 결합(E_b)을 올리는 축합 고리가 흡수를 가시광 쪽으로 끌어옴 — HATCN은 이 상충 관계의 투명한 끝");

// candidate table (left)
const tdN = (n) => td[n].neutral.states[0];
const g = (t) => ({ text: t, options: { color: C.good, bold: true } });
const b = (t) => ({ text: t, options: { color: C.bad, bold: true } });
const w = (t) => ({ text: t, options: { color: C.warn, bold: true } });
const p = (t) => ({ text: t, options: { color: INK } });
const hdr = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: C.navy }, align: "center" } });
const cn2 = "HAT-CN2benzo-CN2benzo-CN2benzo", f4 = "HAT-CN-F4benzo-F4benzo";
const rowsT = [
  [hdr("분자"), hdr("E_b (eV)"), hdr("중성 최저 들뜸"), hdr("음이온 가시광 흡수"), hdr("결정화 위험"), hdr("합성")],
  [p("HATCN"), p("1.63"), g(`${tdN("HATCN").eV.toFixed(2)} eV, 흡수 없음`), g("720 nm, f 0.16"), b("92"), g("원팟")],
  [p("(CN)₂benzo ×3"), g("1.74 (PBE0)"), w(`${tdN(cn2).eV.toFixed(2)} eV, f 0.11`), b("697 nm, f 0.48 (3배)"), b("88"), g("원팟")],
  [p("CN-F₄benzo-F₄benzo"), g("1.89 (예측)"), w(`${tdN(f4).eV.toFixed(2)} eV (391 nm)`), b("545 nm 띠 (초록)"), w("69"), w("분리 필요")],
  [p("CNPh ×3 (아릴 HAT)"), p("~1.45 (외삽)"), p("갭 넓음 (GFN2)"), p("–"), w("59"), b("EA 2.96 미달")],
];
s2.addTable(rowsT, { x: 0.5, y: 1.55, w: 7.1, colW: [1.75, 1.05, 1.35, 1.35, 0.75, 0.85],
  fontFace: FONT, fontSize: 10, valign: "middle", rowH: 0.52,
  border: { type: "solid", pt: 0.5, color: "DDE1E7" }, fill: { color: "FFFFFF" } });
s2.addText([
  { text: "결합·EA 개선 후보 17개 존재 (E_b 1.65–1.89 eV), ", options: {} },
  { text: "그러나 TD-DFT 최종 3개 중 HATCN만큼 투명한 분자 없음", options: { bold: true } },
  { text: " · 축합 HAT는 여전히 평평한 원반 (결정화 위험 92 → 67–72에 그침)", options: {} },
], { x: 0.5, y: 4.3, w: 7.1, h: 0.85, fontFace: FONT, fontSize: 11, color: SUB, margin: 0,
  valign: "top" });
s2.addText("E_b: 기존 PBE0 기준 환산 (Ag 위치 보정 +0.42 eV 가정). TD-DFT: PBE0-TDA, 중성 3·음이온 4개 상태.",
  { x: 0.5, y: 5.2, w: 7.1, h: 0.4, fontFace: FONT, fontSize: 9, color: MUTED, margin: 0 });

// formulation bar chart (right): stacked by verdict so each bar takes its class colour
const SHOW = ["HATCN neat", "HATCN:p-bPPhenB 1:1", "HATCN:TPBi 1:1", "HATCN:PO-T2T 1:1",
  "HATCN:DPEPO 1:1", "HATCN:NPB 1:1", "HATCN:TAPC 1:1"];
const LAB = {
  "HATCN neat": "HATCN 단독", "HATCN:p-bPPhenB 1:1": "+ p-bPPhenB", "HATCN:TPBi 1:1": "+ TPBi",
  "HATCN:PO-T2T 1:1": "+ PO-T2T", "HATCN:DPEPO 1:1": "+ DPEPO", "HATCN:NPB 1:1": "+ NPB",
  "HATCN:TAPC 1:1": "+ TAPC",
};
const VERDICT = {
  "HATCN neat": "bad", "HATCN:p-bPPhenB 1:1": "warn", "HATCN:TPBi 1:1": "good",
  "HATCN:PO-T2T 1:1": "good", "HATCN:DPEPO 1:1": "good", "HATCN:NPB 1:1": "bad",
  "HATCN:TAPC 1:1": "bad",
};
const byName = Object.fromEntries(form.map((r) => [r.formulation, r]));
const cats = SHOW.map((s) => LAB[s]);
const ser = (cls, name) => ({ name, labels: cats,
  values: SHOW.map((s) => (VERDICT[s] === cls ? byName[s].N_sat_rel : 0)) });
s2.addChart(pres.charts.BAR, [ser("good", "권장"), ser("warn", "주의"), ser("bad", "제외")], {
  x: 7.9, y: 1.45, w: 4.95, h: 3.75, barDir: "bar", barGrouping: "stacked",
  chartColors: [C.good, C.warn, C.bad],
  showTitle: true, title: "HATCN 1:1 혼합의 Ag 섬 밀도 (단독 = 1)",
  titleFontFace: FONT, titleFontSize: 12, titleColor: INK,
  showLegend: true, legendPos: "b", legendFontFace: FONT, legendFontSize: 9,
  catAxisOrientation: "maxMin", catAxisLabelFontFace: FONT, catAxisLabelFontSize: 10,
  catAxisLabelColor: INK, catGridLine: { style: "none" },
  valAxisMinVal: 0.7, valAxisMaxVal: 1.0, valAxisMajorUnit: 0.1, valAxisLabelFormatCode: "0.0",
  valAxisLabelFontSize: 9, valAxisLabelColor: MUTED, valGridLine: { color: "E5E7EB", size: 0.5 },
  showValue: false,
});
s2.addText([
  { text: "HATCN 단독: ", options: { bold: true, color: C.bad } },
  { text: "결정 영역이 막 전체로 연결", options: { breakLine: true } },
  { text: "NPB · TAPC: ", options: { bold: true, color: C.bad } },
  { text: "HATCN과 전하이동 흡수", options: { breakLine: true } },
  { text: "TPBi: ", options: { bold: true, color: C.good } },
  { text: "E_b 0.46 eV로 포획 확률 ~50%이지만 결정화 위험 17", options: {} },
], { x: 7.95, y: 5.25, w: 4.9, h: 0.95, fontFace: FONT, fontSize: 10, color: SUB, margin: 0,
  valign: "top" });

banner(s2, 6.45, "권장: HATCN:TPBi 1:1 (상대 물질 45–60%) · 다음 설계 목표: 분자량 < 900의 비틀린 아릴 HAT 상대 물질");

pres.writeFile({ fileName: OUT }).then(() => console.log("written", OUT));
