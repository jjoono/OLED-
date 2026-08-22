const pptxgen = require("pptxgenjs");
const D = require("./data.json");

const NAVY="141B34", DEEP="1E2761", TEAL="2EC4B6", ORANGE="FF6B35",
      LIGHT="F7F9FC", MUTED="7A849C", SILVER="C5CEDE", WHITE="FFFFFF",
      INK="1A1F2E", GOLD="F2B134";
const KR="Malgun Gothic", EN="Calibri";

const p = new pptxgen();
p.layout = "LAYOUT_WIDE";           // 13.3 x 7.5
p.author = "Ag electrode study";
p.title  = "Ultrathin Ag top electrode";

const W=13.33, H=7.5;

function darkSlide(){ const s=p.addSlide(); s.background={color:NAVY}; return s; }
function lightSlide(){ const s=p.addSlide(); s.background={color:LIGHT}; return s; }

function title(s, txt, sub, dark){
  s.addText(txt, {x:0.62,y:0.42,w:12.1,h:0.72, fontFace:KR, fontSize:32, bold:true,
    color: dark?WHITE:INK, margin:0, valign:"middle"});
  if(sub) s.addText(sub, {x:0.62,y:1.14,w:12.1,h:0.42, fontFace:KR, fontSize:14,
    color: dark?SILVER:MUTED, margin:0, valign:"middle"});
}
function num(s, n, x, y, c){
  s.addShape(p.ShapeType.ellipse, {x:x,y:y,w:0.42,h:0.42, fill:{color:c}});
  s.addText(String(n), {x:x,y:y,w:0.42,h:0.42, fontFace:EN, fontSize:15, bold:true,
    color:WHITE, align:"center", valign:"middle", margin:0});
}
function card(s,x,y,w,h,fill){ s.addShape(p.ShapeType.roundRect,{x,y,w,h,rectRadius:0.06,
  fill:{color:fill||WHITE}, shadow:{type:"outer",color:"9AA5BC",blur:9,offset:1,angle:90,opacity:0.22}}); }
function stat(s,x,y,w,val,lab,c){
  s.addText(val,{x:x,y:y,w:w,h:0.78,fontFace:EN,fontSize:40,bold:true,color:c,align:"center",margin:0,valign:"middle"});
  s.addText(lab,{x:x,y:y+0.76,w:w,h:0.62,fontFace:KR,fontSize:11,color:MUTED,align:"center",margin:0,valign:"top"});
}
const AX = (extra)=>Object.assign({
  catAxisLabelColor:MUTED, valAxisLabelColor:MUTED,
  catAxisLabelFontSize:10, valAxisLabelFontSize:10,
  catAxisLabelFontFace:EN, valAxisLabelFontFace:EN,
  valGridLine:{color:"E2E7F0",size:1}, catGridLine:{style:"none"},
  showLegend:false, chartColors:[TEAL]
}, extra||{});

/* ---------------------------------------------------------------- 1 title */
{
  const s=darkSlide();
  s.addShape(p.ShapeType.ellipse,{x:9.6,y:-1.5,w:6.4,h:6.4,fill:{color:DEEP},line:{color:DEEP}});
  s.addShape(p.ShapeType.ellipse,{x:11.0,y:3.6,w:3.4,h:3.4,fill:{color:"22305F"},line:{color:"22305F"}});
  s.addText("초박막 Ag 상부전극의 흡수 저감", {x:0.75,y:2.25,w:8.9,h:1.0, fontFace:KR,
    fontSize:38, bold:true, color:WHITE, margin:0, valign:"middle"});
  s.addText("Seed layer engineering — 측정·시뮬레이션 종합", {x:0.75,y:3.32,w:8.9,h:0.5,
    fontFace:KR, fontSize:17, color:TEAL, margin:0, valign:"middle"});
  s.addText("HATCN vs MoOx · Ag 4–12 nm · 절대 T/R 분광 · 소자 광학 모델",
    {x:0.75,y:3.95,w:8.9,h:0.5, fontFace:KR, fontSize:12.5, color:SILVER, margin:0});
  s.addText("2026-08", {x:0.75,y:5.6,w:4,h:0.35, fontFace:EN, fontSize:12, color:MUTED, margin:0});
  s.addNotes("전극 흡수가 MLA 소자 효율을 지배한다는 것을 측정과 모델로 보인 자료.");
}

/* ------------------------------------------------------------ 2 the thesis */
{
  const s=lightSlide();
  title(s,"왜 '흡수'가 핵심 지표인가","MLA는 빠져나가지 못한 빛을 되돌린다 — 광자는 전극을 여러 번 통과한다");
  card(s,0.62,1.85,5.5,4.5);
  s.addText("한 번 통과 차이는 작다", {x:0.95,y:2.08,w:4.9,h:0.4,fontFace:KR,fontSize:15,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.bar,[{name:"one-pass A",labels:["측정 Ag 8nm","이상 Ag 8nm"],values:[4.74,0.94]}],
    AX({x:0.95,y:2.5,w:4.9,h:1.65, barDir:"bar", chartColors:[ORANGE,TEAL], showValue:true,
        dataLabelPosition:"outEnd", dataLabelColor:INK, dataLabelFontSize:11, dataLabelFontFace:EN,
        dataLabelFormatCode:'0.0"%"', valAxisMaxVal:6, catAxisLabelFontFace:KR, catAxisLabelFontSize:10}));
  s.addText("→ 겨우 3.8 %p 차이", {x:0.95,y:4.15,w:4.9,h:0.35,fontFace:KR,fontSize:12,color:MUTED,margin:0});
  s.addText("그런데 EQE는 18 %p 벌어진다", {x:0.95,y:4.62,w:4.9,h:0.4,fontFace:KR,fontSize:15,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.bar,[{name:"EQE",labels:["측정 Ag","이상 Ag"],values:[55,73]}],
    AX({x:0.95,y:5.02,w:4.9,h:1.1, barDir:"bar", chartColors:[ORANGE,TEAL], showValue:true,
        dataLabelPosition:"inEnd", dataLabelColor:WHITE, dataLabelFontSize:11, dataLabelFontFace:EN,
        dataLabelFormatCode:'0"%"', valAxisMaxVal:90, catAxisLabelFontFace:KR, catAxisLabelFontSize:10}));

  card(s,6.45,1.85,6.25,4.5);
  s.addText("증폭 기전 — 평균 통과 횟수", {x:6.8,y:2.08,w:5.6,h:0.4,fontFace:KR,fontSize:15,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.line,[
      {name:"측정 Ag",labels:D.lever_n.map(r=>""+r[0]),values:[]},
    ].slice(0,0).concat([
      {name:"측정",labels:["0.10","0.15","0.20","0.30","0.40","0.50"],values:[59.4,69.9,76.7,84.9,89.8,92.9]},
      {name:"이상",labels:["0.10","0.15","0.20","0.30","0.40","0.50"],values:[74.0,81.9,86.5,91.6,94.5,96.2]}]),
    AX({x:6.8,y:2.5,w:5.6,h:2.5, chartColors:[ORANGE,TEAL], lineSize:3, lineSmooth:true,
        showLegend:true, legendPos:"b", legendFontSize:10, legendFontFace:KR, legendColor:MUTED,
        valAxisMinVal:50, valAxisMaxVal:100, valAxisLabelFormatCode:'0"%"',
        catAxisTitle:"1회 탈출확률 q", catAxisTitleFontSize:10, catAxisTitleColor:MUTED, showCatAxisTitle:true}));
  s.addText([
    {text:"q ≈ 0.10 이면 광자는 평균 5.9회 통과", options:{bullet:true, breakLine:true}},
    {text:"매 통과마다 전극 흡수를 다시 낸다", options:{bullet:true, breakLine:true}},
    {text:"3.8 %p × 약 5회 → 15~20 %p의 EQE 손실", options:{bullet:true, bold:true, color:ORANGE}}],
    {x:6.8,y:5.1,w:5.6,h:1.05, fontFace:KR, fontSize:11.5, color:INK, margin:0, paraSpaceAfter:5});
  s.addNotes("핵심 논지 슬라이드. 단일 통과 흡수 차이가 MLA 재순환으로 5~6배 증폭된다.");
}

/* --------------------------------------------------- 3 absorption = resistivity */
{
  const s=lightSlide();
  title(s,"흡수와 면저항은 같은 물리량","A / A_bulk  =  ε₂ / ε₂_bulk  =  γ / γ_bulk  =  ρ / ρ_bulk");
  card(s,0.62,1.9,12.1,1.25,DEEP);
  s.addText("전자가 산란되는 정도가 저항을 만들고, 같은 산란이 빛을 흡수한다. 4-point probe가 곧 흡수계다.",
    {x:1.0,y:1.9,w:11.3,h:1.25, fontFace:KR, fontSize:15, color:WHITE, margin:0, valign:"middle"});

  card(s,0.62,3.35,6.0,3.0);
  s.addText("실측 검증 — HATCN / Ag", {x:0.95,y:3.55,w:5.4,h:0.35,fontFace:KR,fontSize:14,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.bar,[
    {name:"광학 ε₂/bulk",labels:["Ag 5","Ag 7","Ag 8"],values:[8.65,4.56,3.79]},
    {name:"전기 ρ/bulk",labels:["Ag 5","Ag 7","Ag 8"],values:[7.33,5.59,4.58]}],
    AX({x:0.95,y:3.95,w:5.4,h:2.25, chartColors:[TEAL,DEEP], showValue:true, dataLabelPosition:"outEnd",
        dataLabelColor:INK, dataLabelFontSize:9.5, dataLabelFontFace:EN, dataLabelFormatCode:'0.0"x"',
        showLegend:true, legendPos:"b", legendFontSize:9.5, legendFontFace:KR, legendColor:MUTED,
        valAxisMaxVal:11, catAxisLabelFontFace:EN}));

  card(s,6.9,3.35,5.82,3.0);
  s.addText("함의", {x:7.22,y:3.55,w:5.2,h:0.35,fontFace:KR,fontSize:14,bold:true,color:INK,margin:0});
  s.addText([
    {text:"두 독립 측정이 20 % 이내로 일치", options:{bullet:true, breakLine:true}},
    {text:"면저항만으로 흡수를 예측할 수 있다", options:{bullet:true, breakLine:true}},
    {text:"최적화 루프가 광학 측정 없이 돌아간다", options:{bullet:true, breakLine:true}},
    {text:"개선 여지 = Rs 비율 그 자체", options:{bullet:true, bold:true, color:DEEP}}],
    {x:7.22,y:3.98,w:5.2,h:1.5, fontFace:KR, fontSize:12, color:INK, margin:0, paraSpaceAfter:7});
  stat(s,7.22,5.35,2.5,"4.7×","Ag 8 nm — 이상적 벌크 대비",ORANGE);
  stat(s,9.9,5.35,2.5,"9.1","Ω/sq  현재 면저항",DEEP);
  s.addNotes("Rs가 흡수계라는 것이 이 프로젝트의 방법론적 토대.");
}

/* --------------------------------------------------------- 4 angle dependence */
{
  const s=lightSlide();
  title(s,"손실은 고각도에서 결정된다","MLA가 회수하려는 기판 모드가 바로 그 각도에 있다");
  card(s,0.62,1.9,7.4,4.5);
  s.addChart(p.ChartType.line,[
    {name:"측정 Ag 8nm",labels:D.ang.map(r=>r[0]+"°"),values:D.ang.map(r=>r[1])},
    {name:"이상 Ag 8nm",labels:D.ang.map(r=>r[0]+"°"),values:D.ang.map(r=>r[2])}],
    AX({x:0.9,y:2.2,w:6.85,h:3.9, chartColors:[ORANGE,TEAL], lineSize:3.5, lineSmooth:true,
        showLegend:true, legendPos:"b", legendFontSize:10.5, legendFontFace:KR, legendColor:MUTED,
        showValue:true, dataLabelPosition:"t", dataLabelColor:MUTED, dataLabelFontSize:9,
        dataLabelFontFace:EN, dataLabelFormatCode:'0.0',
        valAxisLabelFormatCode:'0"%"', catAxisTitle:"유기물 내부 각도", showCatAxisTitle:true,
        catAxisTitleFontSize:10, catAxisTitleColor:MUTED}));
  card(s,8.3,1.9,4.42,4.5);
  s.addText("읽을 것", {x:8.62,y:2.12,w:3.8,h:0.35,fontFace:KR,fontSize:14,bold:true,color:INK,margin:0});
  s.addText([
    {text:"수직 입사에서는 4.7 % — 크지 않다", options:{bullet:true, breakLine:true}},
    {text:"66° 에서 16.3 % 로 3.4배 급증", options:{bullet:true, breakLine:true}},
    {text:"유기물 1.8 → 기판 1.65 임계각이 66.5°", options:{bullet:true, breakLine:true}},
    {text:"이상적 Ag도 같은 배수로 증가하지만 절대량이 1/5", options:{bullet:true}}],
    {x:8.62,y:2.55,w:3.8,h:2.0, fontFace:KR, fontSize:11.5, color:INK, margin:0, paraSpaceAfter:8});
  s.addShape(p.ShapeType.roundRect,{x:8.62,y:4.75,w:3.8,h:1.4,rectRadius:0.06,fill:{color:"FFF1EC"}});
  s.addText("정면 흡수만 보면 문제를 놓친다.\n각도 적분이 필수.",
    {x:8.82,y:4.75,w:3.4,h:1.4, fontFace:KR, fontSize:12.5, bold:true, color:ORANGE, margin:0, valign:"middle"});
  s.addNotes("고각도 흡수가 MLA 효율을 지배한다.");
}


/* ------------------------------------------------------- 5 which metal */
{
  const s=lightSlide();
  title(s,"금속 선택 — Ag를 대체할 후보는 없다","모든 후보에 벌크 상수(최선의 경우)를 주고 비교. 여러분 Ag만 실측(열화된) 값");
  card(s,0.62,1.9,6.6,4.5);
  s.addText("ε₂ — 흡수를 곱하는 항 (550 nm)", {x:0.95,y:2.1,w:6.0,h:0.35,fontFace:KR,fontSize:14,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.bar,[{name:"eps2",
      labels:D.metals.map(r=>r[0]),
      values:D.metals.map(r=>+(2*r[1]*r[2]).toFixed(2))}],
    AX({x:0.95,y:2.5,w:6.0,h:3.65, barDir:"bar",
        chartColors:[TEAL,GOLD,ORANGE,DEEP,ORANGE,GOLD],
        showValue:true, dataLabelPosition:"outEnd", dataLabelColor:INK,
        dataLabelFontSize:10, dataLabelFontFace:EN, dataLabelFormatCode:'0.00',
        valAxisMaxVal:15, catAxisLabelFontFace:EN, catAxisLabelFontSize:10}));
  card(s,7.5,1.9,5.22,4.5);
  s.addText("소자 EQE (MLA 포함)", {x:7.82,y:2.1,w:4.6,h:0.35,fontFace:KR,fontSize:14,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.bar,[{name:"EQE",
      labels:D.metal_eqe.map(r=>r[0]),
      values:D.metal_eqe.map(r=>+r[1].toFixed(1))}],
    AX({x:7.82,y:2.5,w:4.6,h:2.6, barDir:"bar",
        chartColors:[ORANGE,TEAL,"C0392B",GOLD,"C0392B","C0392B"],
        showValue:true, dataLabelPosition:"inEnd", dataLabelColor:WHITE,
        dataLabelFontSize:10, dataLabelFontFace:EN, dataLabelFormatCode:'0"%"',
        valAxisMaxVal:85, catAxisLabelFontFace:EN, catAxisLabelFontSize:9.5}));
  s.addShape(p.ShapeType.roundRect,{x:7.82,y:5.2,w:4.6,h:1.0,rectRadius:0.06,fill:{color:"FFF1EC"}});
  s.addText("완벽히 평평한 Al 3 nm 도 거친 Ag 8 nm 보다 17 %p 나쁘다",
    {x:8.0,y:5.2,w:4.24,h:1.0, fontFace:KR, fontSize:12, bold:true, color:ORANGE, margin:0, valign:"middle"});
  s.addNotes("Al은 ε2가 Ag의 8배. wetting 이득으로 상쇄 불가. Au만 적색에서 검토 가치.");
}

/* ------------------------------------------------------ 6 measured dataset */
{
  const s=lightSlide();
  title(s,"확보한 실측 데이터","2 seed × 7 두께 = 15 시료, 절대 T/R 분광 + 4-point probe");
  card(s,0.62,1.85,5.9,4.55);
  s.addText("면저항 — 두 seed 비교", {x:0.95,y:2.05,w:5.3,h:0.35,fontFace:KR,fontSize:14,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.line,[
    {name:"HATCN 5nm",labels:D.rs.d.map(String),values:D.rs.HATCN},
    {name:"MoOx 5nm", labels:D.rs.d.map(String),values:D.rs.MoOx}],
    AX({x:0.95,y:2.45,w:5.3,h:3.7, chartColors:[TEAL,ORANGE], lineSize:3.5, lineSmooth:false,
        lineDataSymbol:"circle", lineDataSymbolSize:7,
        showLegend:true, legendPos:"b", legendFontSize:10, legendFontFace:EN, legendColor:MUTED,
        valAxisMaxVal:145, valAxisTitle:"Rs (Ω/sq)", showValAxisTitle:true,
        valAxisTitleFontSize:10, valAxisTitleColor:MUTED,
        catAxisTitle:"Ag 두께 (nm)", showCatAxisTitle:true,
        catAxisTitleFontSize:10, catAxisTitleColor:MUTED}));
  card(s,6.8,1.85,5.92,4.55);
  s.addText("데이터 구성", {x:7.12,y:2.05,w:5.3,h:0.35,fontFace:KR,fontSize:14,bold:true,color:INK,margin:0});
  const rows=[["절대 투과도 T","0° / 180°, 350–850 nm, 2 nm"],
              ["절대 반사도 R","6° / 12°, 뒷면 보정 84.3 %"],
              ["흡수도 A","1 − T − R, σ = 0.15 %p"],
              ["면저항","15 시료, 4-point probe"],
              ["추출 n, k","Ag 5·7·8 nm, 400–800 nm"],
              ["DFT 결합에너지","seed 후보 20종"],
              ["소자 광학 모델","TMM + MLA 재순환"]];
  rows.forEach((r,i)=>{
    const y=2.5+i*0.55;
    num(s,i+1,7.12,y,i<5?DEEP:TEAL);
    s.addText(r[0],{x:7.66,y:y-0.02,w:2.0,h:0.28,fontFace:KR,fontSize:11.5,bold:true,color:INK,margin:0});
    s.addText(r[1],{x:7.66,y:y+0.2,w:4.7,h:0.26,fontFace:EN,fontSize:9.5,color:MUTED,margin:0});
  });
  s.addNotes("측정 조건과 보정은 notes/HANDOFF_TR_20260820.md 에 전부 기록.");
}

/* --------------------------------------------------------- 7 A spectra */
{
  const s=lightSlide();
  title(s,"흡수 스펙트럼 — seed가 만드는 차이","glass / seed 5 nm / Ag, 절대 T·R 로부터 A = 1 − T − R");
  card(s,0.62,1.85,12.1,4.55);
  const st=D.spec, step=8;
  const idx=st.lam.map((_,i)=>i).filter(i=>i%step===0);
  s.addChart(p.ChartType.line,[
    {name:"HATCN / Ag 5",labels:idx.map(i=>st.lam[i].toFixed(0)),values:idx.map(i=>st.HATCN5_Ag5[i])},
    {name:"HATCN / Ag 8",labels:idx.map(i=>st.lam[i].toFixed(0)),values:idx.map(i=>st.HATCN5_Ag8[i])},
    {name:"MoOx / Ag 5", labels:idx.map(i=>st.lam[i].toFixed(0)),values:idx.map(i=>st.MoOx5_Ag5[i])},
    {name:"MoOx / Ag 8", labels:idx.map(i=>st.lam[i].toFixed(0)),values:idx.map(i=>st.MoOx5_Ag8[i])}],
    AX({x:0.95,y:2.2,w:11.45,h:4.0, chartColors:[TEAL,DEEP,ORANGE,"C0392B"],
        lineSize:3, lineSmooth:true,
        showLegend:true, legendPos:"b", legendFontSize:11, legendFontFace:EN, legendColor:MUTED,
        valAxisMaxVal:32, valAxisLabelFormatCode:'0"%"',
        valAxisTitle:"흡수도", showValAxisTitle:true, valAxisTitleFontSize:10, valAxisTitleColor:MUTED,
        catAxisTitle:"파장 (nm)", showCatAxisTitle:true, catAxisTitleFontSize:10, catAxisTitleColor:MUTED,
        catAxisLabelFrequency:5}));
  s.addNotes("HATCN 위 Ag는 전 파장에서 MoOx의 1/3 수준.");
}

/* ------------------------------------------------------------ 8 n,k / eps1 */
{
  const s=lightSlide();
  title(s,"추출된 광학상수 — 벌크 은에 근접","T·R 파장별 역산. k는 벌크와 거의 같고, 차이는 전부 n 에 있다");
  card(s,0.62,1.85,6.2,4.55);
  s.addText("n (550 nm 부근)", {x:0.95,y:2.05,w:5.6,h:0.32,fontFace:KR,fontSize:13.5,bold:true,color:INK,margin:0});
  const li=D.ag8.lam.map((_,i)=>i).filter(i=>i%25===0 && D.ag8.lam[i]<=700);
  s.addChart(p.ChartType.line,[
    {name:"Ag 5 nm",labels:li.map(i=>D.ag5.lam[i].toFixed(0)),values:li.map(i=>D.ag5.n[i])},
    {name:"Ag 8 nm",labels:li.map(i=>D.ag8.lam[i].toFixed(0)),values:li.map(i=>D.ag8.n[i])},
    {name:"McPeak 벌크",labels:li.map(i=>D.ag8.lam[i].toFixed(0)),
     values:li.map(i=>{const j=D.mcpeak.lam.indexOf(D.ag8.lam[i]); return j>=0?D.mcpeak.n[j]:null;})}],
    AX({x:0.95,y:2.42,w:5.6,h:1.85, chartColors:[GOLD,DEEP,TEAL], lineSize:3, lineSmooth:true,
        showLegend:true, legendPos:"b", legendFontSize:9.5, legendFontFace:EN, legendColor:MUTED,
        valAxisMaxVal:0.8, catAxisLabelFontSize:9}));
  s.addText("k (전 파장)", {x:0.95,y:4.35,w:5.6,h:0.32,fontFace:KR,fontSize:13.5,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.line,[
    {name:"Ag 8 nm",labels:li.map(i=>D.ag8.lam[i].toFixed(0)),values:li.map(i=>D.ag8.k[i])},
    {name:"McPeak 벌크",labels:li.map(i=>D.ag8.lam[i].toFixed(0)),
     values:li.map(i=>{const j=D.mcpeak.lam.indexOf(D.ag8.lam[i]); return j>=0?D.mcpeak.k[j]:null;})}],
    AX({x:0.95,y:4.7,w:5.6,h:1.5, chartColors:[DEEP,TEAL], lineSize:3, lineSmooth:true,
        showLegend:true, legendPos:"b", legendFontSize:9.5, legendFontFace:EN, legendColor:MUTED,
        catAxisLabelFontSize:9}));
  card(s,7.1,1.85,5.62,4.55);
  s.addText("ε₁ — 금속성의 지표", {x:7.42,y:2.05,w:5.0,h:0.32,fontFace:KR,fontSize:13.5,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.bar,[
    {name:"HATCN",labels:D.eps1.HATCN.map(r=>r[0]+" nm"),values:D.eps1.HATCN.map(r=>r[1])},
    {name:"MoOx", labels:D.eps1.HATCN.map(r=>r[0]+" nm"),
     values:D.eps1.HATCN.map(r=>{const m=D.eps1.MoOx.find(x=>x[0]===r[0]); return m?m[1]:null;})}],
    AX({x:7.42,y:2.42,w:5.0,h:2.9, chartColors:[TEAL,ORANGE],
        showLegend:true, legendPos:"b", legendFontSize:10, legendFontFace:EN, legendColor:MUTED,
        valAxisMinVal:-14, valAxisMaxVal:0, catAxisLabelFontSize:9.5, catAxisLabelFontFace:EN}));
  s.addShape(p.ShapeType.roundRect,{x:7.42,y:5.45,w:5.0,h:0.78,rectRadius:0.06,fill:{color:"E8F7F5"}});
  s.addText("벌크 은 ε₁ = "+D.bulk_eps1+"  ·  HATCN 위 Ag 는 4 nm 부터 오차 8 % 이내",
    {x:7.6,y:5.45,w:4.64,h:0.78, fontFace:KR, fontSize:11, bold:true, color:"0E8074", margin:0, valign:"middle"});
  s.addNotes("k는 벌크와 같고 n만 5~12배. 즉 벌크 밀도의 진짜 은인데 전자 산란만 많다.");
}

/* ------------------------------------------------------- 9 closure thickness */
{
  const s=lightSlide();
  title(s,"닫힘 두께 — 두 독립 측정이 일치한다","전기(면저항)와 광학(흡수)이 서로를 모른 채 같은 답을 준다");
  card(s,0.62,1.85,5.85,4.55);
  s.addText("① 전기 — 닫힌 막의 ρ = ρ₀ + C/d 선에서 벗어나는 지점",
    {x:0.95,y:2.05,w:5.25,h:0.5,fontFace:KR,fontSize:12,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.bar,[
    {name:"HATCN",labels:["4 nm","5 nm","6 nm","7 nm","8 nm"],values:[70,10,21,4,-8]},
    {name:"MoOx", labels:["4 nm","5 nm","6 nm","7 nm","8 nm"],values:[316,109,78,4,-9]}],
    AX({x:0.95,y:2.62,w:5.25,h:3.5, chartColors:[TEAL,ORANGE],
        showValue:true, dataLabelPosition:"outEnd", dataLabelColor:INK,
        dataLabelFontSize:9, dataLabelFontFace:EN, dataLabelFormatCode:'+0"%";-0"%"',
        showLegend:true, legendPos:"b", legendFontSize:10, legendFontFace:EN, legendColor:MUTED,
        valAxisMaxVal:360, valAxisTitle:"닫힌막 선 대비 초과 저항", showValAxisTitle:true,
        valAxisTitleFontSize:9.5, valAxisTitleColor:MUTED, catAxisLabelFontFace:EN}));
  card(s,6.75,1.85,5.97,4.55);
  s.addText("② 광학 — 흡수도의 계단",
    {x:7.07,y:2.05,w:5.35,h:0.32,fontFace:KR,fontSize:12,bold:true,color:INK,margin:0});
  s.addChart(p.ChartType.line,[
    {name:"HATCN",labels:D.a550.HATCN.map(r=>r[0]+""),values:D.a550.HATCN.map(r=>r[1])},
    {name:"MoOx", labels:D.a550.HATCN.map(r=>r[0]+""),
     values:D.a550.HATCN.map(r=>{const m=D.a550.MoOx.find(x=>x[0]===r[0]); return m?m[1]:null;})}],
    AX({x:7.07,y:2.42,w:5.35,h:2.55, chartColors:[TEAL,ORANGE], lineSize:3.5,
        lineDataSymbol:"circle", lineDataSymbolSize:7, lineSmooth:false,
        showLegend:true, legendPos:"b", legendFontSize:10, legendFontFace:EN, legendColor:MUTED,
        valAxisMaxVal:32, valAxisLabelFormatCode:'0"%"',
        catAxisTitle:"Ag 두께 (nm)", showCatAxisTitle:true,
        catAxisTitleFontSize:9.5, catAxisTitleColor:MUTED}));
  s.addText("MoOx 는 4–6 nm 에서 28 % 평평 → 7 nm 에서 23 % 로 계단.  HATCN 은 4 nm 부터 단조 감소.\n※ MoOx 12 nm 은 반사도 파일 결측.",
    {x:7.07,y:5.05,w:5.35,h:0.5, fontFace:KR, fontSize:10.5, color:MUTED, margin:0});
  s.addShape(p.ShapeType.roundRect,{x:7.07,y:5.6,w:2.55,h:0.75,rectRadius:0.06,fill:{color:"E8F7F5"}});
  s.addText("HATCN  5 nm",{x:7.07,y:5.6,w:2.55,h:0.75,fontFace:KR,fontSize:15,bold:true,
    color:"0E8074",align:"center",valign:"middle",margin:0});
  s.addShape(p.ShapeType.roundRect,{x:9.87,y:5.6,w:2.55,h:0.75,rectRadius:0.06,fill:{color:"FFF1EC"}});
  s.addText("MoOx  7 nm",{x:9.87,y:5.6,w:2.55,h:0.75,fontFace:KR,fontSize:15,bold:true,
    color:ORANGE,align:"center",valign:"middle",margin:0});
  s.addNotes("HATCN이 Ag 2 nm를 절약한다. 전기·광학 독립 확인.");
}

/* ------------------------------------------------- 10 device A vs thickness */
{
  const s=lightSlide();
  title(s,"소자 흡수 — 최적 두께와 남은 여지","organic 1.8 / Ag / CPL 2.1 / air, CPL 두께는 매번 재최적화");
  card(s,0.62,1.85,7.9,4.55);
  s.addChart(p.ChartType.line,[
    {name:"HATCN seed",labels:D.adev.HATCN.map(r=>r[0]+""),values:D.adev.HATCN.map(r=>r[1])},
    {name:"MoOx seed", labels:D.adev.HATCN.map(r=>r[0]+""),
     values:D.adev.HATCN.map(r=>{const m=D.adev.MoOx.find(x=>x[0]===r[0]); return m?m[1]:null;})},
    {name:"이상적 벌크 Ag", labels:D.adev.HATCN.map(r=>r[0]+""),
     values:D.adev.HATCN.map(r=>{const m=D.ideal.find(x=>x[0]===r[0]); return m?m[1]:null;})}],
    AX({x:0.95,y:2.2,w:7.25,h:4.0, chartColors:[TEAL,ORANGE,DEEP], lineSize:3.5,
        lineDataSymbol:"circle", lineDataSymbolSize:7, lineSmooth:false,
        showLegend:true, legendPos:"b", legendFontSize:10.5, legendFontFace:KR, legendColor:MUTED,
        valAxisMaxVal:12, valAxisLabelFormatCode:'0"%"',
        valAxisTitle:"소자 one-pass 흡수", showValAxisTitle:true,
        valAxisTitleFontSize:10, valAxisTitleColor:MUTED,
        catAxisTitle:"Ag 두께 (nm)", showCatAxisTitle:true,
        catAxisTitleFontSize:10, catAxisTitleColor:MUTED}));
  card(s,8.8,1.85,3.92,4.55);
  s.addText("핵심", {x:9.1,y:2.05,w:3.4,h:0.32,fontFace:KR,fontSize:14,bold:true,color:INK,margin:0});
  stat(s,9.1,2.45,3.4,"7–8 nm","HATCN 최적 두께",TEAL);
  stat(s,9.1,3.75,3.4,"2.6 %","현재 소자 흡수",DEEP);
  stat(s,9.1,5.05,3.4,"4.7×","이상적 하한 대비 여지",ORANGE);
  s.addNotes("이상적 Ag는 두께에 선형. 실제 막은 7~8 nm에 최소점.");
}

/* --------------------------------------------------------- 11 roadmap to 70 */
{
  const s=darkSlide();
  title(s,"70 % 달성 경로","단일 조치로는 도달하지 못한다 — 기판·AR·Ag 품질을 병행해야 한다",true);
  const rm=D.roadmap;
  s.addChart(p.ChartType.bar,[{name:"EQE",labels:rm.map(r=>r[0]),values:rm.map(r=>r[1])}],
    {x:0.62,y:1.95,w:7.6,h:4.5, barDir:"bar",
     chartColors:[ORANGE,"D98032",GOLD,"8FBF6E",TEAL,"1E9E93"],
     showValue:true, dataLabelPosition:"inEnd", dataLabelColor:WHITE,
     dataLabelFontSize:12, dataLabelFontFace:EN, dataLabelFormatCode:'0.0"%"',
     showLegend:false, valAxisMaxVal:85,
     catAxisLabelColor:SILVER, valAxisLabelColor:SILVER,
     catAxisLabelFontSize:11, valAxisLabelFontSize:10,
     catAxisLabelFontFace:KR, valAxisLabelFontFace:EN,
     valGridLine:{color:"2A3358",size:1}, catGridLine:{style:"none"}});
  s.addShape(p.ShapeType.roundRect,{x:8.55,y:1.95,w:4.17,h:2.1,rectRadius:0.07,fill:{color:DEEP}});
  s.addText("지렛대 크기", {x:8.85,y:2.15,w:3.6,h:0.3,fontFace:KR,fontSize:13,bold:true,color:WHITE,margin:0});
  s.addText([
    {text:"기판 1.65 → 1.80      + 6.7 %p", options:{bullet:true, breakLine:true}},
    {text:"MLA 반사방지         + 1.3 %p", options:{bullet:true, breakLine:true}},
    {text:"Ag Rs 9.1 → 5.3      + 7.3 %p", options:{bullet:true, breakLine:true}},
    {text:"CPL 굴절률           ± 1.5 %p", options:{bullet:true}}],
    {x:8.85,y:2.5,w:3.6,h:1.4, fontFace:KR, fontSize:11, color:SILVER, margin:0, paraSpaceAfter:6});
  s.addShape(p.ShapeType.roundRect,{x:8.55,y:4.25,w:4.17,h:2.2,rectRadius:0.07,fill:{color:"1B2547"}});
  s.addText("병행의 효과", {x:8.85,y:4.45,w:3.6,h:0.3,fontFace:KR,fontSize:13,bold:true,color:TEAL,margin:0});
  s.addText("기판과 AR을 먼저 확보하면 Ag 목표가 완화된다.\n\nRs 4.4 → 5.3 Ω/sq\n개선 부담 2.1배 → 1.7배",
    {x:8.85,y:4.8,w:3.6,h:1.5, fontFace:KR, fontSize:11.5, color:WHITE, margin:0});
  s.addNotes("여러분 MATLAB PSO 두 점(55%, 73%)에 보정한 재순환 모델 기준.");
}

/* ------------------------------------------- 12 seed engineering direction */
{
  const s=lightSlide();
  title(s,"Seed layer engineering — 다음 단계","목표는 닫힘 두께가 아니라 ρ, 즉 결정립 크기와 계면 경면성");
  const cards=[
    ["1","확인된 것","HATCN 이 MoOx 대비 닫힘을 2 nm 앞당기고 흡수를 1/3 로 낮춘다. DFT 결합에너지 1.03 vs 0.45 eV 와 순서가 일치.",TEAL],
    ["2","진짜 지렛대","닫힘 이후에도 ρ 가 벌크의 4.6배. 표면 산란 3.88, 결정립계 1.81 µΩ·cm. 두 항을 모두 줄여야 한다.",DEEP],
    ["3","공정 변수","증착 속도 상향(침투 억제), 기판 냉각(결정립 성장), seed 두께 5 nm 고정, 즉시 캡핑.",GOLD],
    ["4","다음 측정","검출기 각스캔으로 산란 분리, 엘립소·XRR 로 두께 확정, Ag 2·3 nm 로 퍼콜레이션 하한.",ORANGE]];
  cards.forEach((c,i)=>{
    const x=0.62+i*3.06, w=2.88;
    card(s,x,1.9,w,3.5);
    num(s,c[0],x+0.28,2.15,c[3]);
    s.addText(c[1],{x:x+0.28,y:2.72,w:w-0.56,h:0.35,fontFace:KR,fontSize:13.5,bold:true,color:INK,margin:0});
    s.addText(c[2],{x:x+0.28,y:3.12,w:w-0.56,h:2.1,fontFace:KR,fontSize:10.5,color:MUTED,margin:0});
  });
  card(s,0.62,5.62,12.1,1.25,DEEP);
  s.addText("Rs 는 흡수계다 — 광학 측정 없이 4-point probe 만으로 최적화 루프를 돌릴 수 있다",
    {x:1.0,y:5.62,w:11.3,h:1.25, fontFace:KR, fontSize:15, bold:true, color:WHITE, margin:0, valign:"middle"});
  s.addNotes("결론 슬라이드. 닫힘이 아니라 미세구조가 목표라는 전환.");
}
p.writeFile({fileName:"Ag_electrode_summary.pptx"}).then(f=>console.log("ok",f));
