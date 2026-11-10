# HANDOFF — Ag seed layer / thin-Ag electrode 프로젝트 인수인계

작성: 2026-07-06 (클라우드 세션 이관용)
사용자: jhk4733@gmail.com — top-emitting OLED, ultrahigh EQE 목표.
목표: 6~15 nm 열증착 Ag를 island 없이 균일하게 올릴 seed layer를 DFT로 사전 선별 + 광학(투과/반사/흡수) 정량.

---

## 0. 30초 요약 (결론만)

- **Ag wetting 서열(DFT 결합에너지 E_b)**: 청정 Al 2.24 ≫ **HATCN 1.03** > TPBi 0.89 ≈ p-bPPhenB 0.87 > Cs₂CO₃ 0.90(클러스터, 상한) > B3PyMPM 0.63 > Bphen 0.49 > MoOx결손 0.44 ≈ **AlOx 0.42** > MoO₃ 0.28 ≈ LiF 0.25 ≈ TAPC 0.25 ≈ TCTA 0.28 > Liq 0.17.
- **핵심 반전**: 청정 Al은 최강이지만 실공정에서 **AlOx로 산화(0.42)** → 유기물보다 나쁨. SEM에서 Al/Ag 다공성 vs HATCN/Ag 연속막으로 실증됨.
- **확산장벽(#2, DFT)**: HATCN 0.29 > AlOx 0.17 > TAPC ~0 ≈ 청정Al ~0. (TPBi/B3PyMPM은 클러스터 한계로 미확정 — slab 필요)
- **좋은 seed = 강한 결합 × 높은 확산장벽 × 산화 무관 × 가시광 무흡수** → **HATCN이 유일하게 다 만족**. phen계(p-bPPhenB)는 동급이나 실증·열안정성 강점.
- **HATCN 크랙**: 결정화 수축 인장응력 → 임계두께 초과 시 균열. 표준(≤10nm)에선 미보고. seed는 3~5nm면 충분(사이트 화학은 표면 1분자층에서 포화). ≥10~30nm에서 크랙.
- **광학**: 크랙의 투과도 영향 ~0(계산 −0.1%p). 손실은 haze·면저항·국소 dewetting 쪽. n-도핑 ETL 가시광 흡수는 작음(폴라론은 주로 NIR).

---

## 1. 방법론 (재현용)

**환경**: Miniforge conda env `qc` — psi4 1.11, xtb 6.7.1, ase, rdkit, scipy, openpyxl, matplotlib, simple-dftd3.
클라우드 재구축: `conda create -n qc -c conda-forge python=3.11 psi4 xtb ase rdkit scipy openpyxl matplotlib simple-dftd3 dftd3-python requests pypdf`

**디스크립터**:
1. E_b(#1) = Ag 원자 결합에너지. 성장모드 지표(Bauer): E_b vs Ag–Ag 응집(Ag₂ 1.86 eV DFT / 1.65 exp, 벌크 2.95). E_b ≪ Ag–Ag → Volmer-Weber island.
2. 확산장벽(#2) = adatom 사이트간 hopping 장벽. pinning 지표.

**계산 프로토콜**:
- 구조: RDKit → GFN2-xTB opt (기하 생성 전용; GFN2는 Ag–Ag 과대결합 5eV라 에너지엔 안 씀).
- 에너지: **PBE-D3(BJ)/def2-SVP, UKS, counterpoise(CP) 보정**, Ag/Mo/Cs에 def2 ECP.
- 대형 유기물(TPBi, p-bPPhenB)은 DFT 거리 스캔으로 xTB 기하가 최소점임 확인.
- SCF 잘 안 붙는 계(이온성 산화물, 도핑, displaced): level_shift 1.0 + damping 15%, 안되면 PLAIN 우선/SOSCF fallback.

**한계(중요)**: 전부 **단분자/소형 클러스터 모델**. 주기 slab 아님. → 이온성 고체(LiF/MoOx/Cs₂CO₃/AlOx) 절대값은 과대·과소 여지. 순위·군구분(강 0.9~1.0 / 중 0.4 / 약 0.25)은 신뢰. 정밀값은 QE slab+NEB 필요(다음 단계).

---

## 2. 확정 결과 (전부 ALL_RESULTS.json에 통합)

### 2.1 Ag 결합에너지 E_b (eV, PBE-D3BJ/def2-SVP, CP)
| 표면(사이트) | E_b | 파일 |
|---|---|---|
| 청정 Al(금속, Al13) | 2.24 | dimer_BE / al2o3 |
| Ag–Cu / Ag–Ag / Ag–Mg (이합체) | 1.93 / 1.86 / 0.96 | dimer_BE.json |
| HATCN (nitrile N) | 1.03 | psi4_binding_refined |
| Cs₂CO₃ (클러스터, 상한) | 0.90 | bphen_cs_binding |
| TPBi (benzimidazole N) | 0.89 | psi4_binding_refined |
| p-bPPhenB (phen N,N 킬레이트) | 0.87 | psi4_binding_refined |
| F4TCNQ (nitrile) | 0.97 | newcand_binding |
| B3PyMPM (pyridyl N) | 0.63 | b3pympm_binding |
| Bphen (phen 킬레이트) | 0.49 | bphen_cs_binding |
| MoOx 산소결손(Mo3O8) | 0.44 | psi4_binding_refined |
| AlOx 화학양론(Al4O6) | 0.42 | al2o3_binding |
| AlOx Al과잉(Al10O10, 과대) | 0.91 | al2o3_binding |
| MoO3 화학양론(Mo3O9) | 0.28 | psi4_binding_refined |
| TCTA(카바졸) / LiF / TAPC(아민 π) | 0.28 / 0.25 / 0.25 | htl/psi4 |
| Liq / benzene π / HATCN면흡착 | 0.17 / 0.17 / ~0 | newcand/htl |

### 2.2 확산장벽 (eV)
| 표면 | 장벽 | 신뢰도 |
|---|---|---|
| HATCN(nitrile) | 0.29 | 확정 |
| AlOx | 0.17 | 확정 |
| TAPC / 청정Al | ~0 | 확정 |
| TPBi / B3PyMPM | 미확정(클러스터로 브릿지 SCF 실패, slab 필요) | — |

### 2.3 금속 seed 비교 결론
Al·Cu만 E_b > Ag–Ag이라 열역학적 wetting 가능하나, 실공정 산화(AlOx 0.42)·가시광 흡수로 실패. Mg는 seed 아닌 Ag:Mg 공증착용.

---

## 3. 광학 (TMM, Johnson&Christy + 사용자 측정 nk)

- **이상적 평활 Ag** (스크립트 18,19): 8nm 상대T 76.5% / 15nm 48.5%(@평균). 손실은 대부분 **반사**(흡수 ~2%).
- **Island막(MG 유효매질)** (19,20): 장파장 T↑(불연속→반사↓)이나 청색 LSPR **흡수 25%**(회수불가). 손실=흡수+산란+(조밀시)반사.
- **측정 Ag15nm/HATCN** (21): 상대T ~52%@550, **실두께 13.5nm 확정**(tooling factor 0.9). nk McPeak보다 n큼=박막 크기효과+피팅.
- **MgAg 8nm** (28): 상대T 76.7%@550, 흡수 11.8%(순수Ag 4.4%의 2.7배). Mg합금은 wetting↑ 대가로 흡수↑.
- **HATCN 크랙 광학** (29): 크랙부 ΔT ±1.6%p 이내, 면적20%여도 −0.13%p. **투과도 영향 무시가능**.
- **입력 nk 위치**: 사용자 측정 nk는 `.mat`/`.xlsx`에 있음 (아래 4절). nk_JH_total.mat에 Ag_Palik/Ag_McPeak/HATCN/Bphen/BPhen_CS(=Bphen:Cs)/BCP_Li/3TPYMB(TPYMB_3)/ZnS 등 다수.

---

## 4. 외부 데이터 파일 위치 (아카이브에 미포함 — 클라우드에서 재업로드 필요)

- `nk_JH_total.mat`: `C:\Users\Junho\Dropbox\Linkstation\Simulation\LosslessEML_single_distribution\nk_JH_total.mat` — 전 재료 nk 라이브러리(400-700nm 301pt). material 구조체.
- Ag15/MgAg 엘립소: `...\uploads\...\80fe6f6b-Ag15nm_MgAg_8nm_260703.xlsx` 및 `C:\Users\Junho\Dropbox\Linkstation\Co-work\이선정\Ag15nm, MgAg 8nm_260703.xlsx`
  - 시트 260703: n#1=HATCN30, #2/#3=HATCN+Ag15nm, #4/#5=HATCN+MgAg8nm, #6=HATCN+Mg150Å
- SEM: `C:\Users\Junho\Dropbox\[2023]...\Thin Ag film analysis\[SEM]20231030\*.tif` (Al_Ag12, Al_Ag25, HATCN_Ag12 등, 16bit)
- 참고논문: Kim AFM 2015 (ZnS/Ag/ZnS + Cs2CO3/Ag/ZnS TFOLED), Downloads 폴더.

---

## 5. 소자 설계 검토 결과 (사용자 실제 스택)

구조: `Glass/HATCN(3)/Ag(6)/HATCN(3)/TAPC(100)/TCTA(10)/TCTA:B3PyMPM:Ir(ppy)2acac(25,1:1:0.1)/B3PyMPM(10)/3TPYMB(60)/3TPYMB:Cs2CO3(10,20wt%)/Ag(100)`
- **Bottom-emission** (하부 Ag6nm 반투명 양극 DMD, 상부 Ag100nm 반사 음극). 극성·주입 경로 정상.
- 관문1: **하부 6nm Ag가 HATCN 3nm seed 위 연속막인가**(면저항으로 판정). HATCN 3nm은 피복 하한 근처 — 안되면 4nm.
- 관문2: TAPC 100nm(광학 스페이서 추정) 전압 트레이드오프 → TMM으로 최적두께 확인 권장.
- B3PyMPM→3TPYMB 전자주입: LUMO 단차 ~0.2-0.3eV, OK. 3TPYMB 60nm 저이동도=전압↑.
- 3TPYMB:Cs2CO3 n-도핑: deep LUMO(-3.3)+borane라 도핑 잘됨. 20wt%는 과함, 10wt%로 낮춰도 됨.
- 순수 3TPYMB 100nm spacer 회피(전압) → 전구간 저농도 도핑(1wt% ETL+10wt% EIL) 제안. 1wt%도 전도 충분. 흡수는 폴라론이 NIR이라 가시광 k 작음(≲0.01).

---

## 6. 다음 단계 (우선순위)

1. **QE 주기 slab + NEB** (3~4단계 방법론): TPBi/B3PyMPM 확산장벽, 이온성 고체(Cs2CO3/MoOx결손/AlOx/LiF) 절대값 재검증. 클라우드 HPC에 적합. WSL+QuantumESPRESSO 또는 VASP.
2. **소자 TMM 광학 시뮬**: 측정 nk로 bottom-emission 스택 outcoupling·spacer두께·하부Ag두께 최적화.
3. **실험 검증 DOE**: HATCN 두께시리즈(2/4/6/10nm)×Ag6nm → 면저항+AFM으로 seed window·크랙 임계두께 확정.
4. **엘립소 개선**: 도핑 ETL nk를 T&R 동시피팅+물리 오실레이터+두께고정으로 재측정(B-spline 단독 불신).
5. **ZnS seed/CPL**: Ag–S 1.60eV로 강력. 표면에너지 76(최고). Kim2015 실증. 배치별(음극아래 3nm/캡핑 22nm) 검토.

---

## 7. 파일 맵
- `scripts/01~29`: 전체 파이프라인 (번호=작업순서). 각 헤더에 목적 주석.
- `runs/*.json`: 개별 결과. `ALL_RESULTS.json`: 통합.
- `REPORT.md`: 1·2차 스크리닝 정식 보고서.
- `structures/*.xyz`: 48개 최적화 구조.
- `structures_viewer.html`: 인터랙티브 3D 뷰어(19구조, 결합에너지 색상).
- `*.png/*.csv`: 광학 그래프·데이터.
- 메모리: `~/.claude/projects/.../memory/` (oled-unityeqe-project.md, qc-toolchain-setup.md).
