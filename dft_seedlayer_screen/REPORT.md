# Seed layer 스크리닝: 6 nm 열증착 Ag의 wetting 예측 (DFT)

**날짜**: 2026-07-04
**후보**: p-bPPhenB, HATCN, MoOx, LiF, TPBi
**목적**: 6 nm Ag 열증착 시 Volmer–Weber island 성장을 억제할 seed layer의 사전 선별

---

## 1. 방법

로컬 워크스테이션에 계산 환경을 구축했다 (Miniforge → conda env `qc`: Psi4 1.11, xTB 6.7.1, ASE, RDKit).

**디스크립터**: seed 표면 위 Ag 원자 결합에너지 E_b.
금속 성장모드의 1차 지표로, E_b가 Ag–Ag 응집(벌크 2.95 eV/atom, Ag₂ 이량체 1.65 eV 실험값)에 가까울수록 핵생성 밀도가 높고 admolecule 확산이 억제되어 조기 percolation(연속막)에 유리하다. E_b ≪ Ag–Ag이면 흡착원자가 표면을 자유롭게 확산·응집하여 island(Volmer–Weber) 성장한다.

**모델** (분자/클러스터 근사):
| Seed | 모델 | Ag 흡착 사이트 |
|---|---|---|
| HATCN | 단분자 (C₁₈N₁₂) | ① 방향족 코어 면흡착 ② 말단 nitrile N |
| TPBi | 단분자 (C₄₅H₃₀N₆) | benzimidazole 피리딘형 N |
| p-bPPhenB | 단분자 (C₄₂H₂₆N₄) | phenanthroline N,N 킬레이트 포켓 |
| MoOx (x=3) | Mo₃O₉ 클러스터 (표준 기상 MoO₃ 모델) | 링 상부 |
| MoOx (x<3) | Mo₃O₈ (말단 O 제거 = 산소결손) | 결손 사이트 |
| LiF | (LiF)₃₂ 4×4×2 암염 클러스터 (벌크 격자 고정, a=4.026 Å) | 표면 F-top |

**절차**: RDKit 3D 생성 → GFN2-xTB 기하 최적화 → **PBE-D3(BJ)/def2-SVP** (Ag/Mo에 def2-ECP, UKS, counterpoise 보정) 단일점. GFN2 기하가 PBE 기준으로 계통적으로 과밀착이라 Ag–기판 거리를 DFT 강체 스캔으로 재최적화(소형 계). 검증: 동일 수준에서 Ag₂ 결합에너지 1.86 eV (실험 1.65 eV, 오차 +13% — 수용 가능).

주의: GFN2-xTB 자체는 Ag 결합을 심하게 과대평가(Ag₂ 5.0 eV)하므로 **기하 생성에만** 사용하고 에너지는 전부 DFT 값이다.

---

## 2. 결과

PBE-D3(BJ)/def2-SVP, CP 보정된 Ag 원자 결합에너지:

| 순위 | Seed (사이트) | E_b (eV) | 비고 |
|---|---|---|---|
| 1 | **HATCN** (nitrile N) | **1.03** | 강한 Ag–N≡C 결합; 분자당 CN 6개 → 고밀도 트랩 사이트 |
| 2 | **TPBi** (benzimidazole N) | **0.89** | 분자당 접근 가능 N 3개 |
| 3 | **p-bPPhenB** (phen N,N 킬레이트) | **0.87** | 이좌배위; 분자당 포켓 2개 |
| 4 | MoOx 환원면 (Mo₃O₈ 결손) | 0.44 | 화학양론면의 ~2배 |
| 5 | MoO₃ 화학양론면 (Mo₃O₉) | 0.28 | 약함 — 무결함 terrace는 seed 역할 못함 |
| 6 | **LiF**(001) F-top | **0.25** | 순수 물리흡착 |
| – | HATCN 코어 면흡착 | ≈0 | 비작동 사이트 (CN 사이트가 지배) |
| 기준 | Ag–Ag (Ag₂ 이량체) | 1.86 | 실험 1.65 eV |
| 기준 | Ag 벌크 응집에너지 | 2.95 (실험) | |

(TPBi·p-bPPhenB는 이후 DFT 거리 스캔으로 xTB 기하가 이미 최소점임을 확인 — 두 값 모두 수렴값이며 차이 0.02 eV는 방법 오차 이내로 **동급**.)

모든 seed에서 E_b < E_coh(Ag) — 열역학적으로는 어떤 표면에서도 Ag는 Ag 위를 선호한다(3D 성장 경향). 따라서 wetting 제어는 **속도론**(핵생성 밀도·확산 억제)이 관건이고, E_b가 큰 표면일수록 유리하다.

### 예측
- **HATCN, p-bPPhenB, TPBi**: E_b ≈ 0.9–1.0 eV (Ag₂의 절반 이상) → 실온 확산이 크게 억제되어 핵생성 밀도가 높고, 6 nm에서 연속막 형성 가능성이 가장 높은 군.
- **MoOx**: 화학양론 MoO₃ 표면 자체는 약함(0.28 eV). 증착 MoOx의 실험적 seeding 효과는 **산소결손·비화학양론**(0.44 eV, 실제 결손 사이트는 더 강할 것)과 높은 표면에너지에서 온다 → 환원된 MoOx(x≈2.7–2.9) 또는 UV-오존 처리로 극대화.
- **LiF**: 0.25 eV — 명백한 anti-wetting. 단독 seed로는 부적합(island 성장 확정적).

---

## 3. 문헌 대비 검증

| 계 | 본 계산 | 문헌 |
|---|---|---|
| Ag/LiF | 0.25 eV, island 예측 | Ag는 island 형성이 알려져 있고 LiF 위 직접 증착은 dewetting; LiF는 Al seed와 병용해야 함 ([Al-doped Ag cathode, Org. Electron. 2019](https://www.sciencedirect.com/science/article/abs/pii/S1566119919304379)) — **일치** |
| Ag/MoO₃ | terrace 0.28 / 결손 0.44 eV | 1 nm MoO₃ seed로 7 nm percolation, 공기 어닐링·UVO로 표면에너지↑ 시 wetting 개선 ([Nanomaterials 2018](https://www.mdpi.com/2079-4991/8/7/473), [Sol. RRL 2018](https://www.sciencedirect.com/science/article/abs/pii/S0927024818303192)) — 결함 기구로 **정합** |
| Ag/HATCN | nitrile N 1.03 eV | HATCN/Ag(111) 흡착 -2.97~-3.05 eV/molecule (여러 CN 동시 결합), CN기가 Ag와 특이적 상호작용 ([arXiv 2312.08233](https://arxiv.org/pdf/2312.08233)) — 방향 **일치** (분자당 다중 CN ⇒ 흡착합산) |
| Ag/phenanthroline | 킬레이트 0.87 eV | BCP:Ag에서 Ag–phen 배위결합 형성, Ag 확산 억제 보고 ([Org. Electron. 2014](https://www.sciencedirect.com/science/article/abs/pii/S1566119914001608)); **p-bPPhenB/Ag/p-bPPhenB 샌드위치 전극에서 p-bPPhenB가 wetting inducer로 실제 사용됨** — **일치** |
| Ag/TPBi | 0.89 eV | 직접적 wetting 보고는 없음; TPBi/Ag 계면 UPS 연구 존재. 유사 N-헤테로고리 ETL 위 Ag는 중간 수준 상호작용 |
| Ag₂ 기준 | 1.86 eV | 실험 1.65 eV — 방법 검증 (E_b 축) |
| TCNQ EA 기준 | 3.54 eV (wB97X/def2-TZVP) | 실험 3.383 eV (TCNQ⁻ 저온 광전자분광, X.-B. Wang 그룹) — +4.5%, **방법 검증 (EA 축)**. 상대 검증: EA(F4TCNQ)−EA(TCNQ) 0.48 vs 문헌 ~0.55 eV |

기타 문헌 보고 seed: Ge ([Nano Lett. 2009](https://pubs.acs.org/doi/abs/10.1021/nl8027476)), Al ([Micromachines 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9565528/)), PEDOT:PSS, PEI, ZnO, Cu 등. 금속 seed(Ge, Al, Cu)는 강력하지만 가시광 흡수 페널티가 있어 top-emitting 소자에서는 유기/산화물 seed가 유리.

---

## 4. 결론 및 권고

6 nm Ag top electrode용 seed layer 우선순위:

1. **HATCN** — 최강 결합(1.03 eV) + 고밀도 CN 사이트. ETL 위 전자주입층으로 이미 얇게(≤10 nm) 쓰이므로 광학 페널티 최소. 최우선 실험 후보.
2. **p-bPPhenB** — 킬레이트 결합(≥0.87 eV) + 문헌에서 wetting inducer로 실증됨. ETL 겸용 가능(고이동도)이라 소자 구조 단순화.
3. **TPBi** — p-bPPhenB와 동급 결합이나 실증 문헌 부재. ETL로 이미 쓰는 구조라면 추가층 없이 시도 가능.
4. **MoOx** — 단, **비화학양론(환원) 상태** 확보가 조건 (증착 후 처리·저산소 증착). 화학양론 MoO₃면 기대 이하일 것.
5. **LiF** — 단독 사용 배제. 사용하려면 LiF/Al(≤1 nm) 조합으로.

**한계**: 분자/클러스터 모델(주기적 박막 표면·분자 배향·패킹 미반영), Ag 원자 1개 디스크립터(확산장벽·클러스터 단계 미포함), PBE-D3 수준, 대형 유기물 2종은 하한치. 정량 랭킹은 문헌 정합성으로 보강했으며, 절대값보다 **군 구분**(강: 유기 N-사이트 ≈ 0.9–1 eV / 중: 환원 MoOx / 약: LiF·화학양론 MoO₃)이 신뢰 구간이다.

**후속 계산 제안**: (i) Quantum ESPRESSO 주기 slab으로 MoO₃(010) 결손 표면 + Ag, LiF(001) + Ag 검증, (ii) Ag 확산장벽(NEB), (iii) 유기막 표면의 분자 배향별 N-사이트 노출도(MD).

---

## 5. 추가 분석 (2차): 금속 seed 비교, p-bPPhenB 기구, 투명 후보군

**금속 seed와의 비교** (Ag–M 이합체, 동일 수준): Ag–Al 2.24 / Ag–Cu 1.93 / Ag–Ag 1.86 / Ag–Mg 0.96 eV. Al·Cu만 Ag–Ag보다 강해 **열역학적 wetting**이 가능한 유일한 부류(문헌의 Cu·Al seed 실증과 일치). Mg는 seed가 아닌 공증착 합금(속도론적 억제)으로 작동. 단 금속 seed는 가시광 흡수 페널티로 top-emitting 목적에는 불리.

**p-bPPhenB vs TPBi (ACS Photonics 2019, 10.1021/acsphotonics.9b01155 관련)**: 정적 E_b는 0.87 vs 0.89 eV로 동급임을 DFT 거리 스캔까지 확정. 논문의 p-bPPhenB 우수성은 단일원자 흡착에너지가 아니라 (i) 이좌배위 킬레이트의 높은 탈출장벽(두 결합 동시 파단 필요 → 유효 체류시간 수 자릿수 증가), (ii) 성장 후 Ag 표면 원자 킬레이션에 의한 재응집/dewetting 억제(110 °C 열안정성의 기원), (iii) p-bPPhenB/Ag/p-bPPhenB 샌드위치의 상부 캡핑 효과에서 기인하는 것으로 해석됨. 정량 검증에는 탈출장벽(NEB)과 Ag₂/Ag₃ 클러스터-킬레이트 결합 계산이 필요.

**가시광 무흡수 제약 반영 추가 후보**: C60(가시광 흡수), F4TCNQ/F6-TCNNQ(CT 밴드; F4TCNQ nitrile 사이트 계산값 0.97 eV였으나 광학 탈락), NTCDA/PTCDA 제외. 잔존 후보: Phen-NaDPO(phen 킬레이트+P=O), B4PyMPM(수평배향 4-피리딜 N), TmPyPB, Bphen/BCP, PO-T2T — 모두 E_g ≥ 3.3 eV. Liq는 계산 결과 0.17 eV로 탈락(음성 대조). 투명성 제약 하 최상위군은 여전히 HATCN + phen계. MoOx는 wetting을 돕는 산소결손이 가시광 갭 스테이트 흡수를 유발하는 트레이드오프에 주의.

---

## 파일 구성
- `structures/` — 초기 구조 (RDKit/수작업)
- `runs/<계>/xtbopt.xyz` — GFN2 최적화 기하; `runs/psi4_binding_eV.json` (단일점), `runs/psi4_binding_refined_eV.json` (스캔 보정 최종값)
- `scripts/01~08` — 재현용 전체 파이프라인
