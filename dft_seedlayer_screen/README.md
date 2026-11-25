# Ag seed layer screening — DFT + optics

6~15 nm 열증착 Ag를 island 없이 균일하게 올릴 seed layer를 DFT로 사전 선별하고,
그 전극의 광학(투과/반사/흡수)을 정량화하는 프로젝트.

인수인계 문서는 **`HANDOFF.md`** 를 먼저 읽을 것. 이 README는 저장소 구조만 설명한다.

## 문서

| 파일 | 내용 |
|---|---|
| `HANDOFF.md` | 프로젝트 전체 인수인계 (방법론·확정결과·다음단계) |
| `REPORT.md` | 1차 스크리닝 정식 보고서 (p-bPPhenB / HATCN / MoOx / LiF / TPBi) |
| `REPORT2_beyond_hatcn.md` | 2차: "HATCN을 넘는 seed가 있는가" 전수 스크리닝 결론 |
| `LIT_NOTES_beyond_hatcn.md` | 문헌 신규성 판정 노트 (Park & Suh 2018 선행 반영한 정정 포함) |
| `PRIOR_STACKS_HATCN_Ag.md` | HATCN–Ag 계면 선행 스택 목록 + 리뷰어 방어 논리 |

## 스크립트 (`scripts/`, 번호 = 작업 순서)

- `01`–`17` : 1차 DFT 스크리닝 파이프라인 (RDKit → GFN2-xTB → PBE-D3BJ/def2-SVP CP)
- `18`–`21`, `28`–`29` : TMM 광학 (Ag 투과/반사/흡수, 측정 nk, MgAg, 크랙)
- `22`–`27` : AlOx, 확산장벽
- `30`–`31` : 2차 스크리닝 (CuI/CuSCN/포스핀옥사이드/BTD/트리아진/티오펜, Cu 금속사이트)
- `32`–`33` : 금속 seed(Al/Au) 흡수 페널티, 투과도 검증 (독립 TMM 2종 교차검증)
- `34` : HATCN 30nm 거칠기·크랙 시나리오 광학
- `35` : SEM 텍스처 정량분석 (로컬 실행용)

## 결과

- `runs/*.json` — 개별 계산 결과 (스캔 전곡선 포함)
- `ALL_RESULTS.json` — 1차 통합
- `runs/xtb/*/xtbopt.xyz` — xTB 최적화 최종 구조
- `structures/*.xyz` — 최적화 구조 (48+)
- `optics/` — TMM 그래프·CSV
- `structures_viewer.html` — 인터랙티브 3D 뷰어

psi4 원시 로그(`psi4_*.out`, 총 7MB)와 xTB 스크래치는 용량 문제로 미포함.
스크립트를 재실행하면 동일하게 재생성된다.

## 환경 재구축

```bash
conda create -n qc -c conda-forge python=3.11 psi4 xtb ase rdkit scipy \
    openpyxl matplotlib simple-dftd3 dftd3-python requests pypdf
```

검증됨: psi4 1.11, xtb 6.7.1, ase 3.29.0, rdkit 2026.03.4.

## 핵심 결과 요약

Ag 결합에너지 E_b (eV, PBE-D3BJ/def2-SVP, counterpoise):

```
청정Al 2.24 > Ag-Au 2.14 > Ag-Cu 1.93 > Ag-Ag 1.86 > ZnS 1.60
  > HATCN 1.03 > F4TCNQ 0.97 > Cs2CO3 0.90 > TPBi 0.89 > p-bPPhenB 0.87
  > CuI(Cu site) 0.82 > B3PyMPM 0.63 > CuSCN(Cu) 0.57 > Bphen 0.49
  > MoOx결손 0.44 > AlOx 0.42 > s-triazine 0.29 > MoO3 0.28 = BTD 0.28
  > TCTA 0.28 > Me3PO 0.25 = LiF 0.25 = TAPC 0.25 > CuI(I) 0.22
  > Liq 0.17 = thiophene 0.17 > CuSCN(S/N) 0.14/0.10
```

- 증착 가능·투명·비산화 재료 중 **HATCN(1.03 eV)을 넘는 것은 없다**.
  더 강한 두 채널(금속-금속, 설파이드)은 각각 산화 취약(Al→AlOx 0.42)과
  선행보고(ZnS, Kim AFM 2015)로 막혀 있다.
- 금속 seed의 광학 대가: 흡수 구동인자 Im(eps)=2nk 가 Ag 0.43 / Au 2.56 / Al 12.8 (550nm).
  캡핑 포함 상대투과도(Ag 9nm 기준): HATCN 99.0% / Au 2nm 94.4% / Al 1nm 86.4%.
