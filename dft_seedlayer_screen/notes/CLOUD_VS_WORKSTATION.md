# 클라우드 PySCF 대 워크스테이션 Gaussian 16 일치 여부 (2026-11-13)

계산 수준은 둘 다 PBE0-D3(BJ)/def2-SVP다. 클라우드는 `scripts/qclocal.py`, 워크스테이션은 Gaussian 16W C.01 (`PBE1PBE/Def2SVP EmpiricalDispersion=GD3BJ`)이다.

## 1. 같은 구조, 같은 계산: 총에너지

Gaussian 출력(저장소 루트의 `*.out`, 2026-09 업로드)에서 구조를 그대로 읽어 PySCF로 다시 계산했다.
비교 스크립트는 `scripts/168_g16_vs_pyscf.py`다.

| 작업 | 원자 | G16 SCF Done (Ha) | PySCF + D3 (Ha) | 차이 |
|---|---|---|---|---|
| benzene_t0_z0 (Ag 착물, doublet) | 13 | −378.742988 | −378.743125 | −3.7 meV |
| benzene_t0_z0_ref | 13 | −378.742988 | −378.743125 | −3.7 meV |
| pyridine_t0_z0 (Ag 착물) | 12 | −394.760796 | −394.760915 | −3.2 meV |

- 총에너지 차이는 원자당 약 0.3 meV로 일정하다.
- 원인은 PySCF의 density fitting(def2-universal-jkfit)과 적분 격자 차이다.
- 이 오프셋은 E_b를 구할 때 착물과 (분자 + Ag) 사이에서 대부분 상쇄된다.
- Gaussian의 SCF Done에 D3가 포함된다는 것도 확인했다. D3를 빼고 비교하면 차이가 0.1 eV 단위로 커진다.

## 2. 같은 구조, 같은 계산: 결합 에너지 E_b

출처: `runs/common_pbe0_fragment.json`의 `_validation`. 기준값은 워크스테이션 package F의 다중 사이트 Gaussian 결과다(`notes/EB_ED_METHOD_AUDIT.md`).

| 사이트 | Gaussian (워크스테이션) | PySCF (클라우드) | 차이 | 비고 |
|---|---|---|---|---|
| HATCN aza 포켓 | 1.631 | 1.631 | 0.000 | 분자 전체(30/30 원자) |
| B3PyMPM | 0.416 | 0.416 | 0.000 | fragment 32/69 원자 |
| p-bPPhenB 킬레이트 | 0.785 | 0.784 | −0.001 | fragment 32/72 원자 |
| Bphen 킬레이트 | 0.621 | 0.609 | −0.013 | fragment 22/42 원자. 차이는 fragment 절단 때문이다 |

**결론: 같은 구조와 같은 사이트에서는 두 코드가 1 meV 안에서 일치한다. 분자를 잘라 계산할 때만 최대 13 meV 차이가 난다.**

## 3. 숫자가 다른데 코드 문제가 아닌 경우

| 값 | 어디서 | 왜 다른가 |
|---|---|---|
| HATCN 1.631 | Gaussian package F, PBE0로 Ag 위치 최적화, aza 포켓 | 기준값 |
| HATCN 0.604 | Gaussian stage 1 (`runs/binding_energies_pbe0.json`) | 니트릴 포켓 사이트. 사이트가 다르다 |
| HATCN 1.03 | 초기 REPORT와 HANDOFF, PBE-D3BJ, CP 보정 | 범함수(PBE)와 사이트가 다르다 |
| HATCN 1.211 | 클라우드 설계 앵커 (`runs/design/pbe0_anchor.json`) | 같은 PBE0지만 분자는 GFN2 구조이고, Ag를 GFN2 위치에 고정했다. 차이 −0.42 eV는 모두 전자 항이다 |
| HATCN 1.770 | 클라우드, 환원된 HATCN⁻Cs⁺ | 다른 계(아래층에서 전자를 받은 HATCN) |

## 4. 비교할 수 없는 것

- **Bphen:Cs₂CO₃**: 워크스테이션에서는 기준 분자(BphenCs_mol)만 끝났다. 사이트 작업은 멈췄다. 그래서 클라우드 fragment 결과(Cs 사이트 0.886 eV)와 짝지을 Gaussian 값이 없다.
- **설계 앵커 5개, TD-DFT**: 클라우드에서만 계산했다.
