# 로컬 세션 인수인계 (LOCAL HANDOFF)

> 이 문서는 **클라우드 세션 → 로컬 Claude Code 세션** 이전용입니다.
> 저장소 루트를 지정하면 이 문서만 읽고 바로 이어서 작업할 수 있도록 작성했습니다.
> 원본 프로젝트 설명은 `HANDOFF.md`, 결과 해석은 `REPORT*.md`에 있습니다.

**마지막 클라우드 세션 종료 시점 커밋:** `dae5433`
**브랜치:** `claude/project-setup-conda-ftxhnl`

---

## 0. 새 세션에서 가장 먼저 할 일

```bash
git clone <repo> && cd OLED-
bash dft_seedlayer_screen/scripts/bootstrap.sh      # 환경 검증
```

그다음 **§5 미결 과제**의 1번(TCPM 슬랩 Ag 기준값 진단)부터 시작하면 됩니다.
**§4의 "폐기된 주장"을 반드시 먼저 읽으십시오.** 이미 틀린 것으로 밝혀진 결론들이 있고,
그걸 모르면 같은 실수를 반복하게 됩니다.

---

## 1. 왜 로컬로 옮기는가

클라우드 컨테이너가 **스냅샷 복원 방식으로 반복 초기화**되었습니다. 증상: 로컬 `.git`이
과거 커밋으로 되돌아가고, 실행 중 프로세스가 사라지고, `/tmp`가 비워지고, pip 설치분이
소멸. `setsid`/`nohup`/`disown` 모두 무효 — 프로세스를 죽이는 게 아니라 기계를 교체하는
것이기 때문입니다. 세션이 유휴 상태일 때 주로 발생했고, 최소 8회 겪었습니다.

원격 git에 푸시된 것만 살아남았으므로 **체크포인트 구조**(`scripts/_ckpt.py`)를 만들어
계산 단위마다 즉시 커밋·푸시하도록 했습니다. 완료된 계산은 이후 하나도 잃지 않았지만,
**진행 중이던 SCF 하나는 매번 잃었고** 그게 20분짜리라 진도가 나가지 않았습니다.

로컬 환경(14코어 / 32 GB)이 유리한 순서:

| 순위 | 이유 |
|---|---|
| 1 | **롤백 없음** — 이게 압도적입니다 |
| 2 | **RAM 32 GB** — 격자 예산 9 → 20 GB로 올려 셀·CUTOFF 타협 해소 |
| 3 | CPU 2~3배 |
| 4 | 학교 IP → **문헌 열람 가능** (클라우드는 전 호스트 403) |

GPU 없는 것은 무관합니다. CP2K/psi4의 이 워크로드는 GPU를 거의 쓰지 않습니다.

---

## 2. 환경 구축 (Windows → WSL2)

CP2K는 Windows 네이티브 빌드가 사실상 없으므로 **WSL2**를 씁니다. 그러면 클라우드와
완전히 동일한 환경이 되어 스크립트를 한 줄도 고칠 필요가 없습니다.

### 2.1 WSL2

```powershell
wsl --install -d Ubuntu
```

**메모리 상향이 중요합니다.** WSL2는 기본적으로 호스트 RAM의 절반(16 GB)만 씁니다.
`C:\Users\<사용자>\.wslconfig` 생성:

```ini
[wsl2]
memory=24GB
processors=12
swap=8GB
```

적용: `wsl --shutdown` 후 재시작.

### 2.2 conda

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
source ~/miniforge3/etc/profile.d/conda.sh
```

### 2.3 두 환경

`env/qc.yml`, `env/cp2k.yml`에 클라우드 환경을 그대로 내보내 두었습니다.

```bash
conda env create -f dft_seedlayer_screen/env/qc.yml      # psi4 1.11, xtb 6.7.1, ase 3.29, rdkit
conda env create -f dft_seedlayer_screen/env/cp2k.yml    # cp2k 2026.2
```

정확한 재현이 필요 없으면 더 간단히:

```bash
conda create -n qc   -c conda-forge psi4 xtb ase rdkit numpy scipy matplotlib
conda create -n cp2k -c conda-forge cp2k
```

### 2.4 경로 수정이 필요한 곳

스크립트에 클라우드 절대경로가 박혀 있습니다. **`scripts/_ckpt.py` 상단**:

```python
REPO = "/home/user/OLED-"                      # → 로컬 clone 경로로
BRANCH = "claude/project-setup-conda-ftxhnl"   # → 그대로 두거나 로컬 브랜치명으로
```

CP2K 실행 파일 경로가 여러 스크립트에 하드코딩돼 있습니다:

```python
cp2k = "/root/miniforge3/envs/cp2k/bin/cp2k.psmp"
CP2K_DATA_DIR="/root/miniforge3/envs/cp2k/share/cp2k/data"
```

`bootstrap.sh`가 이 경로들을 점검하고 알려줍니다.

### 2.5 절전 해제

슬랩 계산은 몇 시간 단위입니다. 사무실 PC의 절전/최대절전/자동종료를 꺼두십시오.

---

## 3. 저장소 구조

```
dft_seedlayer_screen/
├── scripts/          61개. 번호 순서가 곧 작업 순서
│   ├── _ckpt.py      체크포인트 (계산 단위마다 커밋·푸시)
│   ├── bootstrap.sh  환경 검증
│   ├── 30~31         2차 클러스터 스크리닝
│   ├── 32~34         광학 (TMM, 흡수, 거칠기/크랙)
│   ├── 36            kMC 성장 (Venables 검증)
│   ├── 37            면저항 (FS + MS)
│   ├── 38~41         슬랩 인프라 + HATCN 단분자층 + pocket 분리
│   ├── 42            문헌 수집기 ★ 로컬에서 돌리도록 작성됨
│   ├── 43            F4TCNQ 슬랩
│   ├── 44            TD-DFT 흡수
│   ├── 45            Ag 산화(염 형성) 열역학
│   ├── 46~48         HATCN k 민감도 / 이론 k / 계면 CT 흡수
│   ├── 49,52         F6TCNNQ
│   ├── 50,51,53      크랙 없는 seed 스크리닝 / 킬레이트 / 3D 억셉터
│   ├── 54            후보 EA 순위
│   └── 55            TCPM 슬랩 ★ 진행 중
├── runs/             결과 JSON (이게 진짜 산출물), xtb 최적화 좌표
├── slabs/            CP2K 입력(.inp)·구조(.xyz). .out/.wfn은 gitignore
├── structures/       모든 분자 좌표
├── env/              conda 환경 사양
├── MANUSCRIPT_DRAFT.md   논문 초안 (Figure별 데이터 계획)
├── REPORT*.md            결과 해석
└── LOCAL_HANDOFF.md      이 문서
```

**결과는 `runs/*.json`에 있습니다.** 콘솔 로그(`runs/*.log`)와 CP2K `.out`은
gitignore 대상이고 재생성 가능합니다.

---

## 4. 결과 현황 — 신뢰도별

### 4.1 확립된 것

| 항목 | 값 | 근거 |
|---|---|---|
| 클러스터 E_b 최상위 | HATCN **1.029** eV | PBE-D3BJ/def2-SVP + counterpoise |
| 〃 | F4TCNQ 0.966, CuI(Cu) 0.82, CuSCN(Cu) 0.57, TAPC 0.02 | 동일 프로토콜 |
| 슬랩 단분자층 E_b | HATCN **1.346**, F4TCNQ **1.556** eV | scripts/39, 43 |
| **pocket 기여** | **+0.534 eV** (HATCN, a=15 vs 19) | scripts/41 |
| 격자 오프셋 상쇄 | E_b 편차 0.014 eV vs E_total 0.376 eV | scripts/40 |
| kMC Venables | **χ = 0.284 ± 0.013** (i=1) | scripts/36 |
| TMM 검증 | 독립 구현 2개가 1.9×10⁻¹⁵ 일치 | scripts/33 |
| 금속 seed 흡수 | Al 0.5 nm 단독 **5.43 %** (8 nm Ag 전체보다 큼) | scripts/32 |
| EA (def2-TZVP) | HATCN 3.378, TCNQ 3.536, F4TCNQ 4.020 | scripts/45, 시차 게이트 통과 |
| 염 형성 격차 | HATCN이 F4TCNQ보다 **1.28 eV 불리** | scripts/45, **상대값만** |
| 중성 Ag 추출 | 어느 분자도 bulk/kink/Ag₂에서 못 뜯음 | scripts/45 |
| 조각 E_b | PhCN 0.292, pDCNB 0.379, DMABN 0.287 | scripts/50 |
| **킬레이트 실패** | o-dinitrile **0.130** < 단일 CN 0.292 | scripts/51 |
| **EA↔E_b 대리지표** | **r = 0.962**, 기울기 0.055 eV/eV | scripts/53 |
| 결정화 위험 | HATCN 92, F4TCNQ 88, F6TCNNQ 79 / TCPM·TCPSi 35, SBF2CN 32, mCPCN 23, TPBi 17 | scripts/53 |
| EA (def2-SVP) | HATCN 2.711, TCPM 0.872, SBF2CN 0.540 | scripts/54 |

### 4.2 잠정 — 검증 필요

| 항목 | 문제 |
|---|---|
| F4TCNQ 기체상 EA ≈ 3.9 eV | **기억에서 나온 값.** NIST WebBook 확인 필요. 널리 인용되는 5.24 eV는 **고체상**(UPS/IPES)이며 분극 ~1.2 eV 포함 |
| F6TCNNQ EA 4.598 eV | def2-SVP 차이를 def2-TZVP 기준에 옮긴 **추정치** |
| HATCN k(가시광) ≤ 2×10⁻⁴ | 계산값. **흡수단 위치가 지배 변수** — 0.25 eV 차이로 여유가 145배↔3배로 뒤바뀜 |
| Tg 값들 | 기억 기반. 순위 검증용일 뿐 데이터 아님 |
| 4CzIPN/2CzPN 가시광 흡수 | 화학 지식. 확인 필요 |

### 4.3 폐기된 주장 — **반복하지 말 것**

| 폐기된 주장 | 실제 |
|---|---|
| ~~"HATCN/Ag는 선행 보고 없음"~~ | **Park & Suh, Opt. Express 26, 4979 (2018)** — HATCN(7)/Ag(15–25 nm) |
| ~~"HATCN 흡수 기여 0.00 %p"~~ | k=0 **가정**을 되돌려받은 값. 순환논법 |
| ~~"F4TCNQ는 흡수 때문에 탈락"~~ | TD-DFT는 반대: HATCN 음이온 0.299 > F4TCNQ 0.197 |
| ~~"HATCN이 가장 강한 seed"~~ | 슬랩에서 F4TCNQ 1.556 > HATCN 1.346 |
| ~~"30 nm 막이면 k=10⁻³ 검출 가능"~~ | 실제 흡수 0.07~0.09 %, UV-Vis 재현성 이하 |
| ~~F6TCNNQ 수렴 실패 = Ag 산화 증거~~ | 20 Å에서도 실패. **수렴 안 한 계산은 아무것도 측정하지 않음** |

수정된 계산 버그: Kapustinskii 단위 10배 오류, 기체/고체 EA 혼용,
TMM 부호 규약, 반사식을 투과 산란에 오용, Maxwell-Garnett 특이점,
Mayadas-Shatzkes D=d 가정, TCPM 슬랩 앵커 선택 — 모두 수정 완료.

---

## 5. 미결 과제 (우선순위)

### 1위 — TCPM 슬랩 Ag 기준값 진단 ★ 진행 중

`scripts/55_tcpm_slab.py`. 앵커 버그를 고친 뒤 **r = 4.0 Å에서 E_b = +2.72 eV**가
나옵니다. **물리적으로 불가능합니다** — HATCN의 최대값이 2.2 Å에서 1.346 eV인데
약한 분산력만 있어야 할 4.0 Å에서 두 배가 나올 수 없습니다.

현재 체크포인트 (`runs/tcpm_slab.json`):
```
ml       -214.17348105 Ha   (깨끗한 단분자층)
ag        -36.93630815 Ha   (고립 Ag, 같은 셀)   ← 의심 대상
cx_r4p0  -251.20988866 Ha
```

**가장 의심되는 것: 고립 Ag 기준값.** E_b = E_ml + E_ag − E_cx에서 E_ag가 ~2.5 eV
높으면 정확히 이 증상입니다. 고립 원자에 smearing(1000 K)을 적용하는 구성이라
5s 점유가 잘못 수렴했을 가능성이 있습니다.

`slabs/tcpm/ag.inp`의 설정이 고립 원자에는 과합니다:

```
ADDED_MOS 60          # 원자 1개에 빈 오비탈 60개
SMEAR 1000 K          # 페르미 분포
CELL 13.16 x 13.16 x 26.32 A @ 350 Ry
```

빈 상태 60개에 분수 점유가 퍼지면 에너지가 실제 바닥상태와 달라질 수 있습니다.

**진단 순서:**
1. 다른 셀 크기 2~3개에서 고립 Ag 에너지 재계산 → 안정한지
2. smearing **과 ADDED_MOS를 함께** 제거하고 고정 점유로 비교
   (클라우드에서 SMEAR만 지우고 ADDED_MOS를 남겨 무효한 비교를 한 적 있음 — 둘 다 지울 것)
3. `slabs/tcpm/ag.out`의 점유수·수렴 이력 확인
4. scripts/41(HATCN, 정상 작동)의 Ag 기준값과 대조 — 같은 구성인데 왜 거기선 됐는지
5. 결정적 검증: 같은 셀에서 Ag-N = 10 A 짜리 복합체를 계산.
   그 거리면 **E_b가 0에 가까워야** 하고, 아니면 기준값이 틀린 것이 확정됨

**클라우드에서 이 진단을 시도했다가 실패했습니다.** 4코어로 원자 1개 계산이 240초 안에
안 끝났습니다 — 큰 셀 + 350 Ry + ADDED_MOS 60의 대각화가 무겁습니다.
12코어에 시간 제한이 없으면 문제없이 끝납니다.

### 2위 — TCPM interstitial 자리

현재 스캔은 **on-top(위로 선 nitrile)** 자리입니다. 그 자리에서 adatom–이웃 분자
거리가 **12.5 Å**(HATCN 3.06 Å)이라 pocket이 원리적으로 불가능합니다.

**그러나 TCPM에서 HATCN pocket에 해당하는 건 분자 사이 골(interstitial)입니다.**
격자 13.16 Å에 분자 폭 ~11 Å이라 공간이 있고, 거기라면 adatom이 여러 이웃의 팔에
동시에 닿을 수 있습니다. **이 자리는 시험하지 않았습니다.** pocket이 존재할 수 있는
유일한 경로이므로, 이걸 안 하고 "TCPM에 pocket 없음"이라 결론지으면 안 됩니다.

### 3위 — F6TCNNQ E_b

Ag가 있으면 SCF가 **4.0 eV 떨어진 두 상태(중성 ↔ Ag⁺+음이온)를 무한 왕복**합니다.
20 Å에서도 동일 → 물리적 상호작용이 아니라 범함수/알고리즘 문제.
시도했다 실패한 것: maxiter 500, 감쇠 30/40/60 %, level shift 3종, 초기추측 4종,
far→near 궤도 이어받기, SOSCF.

**해결책: 하이브리드 범함수.** PBE의 비편재화 오차가 원인이라면 PBE0/wB97X가 두
전하 상태를 제대로 분리합니다. 단 **비교 가능성 때문에 HATCN·F4TCNQ·F6TCNNQ 셋을
같은 수준으로 재계산**해야 합니다. 32 GB / 14코어면 현실적입니다.

### 4위 — 문헌 확인 (로컬에서만 가능)

```bash
conda activate qc && python dft_seedlayer_screen/scripts/42_lit_harvest.py
```

OpenAlex 키워드 + 인용 확장 + Unpaywall 합법 OA만 내려받습니다. 출판사 스크래핑은
의도적으로 배제했습니다. 확인해야 할 것:
- F4TCNQ **기체상** EA (§4.2)
- Ag-TCNQF4 형성 조건이 상온 증착과 겹치는지
- HATCN 보고된 n,k / 흡수단 위치

### 5위 — 나머지

- TCPSi EA (미실행, TCPM과 유사할 것으로 예상)
- 슬랩 E_b에 BSSE 보정 없음 → 현재 값은 상한
- kMC `E_d ≈ 0.28·E_b` 가정: 실제 NEB 비율이 45배 분산. **SI에 민감도 분석 필수**
- Ag(111) 슬랩 + Bader 전하 → 염 형성 절대 판정 (현재는 상대값만)

---

## 6. 반복해서 당한 함정들

| 함정 | 대응 |
|---|---|
| `pkill -f <패턴>`이 **자기 명령줄을 매칭**해 셸까지 종료 | PID로 종료 (`ps -C python -o pid=`) |
| psi4 `e_convergence 1e-8`이 **DF 잡음 바닥 아래** → 영원히 수렴 안 함 | 잡음을 실측해서 설정 (26원자에서 7×10⁻⁶ Ha → 1×10⁻⁵) |
| psi4 스크래치 **22 GB** 누적 → 디스크 꽉 참 → PSIO 쓰기 실패 | 계산 사이에 `/tmp/psi.*`, `/tmp/psi4_*.npy` 삭제 |
| CP2K 자유에너지 vs T→0 외삽 혼용 → E_b 0.10 eV 과대 | `"extrapolated to T->0"` 우선 파싱 |
| `argmax(z × [1 if N else -1e9])` — z가 음수면 비질소가 이김 | 명시적 인덱스 전달 + assert |
| xtb가 Ag를 과결합 (킬레이트 2.15 Å vs DFT 3.2 Å) | 구조는 xtb, **에너지는 반드시 DFT** |
| `_ckpt.commit()`이 스테이징된 **다른 파일까지** 커밋 | `git commit --only <path>`로 범위 제한 (미수정) |
| 폴링 루프에서 cwd 리셋 | 절대경로 사용 |

---

## 7. 체크포인트 사용법

```python
from _ckpt import Checkpoint
ck = Checkpoint("dft_seedlayer_screen/runs/myresults.json")

for job in jobs:
    if ck.has(job):
        continue              # 이전 실행에서 완료됨
    ck.put(job, expensive(job))   # 원자적 쓰기 → 커밋 → 푸시
```

**로컬에서는 push를 꺼도 됩니다** (`Checkpoint(..., push=False)`) — 롤백이 없으니까요.
다만 커밋은 유지하는 편이 이력 추적에 좋습니다.

**단위를 작게 잡으십시오.** 6점 스캔을 복합체 단위로 체크포인트하면 스캔 전체를 잃고,
스캔점 단위로 하면 SCF 하나만 잃습니다.

---

## 8. 논문 상태

`MANUSCRIPT_DRAFT.md`에 Figure 1~5 패널별 데이터 계획이 있습니다. 핵심:

- 프레임이 **"HATCN이 최고"가 아니라 "결합력 vs 부식 상충 설계규칙"** 입니다.
  선행논문(Park & Suh 2018) 때문에 물질 최초 제안은 불가능하고, 기전 규명이 기여입니다.
- **Fig 1(b)의 red box는 근거가 없어 철회**했습니다. 산화 임계 EA 값은 알려진 게 없고,
  우리 사이클은 실존하는 Ag-TCNQF4조차 "안 생김"으로 계산합니다. **격차만 인용**하십시오.
- Fig 1(b)를 채우려면 유기물 전반의 EA 계산이 더 필요합니다 (현재 완전한 점이 적음).
- **XPS Ag 3d + N 1s**가 부식 논거를 계산에서 실험으로 격상시킵니다.
  F4TCNQ **양성 대조군 없이는 해석 불가** — 거기서도 신호가 없으면 "깨끗하다"가 아니라
  "측정 감도 부족"입니다.
- ML은 넣지 마십시오. 데이터점 20여 개에 얹으면 상위 저널에서 감점입니다.
