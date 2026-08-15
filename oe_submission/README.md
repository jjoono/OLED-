# Optics Express 투고 패키지

`md2optica.py`가 `manuscript_draft.md` / `supplementary.md`에서 생성한다.
**원고를 고칠 때는 Markdown을 고치고 변환기를 다시 돌릴 것** — .tex 직접 수정은
다음 재생성 때 사라진다.

```
python3 md2optica.py
```

## 파일

| 파일 | 내용 |
|---|---|
| `main.tex` | 본문 (opticajournal.cls 대상) |
| `supplement.tex` | Supplement 1 (독립 article 클래스 — 바로 컴파일 가능) |
| `figures/` | Fig. 1–5, S1–S4 (PDF, 300 dpi 소스에서 생성) |

## 컴파일 방법

`opticajournal.cls`는 Optica 저작권 파일이라 레포에 넣지 않았다. 두 가지 방법:

1. **Overleaf (권장)**: New Project → Templates → "Optica journals" 공식 템플릿
   생성 → 템플릿의 main.tex를 이 `main.tex`로 교체, `figures/` 업로드.
2. 로컬: https://opg.optica.org/submit/templates.cfm 에서 스타일 패키지를 받아
   `opticajournal.cls`를 이 디렉토리에 두고 `pdflatex main` ×2.

`supplement.tex`는 표준 article 클래스라 그대로 `pdflatex supplement` ×2.

## 투고 전 남은 작업 (사람 몫)

1. **저자/소속/이메일** — `main.tex` 상단 `\author` 블록 (마커: `[Affiliation`,
   `[corresponding@email]`)과 Funding/Acknowledgments 자리표시.
2. 컴파일 후 그림 배치/줄바꿈 확인 (figure* 3개: Fig. 1, 2, 5).

완료된 것: 초록은 ~214단어로 축약 완료(장문판은 git 이력에 있음), 참고문헌
17건 전부 Optica 형식("A. B. Author, ``Title,'' J. Abbrev. **vol**, pages
(year). DOI")으로 자동 변환됨 — 서식 규칙은 `md2optica.py`에 있으므로 재생성해도
유지된다.

## 변환 규칙 요약

- `**bold**`→`\textbf`, `` `code` ``→`\texttt`(밑줄 이스케이프), `[n]`→`\cite{rn}`,
  "Fig. N"→`Fig.~\ref{fig:N}` (ref. [17] 자체 그림 인용은 보호됨)
- 유니코드(°, ×, ±, μ, ~)는 `\ensuremath` 매크로로 변환 — pdflatex 안전
- 그림 캡션 블록은 해당 위치의 `figure`/`figure*` 환경으로 변환
- back matter는 `\begin{backmatter}` + `\bmsection` (OE 필수 섹션 구조)
