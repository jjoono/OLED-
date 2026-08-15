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
2. **초록 축약** — 현재 초록은 원고 전체 요약이라 매우 길다(~600단어). OE는
   형식 제한은 없으나 관행상 ~200단어를 권장. 축약본을 만들면 Markdown 쪽
   Abstract를 바꾸고 재생성.
3. **참고문헌 서식** — 현재는 Markdown의 세미콜론 스타일을 그대로 옮겼다
   (성, 약자 순서). Optica 스타일은 "A. B. Author, ..." 순서라 cosmetic 정리
   필요 — 내용(권·쪽·DOI)은 전부 검증 완료 상태.
4. 컴파일 후 그림 배치/줄바꿈 확인 (figure* 3개: Fig. 1, 2, 5).

## 변환 규칙 요약

- `**bold**`→`\textbf`, `` `code` ``→`\texttt`(밑줄 이스케이프), `[n]`→`\cite{rn}`,
  "Fig. N"→`Fig.~\ref{fig:N}` (ref. [17] 자체 그림 인용은 보호됨)
- 유니코드(°, ×, ±, μ, ~)는 `\ensuremath` 매크로로 변환 — pdflatex 안전
- 그림 캡션 블록은 해당 위치의 `figure`/`figure*` 환경으로 변환
- back matter는 `\begin{backmatter}` + `\bmsection` (OE 필수 섹션 구조)
