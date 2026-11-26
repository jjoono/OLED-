# Optics Express 투고 패키지

`md2optica.py`가 `manuscript_draft.md`와 `supplementary.md`에서 생성한다. 원고를 고칠 때는 Markdown을 고치고 변환기를 다시 돌린다(.tex를 직접 고치면 다음 생성 때 사라진다).

```
python3 tools/cite.py tools/manuscript_en.src.md manuscript_draft.md   # 인용 번호·참고문헌 목록 생성
python3 md2optica.py                                                  # main.tex, supplement.tex, figures/
```

## 파일

| 파일 | 내용 |
|---|---|
| `main.tex` | 본문. Optica 공용 템플릿(`optica-article.cls`, `\journal{opticajournal}`, `\articletype{Research Article}`) 대상 |
| `supplement.tex` | Supplement 1. 표준 article 클래스로 바로 컴파일된다 |
| `figures/` | 본문 Fig. 1–5와 Fig. S1–S10의 PDF(벡터, 글꼴 내장). 본문 그림 폭은 5.25 in |

## 컴파일

`optica-article.cls`와 `opticajnl.bst`는 Optica 배포 파일이라 저장소에 넣지 않았다.

1. Overleaf: Optica 공용 템플릿("Optica Publishing Group journal template")으로 프로젝트를 만들고, 템플릿의 본문 파일을 이 `main.tex`로 바꾼 뒤 `figures/`를 올린다.
2. 로컬: Optica Author Resources의 LaTeX 템플릿 묶음에서 클래스 파일을 받아 이 폴더에 두고 `pdflatex main`을 두 번 실행한다.

`supplement.tex`는 `pdflatex supplement`를 두 번 실행한다.

## 양식 점검 결과

- 제목은 문장형 대소문자, 초록은 150단어, 수식은 (1)–(10)으로 번호를 붙였고 본문에서 "Eq. (n)"으로 인용한다.
- 그림은 "Fig. 1(a)" 형식으로 인용하고 번호 순서대로 처음 인용된다(Fig. 1–5, Table 1–3).
- 참고문헌 33건은 Optica 양식(저자 3인 + et al., "제목," 학술지 **권**, 쪽 (연도))이며 첫 인용 순서로 번호가 매겨진다. 검증 기록은 `review/reference_check.md`.
- 뒷부분은 Funding, Acknowledgment, Disclosures, Data availability, Supplemental document 순서이다.

## 투고 전에 사람이 채울 것

1. 저자, 소속, 교신 이메일(`main.tex`의 `\author`, `\address`, `\email`).
2. Funding과 Acknowledgment의 `[...]` 자리.
3. Data availability 문구. 지금은 "요청 시 제공"으로 되어 있다. 자료를 공개 저장소나 Optica의 Dataset으로 올리면 "Data underlying the results presented in this paper are available in Dataset 1, Ref. [n]." 형식으로 바꾸고 참고문헌에 추가한다.
4. 참고문헌 [14](arXiv 사전공개본)의 출판 여부.
5. 커버레터(`review/cover_letter_oe.md`)의 추천 심사위원과 관련 논문 문단.
