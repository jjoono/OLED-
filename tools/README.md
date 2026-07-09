# Gemini 협업(오케스트레이션) 도구

Claude(오케스트레이터)가 작업 중 Gemini를 **보조 엔진**으로 호출하기 위한 도구입니다.
목표는 두 가지입니다.

1. **토큰 절감** — 무거운 조사·초안·대용량 요약을 저렴한 `gemini-2.5-flash`로 오프로딩해,
   상위 오케스트레이터가 소비하는 토큰을 줄입니다.
2. **신뢰도 향상** — 중요한 판단은 `gemini-2.5-pro`로 교차검증(second opinion)하고,
   불일치가 나면 재검토합니다. 최종 사실 판단·통합은 Claude가 유지합니다.

## 사전 준비
- 환경 변수 `GEMINI_API_KEY` 필요 (환경 설정에 이미 등록됨). **코드/커밋에 키를 넣지 마세요.**

## 사용법
```bash
# 저렴한 모델로 조사/요약 오프로딩
echo "질문 또는 긴 텍스트" | python3 tools/gemini_consult.py --model flash

# 고성능 모델로 교차검증
python3 tools/gemini_consult.py --model pro \
  --system "너는 OLED 물리 전문가다" \
  --prompt "다음 주장을 사실검증: ..."

# 파일 통째로 넘겨 요약
python3 tools/gemini_consult.py --model flash --prompt-file draft.md
```
- 출력: 본문은 stdout, 토큰 사용량은 stderr(`[gemini:...] in=.. out=.. total=..`).
- `--json` 으로 원본 응답 확인 가능.

## 모델 별칭
| 별칭 | 실제 모델 | 용도 |
|------|-----------|------|
| `flash` | gemini-2.5-flash | 기본. 조사/초안/요약 오프로딩 |
| `pro` | gemini-2.5-pro | 교차검증/어려운 판단 |
| `flash2` | gemini-2.0-flash | 경량 대안 |

## 운영 원칙 (Claude ↔ Gemini 프로토콜)
1. **오프로딩**: 넓은 조사·긴 초안은 flash에 위임하고, Claude는 결과만 통합.
2. **교차검증**: 수치·사실 주장은 pro로 재확인. Gemini가 근거 없이 덧붙인 내용
   (예: 미확인 수치)은 Claude가 걸러냄 — Gemini 출력을 맹신하지 않음.
3. **키 보안**: 키는 환경 변수로만 참조. 로그/커밋/코드에 노출 금지.
