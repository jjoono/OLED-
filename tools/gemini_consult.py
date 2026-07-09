#!/usr/bin/env python3
"""
gemini_consult.py — Claude가 작업 중 Gemini를 보조 엔진으로 호출하기 위한 CLI.

목적:
  - 무거운 조사/초안/대용량 요약을 저렴한 모델(flash)로 오프로딩해 상위 오케스트레이터
    (Claude)의 토큰 소비를 줄인다.
  - 중요한 판단은 pro 모델로 교차검증(second opinion)해 결과 신뢰도를 높인다.

사용법:
  echo "질문" | python3 tools/gemini_consult.py --model flash
  python3 tools/gemini_consult.py --model pro --system "너는 OLED 물리 전문가다" --prompt "..."
  python3 tools/gemini_consult.py --model flash --prompt-file draft.md

환경:
  GEMINI_API_KEY 환경 변수 필요 (코드에 하드코딩 금지).
"""
import argparse
import json
import os
import sys
import time
import urllib.request
import urllib.error

MODEL_ALIASES = {
    # 최신 pro 우선. 무료 티어에서 pro가 429(quota)면 --fallback 로 하위 모델로 자동 강등.
    "pro": "gemini-3.1-pro-preview",   # 최신 고성능 (유료 티어 권장)
    "pro25": "gemini-2.5-pro",         # 이전 세대 pro
    "flash": "gemini-3.5-flash",       # 최신 flash
    "flash25": "gemini-2.5-flash",     # 무료 티어에서 가장 안정적으로 열려 있음
}

# --fallback 지정 시 429/503이 나면 순서대로 다음 모델을 시도한다.
FALLBACK_CHAIN = [
    "gemini-3.1-pro-preview",
    "gemini-2.5-pro",
    "gemini-2.5-flash",
]

ENDPOINT = "https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent"


def call_gemini(prompt, model="gemini-2.5-flash", system=None,
                temperature=0.4, max_output_tokens=8192, retries=3):
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        sys.exit("ERROR: GEMINI_API_KEY 환경 변수가 없습니다.")

    body = {
        "contents": [{"role": "user", "parts": [{"text": prompt}]}],
        "generationConfig": {
            "temperature": temperature,
            "maxOutputTokens": max_output_tokens,
        },
    }
    if system:
        body["systemInstruction"] = {"parts": [{"text": system}]}

    url = ENDPOINT.format(model=model)
    data = json.dumps(body).encode("utf-8")

    last_err = None
    for attempt in range(retries):
        req = urllib.request.Request(
            url, data=data, method="POST",
            headers={"Content-Type": "application/json", "x-goog-api-key": api_key},
        )
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            last_err = f"HTTP {e.code}: {e.read().decode('utf-8', 'ignore')[:500]}"
            if e.code in (429, 500, 503):  # 일시적 오류 → 백오프 재시도
                time.sleep(2 ** attempt)
                continue
            break
        except (urllib.error.URLError, TimeoutError) as e:
            last_err = str(e)
            time.sleep(2 ** attempt)
    sys.exit(f"ERROR: Gemini 호출 실패 — {last_err}")


def extract(resp):
    """응답에서 텍스트와 토큰 사용량을 뽑아낸다."""
    text_parts = []
    for cand in resp.get("candidates", []):
        for part in cand.get("content", {}).get("parts", []):
            if "text" in part:
                text_parts.append(part["text"])
    usage = resp.get("usageMetadata", {})
    return "".join(text_parts), usage


def main():
    ap = argparse.ArgumentParser(description="Gemini 보조 호출 CLI")
    ap.add_argument("--model", default="flash",
                    help="flash | pro | flash2 또는 정식 모델명 (기본 flash)")
    ap.add_argument("--system", default=None, help="시스템 지시문")
    ap.add_argument("--prompt", default=None, help="프롬프트 (미지정 시 stdin)")
    ap.add_argument("--prompt-file", default=None, help="프롬프트 파일 경로")
    ap.add_argument("--temperature", type=float, default=0.4)
    ap.add_argument("--max-tokens", type=int, default=8192)
    ap.add_argument("--fallback", action="store_true",
                    help="429/503 시 하위 모델로 자동 강등 (무료 티어 대비)")
    ap.add_argument("--json", action="store_true", help="원본 JSON 출력")
    args = ap.parse_args()

    if args.prompt_file:
        with open(args.prompt_file, "r", encoding="utf-8") as f:
            prompt = f.read()
    elif args.prompt:
        prompt = args.prompt
    else:
        prompt = sys.stdin.read()

    if not prompt.strip():
        sys.exit("ERROR: 프롬프트가 비어 있습니다.")

    model = MODEL_ALIASES.get(args.model, args.model)

    if args.fallback:
        # 지정 모델부터 시작해, 429/503이면 체인의 다음(하위) 모델로 자동 강등.
        chain = [model] + [m for m in FALLBACK_CHAIN if m != model]
        resp = None
        for m in chain:
            try:
                resp = call_gemini(prompt, model=m, system=args.system,
                                   temperature=args.temperature,
                                   max_output_tokens=args.max_tokens, retries=1)
                model = m
                break
            except SystemExit as e:
                if "429" in str(e) or "503" in str(e):
                    sys.stderr.write(f"[fallback] {m} 사용 불가 → 다음 모델 시도\n")
                    continue
                raise
        if resp is None:
            sys.exit("ERROR: 폴백 체인의 모든 모델이 사용 불가(quota/unavailable).")
    else:
        resp = call_gemini(prompt, model=model, system=args.system,
                           temperature=args.temperature, max_output_tokens=args.max_tokens)

    if args.json:
        print(json.dumps(resp, ensure_ascii=False, indent=2))
        return

    text, usage = extract(resp)
    print(text.rstrip())
    # 토큰 사용량은 stderr로 (파이프 오염 방지)
    if usage:
        sys.stderr.write(
            f"\n[gemini:{model}] "
            f"in={usage.get('promptTokenCount','?')} "
            f"out={usage.get('candidatesTokenCount','?')} "
            f"total={usage.get('totalTokenCount','?')} tokens\n"
        )


if __name__ == "__main__":
    main()
