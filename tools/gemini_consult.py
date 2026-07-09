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
    "flash": "gemini-2.5-flash",   # 빠르고 저렴 — 조사/초안/요약 오프로딩용
    "pro": "gemini-2.5-pro",       # 고성능 — 교차검증/어려운 판단용
    "flash2": "gemini-2.0-flash",
}

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
