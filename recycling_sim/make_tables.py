"""make_tables.py -- markdown tables from results/summary.json for the report."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(HERE, "results", "summary.json")))
D = S["designs"]

LABELS = {
    "ref_narrow528": "기준: 강공진 TEOLED (Ag/ITO 하부, thin Ag + NPB CPL), DE 최적",
    "ref_narrow528_fine": "기준: 강공진 TEOLED, 정밀격자 재최적",
    "ref_narrow528_ext": "기준: 강공진 TEOLED, 두께 범위 확장",
    "refdbr3_narrow528": "기준: 강공진 + 상부 DBR(3쌍) 캡",
    "refdbr3_narrow528_fine": "기준: 강공진 + 상부 DBR 캡, 정밀격자 재최적",
    "ref_ref_bulkAg": "기준: 강공진, thin Ag 벌크 광학상수(상한)",
    "prop3_narrow528": "제안: 사용자 구조 그대로 (HATCN/Ag 10 nm 실측 nk), DE 최적",
    "prop3_narrow528_ext": "제안: 동일, 두께 범위 확장",
    "prop5_narrow528": "제안: DBR 5쌍",
    "prop3_narrow528_Ag8": "제안: thin Ag 8 nm(seed 가정, nk는 HATCN 위 15 nm 실측값)",
    "prop3_narrow528_Ag12": "제안: thin Ag nk = 12 nm 실측(n=0.092)",
    "prop3_narrow528_Ag12_8": "제안: Ag12 nk + 7-8 nm",
    "prop3_narrow528_ODR": "제안: 반사판 아래 SiO2 800 nm (ODR)",
    "prop3_narrow528_lowkITO": "제안: ITO k ×0.3",
    "prop3_narrow528_IZO": "제안: 하부 IZO",
    "prop3_narrow528_lowHTL": "제안: 저굴절 HTL (n=1.55)",
    "prop3_narrow528_Ag12ODR": "제안: Ag12 nk + ODR",
    "prop3_narrow528_real": "제안: 개선 (Ag12 nk로 7-8 nm, ODR, low-k ITO) — 낙관적",
    "prop3_narrow528_cons": "제안: 개선 (HATCN 위 Ag nk로 7-8 nm, ODR, low-k ITO) — 보수적",
    "prop3_narrow528_real_lowHTL": "제안: 현실적 개선 + 저굴절 HTL",
    "prop3_narrow528_bulkAg": "제안: thin Ag 벌크 nk (도달 불가 상한)",
    "prop3_narrow528_best": "제안: 벌크 Ag + ODR + low-k ITO + 저굴절 HTL (상한)",
    "prop3_narrow528_IZOtop": "제안: 상부 양극을 무손상 IZO로 (가정) + ODR",
    "prop3_narrow528_asl1": "제안: SL 흡수 1%/왕복",
    "prop3_narrow528_asl3": "제안: SL 흡수 3%/왕복",
    "prophyb_narrow528": "하이브리드: 하부 TiO2/SiO2 DBR 2쌍 추가",
}


def row(t):
    d = D[t]
    lab = LABELS.get(t, t)
    cells = [lab, f"{d['CE']:.0f}", f"{100*d['eqe']:.1f}",
             f"{100*d.get('eqe_within_20deg', float('nan')):.1f}", f"{100*d.get('eqe_within_30deg', float('nan')):.1f}"]
    if "p_avg" in d:
        cells += [f"{d['p_avg']:.3f}", f"{d.get('Aprime_eff', float('nan')):.3f}", f"{d.get('T0_peak', float('nan')):.2f}",
                  f"{d.get('G_radiance_gain', float('nan')):.2f}"]
    else:
        cells += ["-", "-", "-", "-"]
    return "| " + " | ".join(cells) + " |"


def table(tags, title):
    print(f"\n**{title}**\n")
    print("| 설계 | 정면 cd/A | 총 EQE (%) | EQE ±20° (%) | EQE ±30° (%) | p | A′ | T₀ | G |")
    print("|---|---|---|---|---|---|---|---|---|")
    for t in tags:
        if t in D:
            print(row(t))


if __name__ == "__main__":
    green = [t for t in D if t.endswith("narrow528") or "narrow528_" in t and not t.endswith("_40")]
    order = ["ref_narrow528", "ref_narrow528_fine", "ref_narrow528_ext", "refdbr3_narrow528_fine", "refdbr3_narrow528",
             "ref_ref_bulkAg", "prop3_narrow528", "prop3_narrow528_ext", "prop5_narrow528", "prop3_narrow528_Ag8",
             "prop3_narrow528_Ag12", "prop3_narrow528_Ag12_8", "prop3_narrow528_ODR", "prop3_narrow528_lowkITO",
             "prop3_narrow528_IZO", "prop3_narrow528_lowHTL", "prop3_narrow528_Ag12ODR", "prop3_narrow528_cons",
             "prop3_narrow528_real",
             "prop3_narrow528_real_lowHTL", "prop3_narrow528_asl1", "prop3_narrow528_asl3", "prophyb_narrow528",
             "prop3_narrow528_bulkAg", "prop3_narrow528_best", "prop3_narrow528_IZOtop"]
    table(order, "녹색 협대역(528 nm / 25 nm) — 모든 설계")
    table(["ref_narrow528_40", "prop3_narrow528_40", "ref_irppy2acac", "ref_irppy2acac_DE", "prop3_irppy2acac"], "선폭 효과")
    table(["ref_blue_narrow_DE", "prop3_blue_narrow_DE", "prop3_blue_narrow_cons",
           "ref_vdabna_hf_DE", "prop3_vdabna_hf_DE", "prop3_vdabna_hf_best",
           "ref_red_narrow_DE", "ref_red_narrow_DEfine", "prop3_red_narrow_DE", "prop3_red_narrow_cons",
           "prop3_red_narrow_best"], "청색·적색")
    print("\n**정면 스펙트럼 CIE와 BI**\n")
    for t in sorted(D):
        if "CIE_onaxis" in D[t] and any(k in t for k in ("blue", "vdabna", "red", "narrow528_fine", "prop3_narrow528_ext", "cons")):
            x, y = D[t]["CIE_onaxis"]
            print(f"- {t}: CIE(x,y)=({x:.3f},{y:.3f}), cd/A={D[t]['CE']:.1f}, BI(cd/A/y)={D[t]['CE']/y:.0f}")

    print("\n**편광 재순환 (반사형 vs 흡수형 편광자, 정면 cd/A in one polarisation)**\n")
    for t in ["prop3_narrow528_ext", "prop3_narrow528_cons", "prop3_narrow528_real", "prop3_narrow528_best"]:
        if t in D and "gain" in D[t]:
            d = D[t]
            print(f"- {LABELS.get(t,t)}: 흡수형 {d['CE_pol_absorbing']:.0f} → 반사형 {d['CE_pol_reflective']:.0f} cd/A (×{d['gain']:.2f})")
    print("\n**손실 채널 분해 (기하 고정, 손실을 하나씩 제거)**\n")
    for t, rows in S.get("waterfall", {}).items():
        print(f"- {LABELS.get(t,t)}")
        for r in rows:
            print(f"  - {r['step']}: {r['CE']:.0f} cd/A, EQE {100*r['eqe']:.1f}%, p={r['p']:.3f}, a_top={r['a_top']:.3f}, R_b={r['rb']:.3f}")
    print("\n**MC: 산란층 파라미터와 블러 (재순환 성분)**\n")
    print("| 설계 | d (µm) | μs (1/µm) | g | μa (1/µm) | μs′d | η_rec 기여 EQE | 재순환 CE | r50 (µm) | r90 (µm) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for t, rows in S.get("mc_blur", {}).items():
        for r in rows:
            print(f"| {t} | {r['d_um']} | {r['mus']} | {r['g']} | {r['mua']} | {r['Sred']:.2f} | {100*r['eqe_rec']:.1f}% | "
                  f"{r['CE_rec_rad']:.0f} | {r['r50_um']:.1f} | {r['r90_um']:.1f} |")
