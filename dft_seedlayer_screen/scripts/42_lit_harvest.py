"""Systematic literature harvest for the Ag seed-layer novelty search.

Run this ON YOUR OWN MACHINE. The cloud session it was written in cannot reach
any scholarly API (the network policy answers 403 to openalex, crossref,
semantic scholar, unpaywall and arxiv alike), which is exactly why the earlier
keyword search missed Park & Suh 2018.

    python 42_lit_harvest.py

Everything is configured in CONFIG below. Only `requests` is needed:
    pip install requests

WHAT IT DOES
  1. OpenAlex keyword search over the query list.
  2. Citation expansion around the seed papers -- every work that CITES them and
     every work they REFERENCE. This is the part plain keyword search cannot do,
     and it is where the misses come from: the HATCN/Ag literature is written in
     "hole injection layer" language, so "seed layer" queries never reach it.
  3. Unpaywall lookup for a LEGAL open-access copy of each hit.
  4. Downloads the open-access PDFs into pdfs/.
  5. Writes results.csv (everything found, with a relevance flag) and
     manual_fetch.html -- a clickable list of the paywalled ones for you to grab
     through your institutional access.

WHAT IT DELIBERATELY DOES NOT DO
  It does not log in to publishers or pull paywalled PDFs through an
  institutional proxy. That breaches publisher terms and gets university IP
  ranges cut off. The paywalled items are listed for manual download instead;
  in this field expect roughly half the hits to be open access.

POLITENESS
  OpenAlex and Unpaywall both ask for an email to put you in their fast "polite"
  pool. Put a real one in CONFIG. The script sleeps between calls; do not remove
  that -- these are free services.
"""
import os, re, csv, json, time, html
import requests

# ============================== CONFIG =====================================
CONFIG = {
    # REQUIRED: your real email. OpenAlex/Unpaywall use it for the polite pool.
    "email": "jhk4733@gmail.com",

    "outdir": "lit_harvest",

    # Papers whose citation neighbourhood we want mapped. These are the ones
    # already known to matter; the harvest grows outward from them.
    "seed_dois": [
        "10.1364/OE.26.004979",       # Park & Suh 2018, HATCN(7)/Ag(15-25), OE
        "10.1002/adfm.201502542",     # Kim 2015, ZnS/Cs2CO3/Ag TFOLED, AFM
        "10.1002/adom.201300241",     # Schwab 2013, ultra-thin wetting layer top electrode
        "10.1016/j.orgel.2013.09.010",# ZnS/Ag/MoO3 "quasi-perfect" Ag (verify this DOI)
    ],

    # Keyword queries. Deliberately mixes the "seed/wetting" vocabulary with the
    # "injection layer" vocabulary, because the relevant prior art uses both.
    "queries": [
        "HATCN silver electrode",
        "hexaazatriphenylene hexacarbonitrile silver",
        "ultrathin silver transparent electrode seed layer",
        "silver wetting layer organic light emitting diode",
        "nucleation layer ultrathin silver film percolation",
        "dielectric metal dielectric transparent electrode OLED",
        "top emitting OLED semitransparent silver anode",
        "silver film organic semiconductor Volmer-Weber growth",
        "seed layer thin silver sheet resistance transmittance",
        "hole injection layer thin silver top electrode",
    ],

    # Words that mark a hit as worth reading first.
    "priority_terms": [
        "hatcn", "hexaazatriphenylene", "seed", "wetting", "nucleation",
        "percolation", "ultrathin", "ultra-thin", "coverage", "volmer",
    ],

    "max_per_query": 200,        # OpenAlex pages of 200
    "download_oa": True,
    "sleep": 0.4,                # seconds between API calls
}
# ===========================================================================

OA = "https://api.openalex.org"
UP = "https://api.unpaywall.org/v2"


def _get(url, **params):
    params.setdefault("mailto", CONFIG["email"])
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=40)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503):
                time.sleep(2 ** attempt)
                continue
            return None
        except requests.RequestException:
            time.sleep(2 ** attempt)
    return None


def norm_doi(d):
    if not d:
        return None
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", d.strip()).lower()


def record(w):
    """Flatten one OpenAlex work."""
    inv = w.get("abstract_inverted_index") or {}
    abstract = ""
    if inv:
        pos = {}
        for word, idxs in inv.items():
            for i in idxs:
                pos[i] = word
        abstract = " ".join(pos[i] for i in sorted(pos))
    pl = (w.get("primary_location") or {})
    src = (pl.get("source") or {})
    return {
        "doi": norm_doi(w.get("doi")),
        "year": w.get("publication_year"),
        "title": (w.get("title") or "").replace("\n", " "),
        "journal": src.get("display_name") or "",
        "authors": "; ".join(
            (a.get("author") or {}).get("display_name", "")
            for a in (w.get("authorships") or [])[:6]),
        "cited_by": w.get("cited_by_count", 0),
        "is_oa": bool((w.get("open_access") or {}).get("is_oa")),
        "openalex": w.get("id", ""),
        "abstract": abstract[:1500],
        "found_via": "",
    }


def search(query, found):
    print(f"  query: {query}", flush=True)
    cursor, n = "*", 0
    while cursor and n < CONFIG["max_per_query"]:
        d = _get(f"{OA}/works", search=query, per_page=100, cursor=cursor)
        time.sleep(CONFIG["sleep"])
        if not d:
            break
        for w in d.get("results", []):
            r = record(w)
            if not r["doi"]:
                continue
            n += 1
            if r["doi"] not in found:
                r["found_via"] = f"query:{query[:30]}"
                found[r["doi"]] = r
        cursor = (d.get("meta") or {}).get("next_cursor")
    print(f"    -> {n} works", flush=True)


def expand_citations(doi, found):
    """Everything that cites this work, and everything it references."""
    w = _get(f"{OA}/works/https://doi.org/{doi}")
    time.sleep(CONFIG["sleep"])
    if not w:
        print(f"  [!] seed not found: {doi}", flush=True)
        return
    title = (w.get("title") or "")[:60]
    print(f"  seed: {doi}  {title}", flush=True)

    # citing works
    cursor, n = "*", 0
    while cursor:
        d = _get(f"{OA}/works", filter=f"cites:{w['id'].split('/')[-1]}",
                 per_page=100, cursor=cursor)
        time.sleep(CONFIG["sleep"])
        if not d:
            break
        for x in d.get("results", []):
            r = record(x)
            if r["doi"] and r["doi"] not in found:
                r["found_via"] = f"cites:{doi}"
                found[r["doi"]] = r
                n += 1
        cursor = (d.get("meta") or {}).get("next_cursor")
    print(f"    citing: {n} new", flush=True)

    # referenced works
    refs = w.get("referenced_works") or []
    m = 0
    for i in range(0, len(refs), 50):
        ids = "|".join(u.split("/")[-1] for u in refs[i:i + 50])
        d = _get(f"{OA}/works", filter=f"openalex_id:{ids}", per_page=100)
        time.sleep(CONFIG["sleep"])
        if not d:
            continue
        for x in d.get("results", []):
            r = record(x)
            if r["doi"] and r["doi"] not in found:
                r["found_via"] = f"ref_of:{doi}"
                found[r["doi"]] = r
                m += 1
    print(f"    referenced: {m} new", flush=True)


def oa_pdf_url(doi):
    d = _get(f"{UP}/{doi}", email=CONFIG["email"])
    time.sleep(CONFIG["sleep"])
    if not d or not d.get("is_oa"):
        return None
    best = d.get("best_oa_location") or {}
    return best.get("url_for_pdf") or best.get("url")


def safe_name(r):
    a = (r["authors"].split(";")[0] or "unknown").split()[-1]
    t = re.sub(r"[^A-Za-z0-9]+", "_", r["title"])[:60].strip("_")
    return f"{r['year']}_{a}_{t}.pdf"


def main():
    out = CONFIG["outdir"]
    os.makedirs(os.path.join(out, "pdfs"), exist_ok=True)
    if "@" not in CONFIG["email"]:
        raise SystemExit("set a real email in CONFIG['email'] first")

    found = {}
    print("== keyword search ==")
    for q in CONFIG["queries"]:
        search(q, found)
    print("\n== citation expansion ==")
    for doi in CONFIG["seed_dois"]:
        expand_citations(doi, found)

    print(f"\n{len(found)} unique works")

    # flag the ones to read first
    for r in found.values():
        blob = (r["title"] + " " + r["abstract"]).lower()
        hits = [t for t in CONFIG["priority_terms"] if t in blob]
        ag = ("silver" in blob or " ag " in blob or "ag/" in blob)
        r["priority"] = "HIGH" if (ag and len(hits) >= 2) else (
            "MED" if (ag and hits) else "low")
        r["matched_terms"] = ",".join(hits)

    order = {"HIGH": 0, "MED": 1, "low": 2}
    rows = sorted(found.values(),
                  key=lambda r: (order[r["priority"]], -(r["cited_by"] or 0)))

    print("== unpaywall lookup + download ==")
    n_dl = 0
    for r in rows:
        r["oa_pdf"] = ""
        r["local_pdf"] = ""
        if r["priority"] == "low":
            continue                       # only chase the ones worth reading
        url = oa_pdf_url(r["doi"])
        if not url:
            continue
        r["oa_pdf"] = url
        if not CONFIG["download_oa"]:
            continue
        path = os.path.join(out, "pdfs", safe_name(r))
        if os.path.exists(path):
            r["local_pdf"] = path
            continue
        try:
            resp = requests.get(url, timeout=60,
                                headers={"User-Agent": f"lit-harvest ({CONFIG['email']})"})
            if resp.status_code == 200 and resp.content[:4] == b"%PDF":
                open(path, "wb").write(resp.content)
                r["local_pdf"] = path
                n_dl += 1
                print(f"  [pdf] {os.path.basename(path)}", flush=True)
        except requests.RequestException:
            pass
        time.sleep(CONFIG["sleep"])

    cols = ["priority", "year", "title", "journal", "authors", "cited_by",
            "doi", "is_oa", "oa_pdf", "local_pdf", "matched_terms",
            "found_via", "openalex", "abstract"]
    with open(os.path.join(out, "results.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    todo = [r for r in rows if r["priority"] in ("HIGH", "MED") and not r["local_pdf"]]
    with open(os.path.join(out, "manual_fetch.html"), "w", encoding="utf-8") as f:
        f.write("<meta charset='utf-8'><style>body{font:14px/1.5 sans-serif;max-width:60em;"
                "margin:2em auto}li{margin:.6em 0}.h{color:#b00;font-weight:700}</style>")
        f.write(f"<h2>Manual fetch: {len(todo)} papers</h2>")
        f.write("<p>Open through your institutional access and save the PDF into "
                "<code>pdfs/</code>.</p><ol>")
        for r in todo:
            f.write(f"<li><span class='{'h' if r['priority']=='HIGH' else ''}'>"
                    f"[{r['priority']}]</span> {html.escape(r['title'])}<br>"
                    f"<small>{html.escape(r['journal'])} {r['year']} &middot; "
                    f"cited {r['cited_by']} &middot; via {html.escape(r['found_via'])}</small><br>"
                    f"<a href='https://doi.org/{r['doi']}'>{r['doi']}</a></li>")
        f.write("</ol>")

    hi = sum(1 for r in rows if r["priority"] == "HIGH")
    md = sum(1 for r in rows if r["priority"] == "MED")
    print(f"\n  {len(rows)} works   HIGH {hi}   MED {md}")
    print(f"  downloaded {n_dl} open-access PDFs -> {out}/pdfs/")
    print(f"  {len(todo)} need manual fetch    -> open {out}/manual_fetch.html")
    print(f"  full table                       -> {out}/results.csv")
    print("\nUpload results.csv back and I will classify the prior art from it.")


if __name__ == "__main__":
    main()
