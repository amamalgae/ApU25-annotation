#!/usr/bin/env python3
"""InterPro Matches API lookup for the UTEX 25 projected proteome (step 0).

The API returns pre-computed InterProScan 6 results for any sequence already in
UniParc, keyed by the MD5 of the sequence.  Sequences that miss are the input to
a local InterProScan 6 run.

MD5 convention (from the API README): hash the sequence uppercased and with the
trailing ``*`` removed, then send the digest as uppercase hex.

Subcommands, in the order they are meant to be run:

  md5     offline; writes md5_table.tsv.  No network.
  probe   finds the endpoint that answers 200 and dumps ONE raw response so the
          JSON schema can be read off it before any aggregation is written.
  fetch   runs all batches, writes matches_raw.json.gz.
  report  turns matches_raw.json.gz into lookup_status.tsv + unmatched.faa and
          prints the numbers SUMMARY.md needs.

`fetch` paces itself below 3 req/s and retries 429/5xx with exponential backoff.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import statistics
import sys
import time
from pathlib import Path

import requests

BATCH_SIZE = 100          # API limit: 100 MD5 per request
MIN_INTERVAL = 0.40       # seconds between request starts -> 2.5 req/s < 3 req/s
MAX_RETRIES = 4           # on 429/5xx, backing off 2s, 4s, 8s, 16s

# (a) in the brief: the public EBI path is not documented in the API README,
# which only covers self-hosting.  These are tried in order; if none answers,
# `probe` falls back to reading /openapi.json and /docs for the real path.
CANDIDATE_ENDPOINTS = [
    "https://www.ebi.ac.uk/interpro/matches/api/matches",
    "https://www.ebi.ac.uk/interpro/matches/api/",
]
DISCOVERY_URLS = [
    "https://www.ebi.ac.uk/interpro/matches/api/openapi.json",
    "https://www.ebi.ac.uk/interpro/matches/api/docs",
]

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FAA = ROOT / "annotation" / "UTEX25_proteins.faa"
OUT_DIR = ROOT / "step0_out"


# --------------------------------------------------------------------------- #
# FASTA / MD5
# --------------------------------------------------------------------------- #

def read_fasta(path: Path):
    """Yield (header, sequence) preserving file order."""
    header, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header, buf = line[1:], []
            else:
                buf.append(line.strip())
    if header is not None:
        yield header, "".join(buf)


def seq_md5(seq: str) -> tuple[str, str]:
    """Return (normalised sequence, uppercase hex MD5)."""
    core = seq.upper().rstrip("*")
    return core, hashlib.md5(core.encode()).hexdigest().upper()


def load_records(faa: Path):
    """[(protein_id, header, sequence, md5, has_internal_stop)] in file order."""
    records = []
    for header, seq in read_fasta(faa):
        core, digest = seq_md5(seq)
        records.append((header.split()[0], header, core, digest, "*" in core))
    return records


# --------------------------------------------------------------------------- #
# HTTP
# --------------------------------------------------------------------------- #

def post_batch(session, url, md5s, *, timeout=120):
    """POST one batch, retrying 429/5xx.  Returns the requests.Response."""
    delay = 2.0
    for attempt in range(MAX_RETRIES + 1):
        resp = session.post(
            url,
            json={"md5": md5s},
            headers={"Content-Type": "application/json",
                     "Accept": "application/json"},
            timeout=timeout,
        )
        if resp.status_code < 400:
            return resp
        retryable = resp.status_code == 429 or 500 <= resp.status_code < 600
        if not retryable or attempt == MAX_RETRIES:
            return resp
        wait = delay
        retry_after = resp.headers.get("Retry-After")
        if retry_after:
            try:
                wait = max(wait, float(retry_after))
            except ValueError:
                pass
        print(f"  HTTP {resp.status_code}; retry {attempt + 1}/{MAX_RETRIES} "
              f"in {wait:.0f}s", file=sys.stderr)
        time.sleep(wait)
        delay *= 2
    return resp  # unreachable


# --------------------------------------------------------------------------- #
# Response schema
# --------------------------------------------------------------------------- #

def resolve_by_md5(body):
    """Map MD5 -> payload for one response body.

    The response schema is item (b) of the brief and is NOT documented in the
    API README, so this accepts the handful of shapes a FastAPI service of this
    kind plausibly returns and raises loudly on anything else rather than
    guessing.  Confirm it against probe_response.json before trusting a run.
    """
    def from_list(items):
        out = {}
        for item in items:
            if not isinstance(item, dict):
                raise TypeError(f"list element is {type(item).__name__}, not object")
            key = next((item[k] for k in ("md5", "MD5", "sequence_md5", "id")
                        if isinstance(item.get(k), str)), None)
            if key is None:
                raise KeyError(f"no md5-like key in result object: {sorted(item)}")
            out[key.upper()] = item
        return out

    if isinstance(body, list):
        return from_list(body)

    if isinstance(body, dict):
        if not body:
            return {}  # a batch in which no MD5 was found
        # {"MD5HEX": {...}, ...}
        if all(re.fullmatch(r"[0-9A-Fa-f]{32}", k) for k in body):
            return {k.upper(): v for k, v in body.items()}
        for key in ("results", "matches", "data", "items", "sequences"):
            if key in body:
                inner = body[key]
                if isinstance(inner, list):
                    return from_list(inner)
                if isinstance(inner, dict):
                    return {k.upper(): v for k, v in inner.items()}

    raise ValueError(
        "unrecognised response schema; top-level "
        f"{type(body).__name__} with keys "
        f"{sorted(body)[:20] if isinstance(body, dict) else 'n/a'}. "
        "Inspect step0_out/probe_response.json and extend resolve_by_md5()."
    )


ACCESSION_PATTERNS = {
    "interpro": re.compile(r"\bIPR\d{6}\b"),
    "pfam": re.compile(r"\bPF\d{5}\b"),
    "go": re.compile(r"\bGO:\d{7}\b"),
}


def count_accessions(payload):
    """Distinct InterPro / Pfam / GO accessions anywhere in one MD5's payload.

    Deliberately schema-agnostic: it walks every string in the payload and
    matches accession syntax, so it stays correct regardless of where the API
    nests them.
    """
    found = {k: set() for k in ACCESSION_PATTERNS}
    stack = [payload]
    while stack:
        node = stack.pop()
        if isinstance(node, str):
            for name, pattern in ACCESSION_PATTERNS.items():
                found[name].update(pattern.findall(node))
        elif isinstance(node, dict):
            stack.extend(node.keys())
            stack.extend(node.values())
        elif isinstance(node, list):
            stack.extend(node)
    return {k: len(v) for k, v in found.items()}


# Keys that only echo the request back, carrying no result content.
ECHO_KEYS = {"md5", "MD5", "sequence_md5", "id"}


def is_hit(payload) -> bool:
    """Whether an entry returned for an MD5 means "this sequence is in UniParc".

    The API is keyed on UniParc membership, so the presence of an entry is
    taken as the hit signal and a sequence that is in UniParc but has no
    signature matches still counts as a hit (in_uniparc=yes, 0 InterPro
    entries).  Explicitly empty payloads, and entries that only echo the
    submitted MD5 back with no result fields, count as misses.

    UNCONFIRMED: whether the API omits absent MD5 entirely or returns an empty
    entry for them is item (b) of the brief and needs checking against
    step0_out/probe_response.json before a run is trusted.
    """
    if payload is None:
        return False
    if isinstance(payload, (list, dict)) and len(payload) == 0:
        return False
    if isinstance(payload, dict) and not (set(payload) - ECHO_KEYS):
        return False
    return True


# --------------------------------------------------------------------------- #
# Subcommands
# --------------------------------------------------------------------------- #

def cmd_md5(args):
    records = load_records(args.faa)
    OUT_DIR.mkdir(exist_ok=True)
    path = OUT_DIR / "md5_table.tsv"
    with open(path, "w") as fh:
        fh.write("protein_id\tmd5\tlength\thas_internal_stop\n")
        for pid, _hdr, core, digest, stop in records:
            fh.write(f"{pid}\t{digest}\t{len(core)}\t{'yes' if stop else 'no'}\n")
    unique = {r[3] for r in records}
    print(f"fasta           : {args.faa}")
    print(f"sequences       : {len(records)}")
    print(f"unique md5      : {len(unique)}")
    print(f"duplicate seqs  : {len(records) - len(unique)}")
    print(f"internal stops  : {sum(1 for r in records if r[4])}")
    print(f"requests needed : {math.ceil(len(unique) / BATCH_SIZE)}")
    print(f"wrote           : {path}")


def cmd_probe(args):
    records = load_records(args.faa)
    probe_md5 = next(r[3] for r in records if not r[4])  # a clean sequence
    OUT_DIR.mkdir(exist_ok=True)
    session = requests.Session()

    for url in CANDIDATE_ENDPOINTS:
        print(f"POST {url}", file=sys.stderr)
        try:
            resp = session.post(url, json={"md5": [probe_md5]}, timeout=60,
                                headers={"Content-Type": "application/json"})
        except requests.RequestException as exc:
            print(f"  request failed: {exc}", file=sys.stderr)
            continue
        print(f"  -> HTTP {resp.status_code}", file=sys.stderr)
        if resp.status_code != 200:
            continue
        body = resp.json()
        (OUT_DIR / "probe_response.json").write_text(
            json.dumps({"endpoint": resp.url, "md5": probe_md5, "body": body},
                       indent=2, ensure_ascii=False))
        print(f"\nendpoint returning 200 : {resp.url}")
        print(f"top-level type         : {type(body).__name__}")
        if isinstance(body, dict):
            print(f"top-level keys         : {sorted(body)}")
        first = body[0] if isinstance(body, list) and body else body
        if isinstance(first, dict):
            print(f"first element keys     : {sorted(first)}")
        print(f"raw response saved to  : {OUT_DIR / 'probe_response.json'}")
        return 0

    print("\nNo candidate endpoint answered 200. Trying schema discovery:",
          file=sys.stderr)
    for url in DISCOVERY_URLS:
        try:
            resp = session.get(url, timeout=60)
        except requests.RequestException as exc:
            print(f"  GET {url} failed: {exc}", file=sys.stderr)
            continue
        print(f"  GET {url} -> HTTP {resp.status_code}", file=sys.stderr)
        if resp.status_code == 200:
            dest = OUT_DIR / ("openapi.json" if url.endswith(".json") else "docs.html")
            dest.write_bytes(resp.content)
            print(f"  saved {dest}; read the real path out of it", file=sys.stderr)
    return 1


def cmd_fetch(args):
    records = load_records(args.faa)
    md5s = sorted({r[3] for r in records})
    batches = [md5s[i:i + BATCH_SIZE] for i in range(0, len(md5s), BATCH_SIZE)]
    print(f"{len(md5s)} unique MD5 in {len(batches)} requests to {args.endpoint}",
          file=sys.stderr)

    OUT_DIR.mkdir(exist_ok=True)
    session = requests.Session()
    collected, failures = [], 0
    last = 0.0

    for i, batch in enumerate(batches, 1):
        gap = MIN_INTERVAL - (time.monotonic() - last)
        if gap > 0:
            time.sleep(gap)
        last = time.monotonic()
        resp = post_batch(session, args.endpoint, batch)
        entry = {"batch": i, "md5": batch, "status": resp.status_code}
        if resp.status_code == 200:
            entry["body"] = resp.json()
        else:
            entry["error"] = resp.text[:2000]
            failures += 1
            print(f"batch {i}: HTTP {resp.status_code} (gave up)", file=sys.stderr)
        collected.append(entry)
        if i % 10 == 0 or i == len(batches):
            print(f"  {i}/{len(batches)} batches", file=sys.stderr)

    out = OUT_DIR / "matches_raw.json.gz"
    with gzip.open(out, "wt", encoding="utf-8") as fh:
        json.dump({"endpoint": args.endpoint,
                   "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   "batch_size": BATCH_SIZE,
                   "n_unique_md5": len(md5s),
                   "batches": collected}, fh, ensure_ascii=False)
    print(f"wrote {out}; {failures} failed batches", file=sys.stderr)
    return 1 if failures else 0


def quartiles(values):
    if not values:
        return (0.0, 0.0, 0.0)
    ordered = sorted(values)
    if len(ordered) == 1:
        v = float(ordered[0])
        return (v, v, v)
    q1, med, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
    return (q1, med, q3)


def cmd_report(args):
    with gzip.open(OUT_DIR / "matches_raw.json.gz", "rt", encoding="utf-8") as fh:
        raw = json.load(fh)

    by_md5, bad = {}, 0
    for entry in raw["batches"]:
        if entry.get("status") != 200:
            bad += 1
            continue
        by_md5.update(resolve_by_md5(entry["body"]))
    if bad:
        print(f"WARNING: {bad} batches did not return 200; "
              "their MD5 are reported as not looked up", file=sys.stderr)

    records = load_records(args.faa)
    counts = {"interpro": [], "pfam": [], "go": []}
    n_hit = n_hit_clean = 0
    lookup = OUT_DIR / "lookup_status.tsv"
    unmatched = OUT_DIR / "unmatched.faa"

    with open(lookup, "w") as tsv, open(unmatched, "w") as faa:
        tsv.write("protein_id\tmd5\tin_uniparc\n")
        for pid, header, core, digest, stop in records:
            hit = is_hit(by_md5.get(digest))
            tsv.write(f"{pid}\t{digest}\t{'yes' if hit else 'no'}\n")
            if hit:
                n_hit += 1
                if not stop:
                    n_hit_clean += 1
                for key, value in count_accessions(by_md5[digest]).items():
                    counts[key].append(value)
            else:
                faa.write(f">{header}\n")
                for j in range(0, len(core), 60):
                    faa.write(core[j:j + 60] + "\n")

    total = len(records)
    n_stop = sum(1 for r in records if r[4])
    clean = total - n_stop
    print(f"sequences                    : {total}")
    print(f"internal-stop sequences      : {n_stop}")
    print(f"UniParc hits                 : {n_hit} ({100 * n_hit / total:.2f}% of {total})")
    print(f"UniParc hits, clean subset   : {n_hit_clean} "
          f"({100 * n_hit_clean / clean:.2f}% of {clean})")
    print(f"unmatched                    : {total - n_hit}")
    for key in ("interpro", "pfam", "go"):
        q1, med, q3 = quartiles(counts[key])
        print(f"{key:<12} per hit (Q1/med/Q3): {q1:.1f} / {med:.1f} / {q3:.1f}")
    print(f"wrote {lookup} and {unmatched}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--faa", type=Path, default=DEFAULT_FAA)
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("md5").set_defaults(func=cmd_md5)
    sub.add_parser("probe").set_defaults(func=cmd_probe)
    fetch = sub.add_parser("fetch")
    fetch.add_argument("--endpoint", required=True,
                       help="the URL confirmed by `probe`")
    fetch.set_defaults(func=cmd_fetch)
    sub.add_parser("report").set_defaults(func=cmd_report)
    args = parser.parse_args()
    sys.exit(args.func(args) or 0)


if __name__ == "__main__":
    main()
