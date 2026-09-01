#!/usr/bin/env python3
"""Check whether an InterPro release matches the versions the Matches API returned.

The Matches API gives no release number of its own; it reports a version per
signature library.  InterProScan 6 ships one ``databases.json`` per InterPro
release listing exactly those library versions, so comparing the two says which
release ``--interpro`` has to be pinned to for the API results and a local run to
be annotated against the same basis.

``raw/interpro-109.0-databases.json`` is that file, taken from the data the
pipeline downloads:

    nextflow run <interproscan6> -profile docker,test --datadir data \\
        --interpro latest
    cp data/interpro/109.0/databases.json raw/interpro-109.0-databases.json

    python3 scripts/08_ips6_versions.py
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MATCHES = ROOT / "step0_out" / "matches_raw.json.gz"
DEFAULT_DATABASES = ROOT / "raw" / "interpro-109.0-databases.json"


def normalise(name: str) -> str:
    """InterProScan and the API spell some libraries differently.

    'MobiDB-lite' vs 'MobiDB Lite' is the only difference observed, but the
    comparison is done on a punctuation- and case-free key so a new one does not
    silently register as a missing library.
    """
    return name.lower().replace("-", "").replace("_", "").replace(" ", "")


def api_versions(matches: Path) -> dict[str, str]:
    with gzip.open(matches, "rt", encoding="utf-8") as fh:
        raw = json.load(fh)
    versions: dict[str, str] = {}
    for batch in raw["batches"]:
        if batch.get("status") != 200:
            raise SystemExit(f"{matches}: batch {batch.get('batch')} is not 200")
        for result in batch["body"]["results"]:
            for match in result["matches"]:
                release = match["signature"]["signatureLibraryRelease"]
                library, version = release["library"], release["version"]
                seen = versions.setdefault(library, version)
                if seen != version:
                    raise SystemExit(
                        f"{matches}: {library} appears with both {seen!r} and "
                        f"{version!r}; the response is not from one release")
    return versions


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    ap.add_argument("--databases", type=Path, default=DEFAULT_DATABASES)
    args = ap.parse_args()

    api = api_versions(args.matches)
    release = json.loads(args.databases.read_text())
    by_key = {normalise(k): (k, v) for k, v in release.items()}

    print(f"InterPro release in {args.databases.name}: "
          f"{release.get('InterPro', '(not stated)')}")
    print(f"{'library (API)':<22} {'API':<12} {'InterPro release':<12} verdict")
    match = mismatch = missing = 0
    for library, version in sorted(api.items()):
        found = by_key.get(normalise(library))
        if found is None:
            print(f"{library:<22} {version:<12} {'-':<12} NOT LISTED")
            missing += 1
        elif found[1] == version:
            print(f"{library:<22} {version:<12} {found[1]:<12} match")
            match += 1
        else:
            print(f"{library:<22} {version:<12} {found[1]:<12} MISMATCH")
            mismatch += 1
    print(f"\n{match} match / {mismatch} mismatch / {missing} not listed, "
          f"out of {len(api)} libraries in the API response")
    return 0 if mismatch == missing == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
