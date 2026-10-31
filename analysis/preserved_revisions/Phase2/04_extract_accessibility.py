#!/usr/bin/env python3
"""Exact 23-bp Ag1000G Phase 2 accessibility for the FROZEN MosqEdit-R panel.

Reads only slices of /<arm>/is_accessible and /<arm>/pos from the original
Phase 2 accessibility.h5. Does not download VCFs or change candidate sites.

Outputs 530-site rows (515 evaluated + 15 UNKN explicitly unresolved) and a
base-level table; suitable for downstream reproducible R analysis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

try:
    import h5py
except ImportError as exc:
    raise SystemExit("Install HDF5 Python bindings: sudo apt install -y python3-h5py") from exc

ARMS = ("2L", "2R", "3L", "3R", "X")
EXPECTED_LEN = {"2L": 49364325, "2R": 61545105, "3L": 41963435,
                "3R": 53200684, "X": 24393108}
EXPECTED_MANIFEST_SHA = "04687e2ea261caecee18065dcd2523313aca345ac865eb0483261976f703d85c"
EXPECTED_MASK_BYTES = 2005427410

SITE_FIELDS = [
    "population_site_id", "gene_id", "genomic_seqid", "genomic_start",
    "genomic_end", "guide_genomic_strand", "n_bp", "n_accessible",
    "fraction_accessible", "fully_accessible_23bp", "n_pam_bp",
    "n_pam_accessible", "pam_fully_accessible", "n_protospacer_bp",
    "n_protospacer_accessible", "protospacer_fully_accessible", "status",
]
BASE_FIELDS = ["population_site_id", "gene_id", "chromosome", "genomic_pos_1based",
               "position_in_genomic_23mer_1based", "region", "is_accessible"]


def _read_manifest(path, sha, n_expected):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if sha and digest != sha:
        raise ValueError(f"FROZEN MANIFEST SHA256 MISMATCH: {digest} != {sha}; STOP, do not proceed")
    with path.open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != n_expected:
        raise ValueError(f"Expected {n_expected} frozen rows, found {len(rows)}")
    ids = [r["population_site_id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate population_site_id in manifest")
    for r in rows:
        s, e = int(r["genomic_start"]), int(r["genomic_end"])
        if s < 1 or e - s + 1 != 23:
            raise ValueError(f"Site coordinates not 23 bp: {r['population_site_id']}")
        if r["genomic_seqid"] not in (*ARMS, "UNKN"):
            raise ValueError(f"Unknown chromosome: {r['population_site_id']}")
        if r["guide_genomic_strand"] not in ("+", "-"):
            raise ValueError(f"Invalid strand: {r['population_site_id']}")
        if len(r["protospacer_20nt"]) != 20 or len(r["pam"]) != 3:
            raise ValueError(f"Invalid guide or PAM length: {r['population_site_id']}")
    return rows, digest


def run(manifest, h5path, outdir, expected_sha=EXPECTED_MANIFEST_SHA,
        expected_count=530, expected_lengths=None, expected_bytes=EXPECTED_MASK_BYTES):
    expected_lengths = EXPECTED_LEN if expected_lengths is None else expected_lengths
    rows, digest = _read_manifest(manifest, expected_sha, expected_count)
    if not h5path.is_file():
        raise ValueError(f"No accessibility file at: {h5path}")
    if expected_bytes is not None and h5path.stat().st_size != expected_bytes:
        raise ValueError(f"Unexpected HDF5 file size: {h5path.stat().st_size}, expected {expected_bytes}")
    outdir.mkdir(parents=True, exist_ok=True)
    site_file = outdir / "AR1_frozen_site_accessibility_530.csv"
    base_file = outdir / "AR1_frozen_base_accessibility_515sites.csv"
    summary_file = outdir / "AR1_accessibility_extraction_summary.json"
    evaluated = Counter()
    complete = Counter()
    results = []
    bases = []
    with h5py.File(h5path, "r") as f:
        for arm in ARMS:
            if arm not in f:
                raise ValueError(f"Missing chromosome group /{arm}")
            g = f[arm]
            if "is_accessible" not in g or "pos" not in g:
                raise ValueError(f"Missing /{arm}/is_accessible or /{arm}/pos")
            mask, pos = g["is_accessible"], g["pos"]
            length = int(expected_lengths[arm])
            if len(mask.shape) != 1 or len(pos.shape) != 1 or len(mask) != length or len(pos) != length:
                raise ValueError(f"HDF5 chromosome length inconsistent for {arm}: "
                                 f"mask={mask.shape}, pos={pos.shape}, expected={length}")
            # Check position convention before interpreting genomic coordinates.
            for i in (0, length // 2, length - 1):
                observed = int(pos[i])
                if observed != i + 1:
                    raise ValueError(f"HDF5 /{arm}/pos[{i}]={observed}, expected {i+1}; "
                                     "STOP rather than assume an off-by-one correction")
            print(f"{arm}: length={length:,}, mask_dtype={mask.dtype}, pos_1based=confirmed", flush=True)
        for r in rows:
            chrom = r["genomic_seqid"]
            s = int(r["genomic_start"])
            e = int(r["genomic_end"])
            strand = r["guide_genomic_strand"]
            common = {"population_site_id": r["population_site_id"],
                      "gene_id": r["gene_id"], "genomic_seqid": chrom,
                      "genomic_start": s, "genomic_end": e,
                      "guide_genomic_strand": strand}
            if chrom == "UNKN":
                results.append({**common, "n_bp": 23, "n_accessible": "",
                                "fraction_accessible": "", "fully_accessible_23bp": "",
                                "n_pam_bp": 3, "n_pam_accessible": "",
                                "pam_fully_accessible": "", "n_protospacer_bp": 20,
                                "n_protospacer_accessible": "",
                                "protospacer_fully_accessible": "",
                                "status": "UNKN_COORDINATES_NOT_EVALUABLE"})
                continue
            mask = f[chrom]["is_accessible"]
            pos = f[chrom]["pos"]
            if e > len(mask):
                raise ValueError(f"Out-of-range site {r['population_site_id']} on {chrom}: {s}-{e}")
            if int(pos[s-1]) != s or int(pos[e-1]) != e:
                raise ValueError(f"Position check failed: {r['population_site_id']}")
            values = [int(x) for x in mask[s-1:e]]
            if len(values) != 23 or any(x not in (0, 1) for x in values):
                raise ValueError(f"Nonbinary/unexpected accessibility mask for {r['population_site_id']}")
            # All manifest intervals are [protospacer(20nt), PAM(3nt)] in guide
            # orientation. On reverse genomic strand the PAM occupies the LEFT
            # (lower-coordinate) 3 bp, unlike the forward genomic strand.
            pam_idx = set(range(20, 23)) if strand == "+" else set(range(0, 3))
            proto_idx = set(range(23)) - pam_idx
            n_acc = sum(values)
            n_pam = sum(values[i] for i in pam_idx)
            n_proto = sum(values[i] for i in proto_idx)
            evaluated[chrom] += 1
            complete[chrom] += int(n_acc == 23)
            results.append({**common, "n_bp": 23, "n_accessible": n_acc,
                            "fraction_accessible": f"{n_acc/23:.8f}",
                            "fully_accessible_23bp": int(n_acc == 23),
                            "n_pam_bp": 3, "n_pam_accessible": n_pam,
                            "pam_fully_accessible": int(n_pam == 3),
                            "n_protospacer_bp": 20, "n_protospacer_accessible": n_proto,
                            "protospacer_fully_accessible": int(n_proto == 20),
                            "status": "EVALUATED_AG1000G_PHASE2_SITE_MASK"})
            for i, v in enumerate(values):
                bases.append({"population_site_id": r["population_site_id"],
                              "gene_id": r["gene_id"], "chromosome": chrom,
                              "genomic_pos_1based": s+i,
                              "position_in_genomic_23mer_1based": i+1,
                              "region": "PAM" if i in pam_idx else "PROTOSPACER",
                              "is_accessible": v})
    if sum(evaluated.values()) != 515 or sum(r["genomic_seqid"] == "UNKN" for r in rows) != 15:
        raise ValueError("Unexpected chromosome/UNKN breakdown; no outputs written")
    for path, fields, data in ((site_file, SITE_FIELDS, results),
                               (base_file, BASE_FIELDS, bases)):
        with path.open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(data)
    summary = {"source": "Ag1000G Phase 2 AR1 accessibility.h5",
               "mask": "is_accessible; aggregate genomic site-level mask, not sample-specific genotype QC",
               "manifest_sha256": digest,
               "manifest_n_sites": len(rows), "n_evaluated": sum(evaluated.values()),
               "n_unresolved_UNKN": 15,
               "n_bp_evaluated": len(bases),
               "per_arm_n_sites": dict(evaluated),
               "per_arm_n_fully_accessible_sites": dict(complete),
               "n_fully_accessible_sites": sum(complete.values()),
               "coordinate_convention": "manifest and HDF5 pos 1-based inclusive; Python array indices 0-based",
               "pam_mapping": "forward: upper three positions; reverse: lower three positions",
               "important_limitations": [
                   "Accessibility does not establish a reference genotype in every individual",
                   "No Phase 2 SNP record is not equivalent to a verified invariant guide",
                   "Individual genotype missingness, allele frequencies, PAM/protospacer variants and haplotypes require separate evaluation",
                   "15 UNKN sites must not be classified as conserved or inaccessible from this file",
               ]}
    summary_file.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Processed {sum(evaluated.values())} target sites and {len(bases)} target bases", flush=True)
    print(f"UNKN sites kept as explicitly unevaluable: 15", flush=True)
    print(f"Fully accessible 23-bp targets: {sum(complete.values())}/{sum(evaluated.values())}", flush=True)
    print(f"Outputs:\n  {site_file}\n  {base_file}\n  {summary_file}", flush=True)
    return summary


def main():
    root = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, default=root / "13_population_site_manifest_FROZEN.csv")
    p.add_argument("--accessibility", type=Path, default=root.parent / "downloads" / "accessibility.h5")
    p.add_argument("--out", type=Path, default=root / "output")
    args = p.parse_args()
    try:
        run(args.manifest, args.accessibility, args.out)
    except (ValueError, OSError, KeyError) as e:
        print(f"STOP: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
