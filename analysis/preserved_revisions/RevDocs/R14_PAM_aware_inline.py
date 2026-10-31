#!/usr/bin/env python3
import csv, hashlib, json, sys
from pathlib import Path
from collections import Counter, defaultdict

EXPECTED_SHA = "04687e2ea261caecee18065dcd2523313aca345ac865eb0483261976f703d85c"
N_TOTAL, N_MAPPED, N_SAM = 530, 515, 35985
RC = str.maketrans("ACGTN", "TGCAN")

def revcomp(s):
    return s.translate(RC)[::-1]

def audit(manifest, fasta, sam, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    assert hashlib.sha256(manifest.read_bytes()).hexdigest() == EXPECTED_SHA, "Frozen SHA mismatch"

    with manifest.open(newline="") as f:
        data = list(csv.DictReader(f))

    assert len(data) == N_TOTAL, "Unexpected frozen panel size"
    assert len(set(r["population_site_id"] for r in data)) == len(data), "Duplicate IDs"

    guides = {
        r["population_site_id"]: r
        for r in data
        if r["genomic_seqid"] in ("2L", "2R", "3L", "3R", "X")
    }
    assert len(guides) == N_MAPPED, "Unexpected mapped guide count"

    idx = {}
    with open(str(fasta) + ".fai") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            idx[p[0]] = tuple(map(int, p[1:5]))

    genome = open(fasta, "rb")

    def fetch(ch, st, en):
        if ch not in idx:
            raise ValueError("Missing reference contig " + ch)
        length, off, lb, lw = idx[ch]
        if st < 1 or en > length or en < st:
            return None
        a, b = st - 1, en - 1
        first = off + (a // lb) * lw + (a % lb)
        last = off + (b // lb) * lw + (b % lb)
        genome.seek(first)
        v = genome.read(last - first + 1)
        v = v.replace(b"\n", b"").replace(b"\r", b"").decode().upper()
        assert len(v) == en - st + 1, "FASTA index/read mismatch"
        return v

    for sid, g in guides.items():
        assert len(g["protospacer_20nt"]) == 20, sid
        assert g["pam"].upper()[1:] == "GG", sid
        assert int(g["genomic_end"]) - int(g["genomic_start"]) + 1 == 23, sid

    cols = [
        "population_site_id", "gene_id", "contig", "start1",
        "strand", "mismatches_20nt", "reference_PAM",
        "canonical_NGG", "intended_on_target"
    ]

    counts = defaultdict(Counter)
    targets = Counter()
    seen = set()
    n = 0

    try:
        with sam.open() as s, \
             (outdir / "R14_all_alignments_PAM_classified.csv").open("w", newline="") as o:

            w = csv.DictWriter(o, cols)
            w.writeheader()

            for line in s:
                if line.startswith("@"):
                    continue

                p = line.rstrip("\n").split("\t")
                if len(p) < 11:
                    raise ValueError("Malformed SAM record")

                sid = p[0]
                flag = int(p[1])

                if flag & 4:
                    continue

                assert sid in guides, "Unknown guide " + sid
                assert p[5] == "20M", "Expected ungapped 20M alignment"

                chrom = p[2]
                pos = int(p[3])
                strand = "-" if flag & 16 else "+"

                sig = (sid, chrom, pos, strand)
                assert sig not in seen, "Duplicate genomic alignment"
                seen.add(sig)

                raw = fetch(chrom, pos, pos + 19)
                assert raw is not None, "Off-reference alignment"

                q = raw if strand == "+" else revcomp(raw)

                mm = sum(
                    x != y for x, y in
                    zip(guides[sid]["protospacer_20nt"].upper(), q)
                )
                assert mm <= 3, "Mismatches exceed Bowtie threshold"

                nm = [
                    int(t[5:]) for t in p[11:]
                    if t.startswith("NM:i:")
                ]
                assert not nm or nm == [mm], "SAM NM disagrees with reference"

                upstream = fetch(chrom, pos - 3, pos - 1) if strand == "-" else None
                pam = (
                    fetch(chrom, pos + 20, pos + 22)
                    if strand == "+"
                    else (revcomp(upstream) if upstream else None)
                )

                canonical = bool(pam and pam[1:] == "GG")
                g = guides[sid]

                intended = (
                    chrom == g["genomic_seqid"]
                    and strand == g["guide_genomic_strand"]
                    and pos == int(g["genomic_start"]) + (3 if strand == "-" else 0)
                )

                if intended:
                    assert mm == 0 and pam == g["pam"].upper(), \
                        "Intended target sequence/PAM mismatch"
                    targets[sid] += 1
                else:
                    counts[sid]["other"] += 1
                    if canonical:
                        counts[sid][f"NGG_{mm}mm"] += 1
                    else:
                        counts[sid]["not_NGG"] += 1

                w.writerow(dict(
                    population_site_id=sid,
                    gene_id=g["gene_id"],
                    contig=chrom,
                    start1=pos,
                    strand=strand,
                    mismatches_20nt=mm,
                    reference_PAM=pam or "",
                    canonical_NGG=int(canonical),
                    intended_on_target=int(intended)
                ))

                n += 1
    finally:
        genome.close()

    assert n == N_SAM, f"Expected {N_SAM} alignments, found {n}"
    assert all(targets[sid] == 1 for sid in guides), \
        "Each guide must have exactly one intended match"

    table = []

    for sid, g in guides.items():
        c = counts[sid]
        row = dict(
            population_site_id=sid,
            gene_id=g["gene_id"],
            on_target=targets[sid],
            other_0to3mm=c["other"],
            other_without_NGG=c["not_NGG"],
            other_NGG_0mm=c["NGG_0mm"],
            other_NGG_1mm=c["NGG_1mm"],
            other_NGG_2mm=c["NGG_2mm"],
            other_NGG_3mm=c["NGG_3mm"]
        )
        row["any_other_NGG_0to3mm"] = int(
            sum(c[f"NGG_{k}mm"] for k in range(4)) > 0
        )
        table.append(row)

    with (outdir / "R14_per_guide_PAM_515.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, table[0].keys())
        w.writeheader()
        w.writerows(table)

    sums = {
        k: sum(r[k] for r in table)
        for k in table[0]
        if k not in ("population_site_id", "gene_id")
    }

    assert sums["other_0to3mm"] + sums["on_target"] == n
    assert (
        sums["other_without_NGG"]
        + sum(sums[f"other_NGG_{k}mm"] for k in range(4))
        == sums["other_0to3mm"]
    )

    result = dict(
        frozen_sites=N_TOTAL,
        mapped_guides=N_MAPPED,
        unresolved_UNKN=N_TOTAL - N_MAPPED,
        observed_alignments=n,
        counts=sums,
        limits=(
            "Reference-only ungapped 0-3 mismatch screen with canonical NGG; "
            "no bulges, indels, non-NGG PAMs, chromatin, or experimental cleavage inference"
        )
    )

    (outdir / "R14_PAM_audit_summary.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )

    print("PASS: Frozen manifest and all 515 intended target sites verified")
    print("Alignments:", n)
    print("Additional alignments:", sums["other_0to3mm"])
    print(
        "Additional NGG-compatible matches (0/1/2/3 mismatches):",
        *[sums[f"other_NGG_{k}mm"] for k in range(4)]
    )
    print(
        "Guides with >=1 additional NGG-compatible match:",
        sums["any_other_NGG_0to3mm"]
    )
    print("Outputs:", outdir)
    print("CAUTION: Sequence matches do not demonstrate off-target cleavage")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(
            "Usage: python3 R14_PAM_aware_inline.py "
            "manifest.csv AgamP4.fa Bowtie.sam output_directory"
        )
    audit(*[Path(x) for x in sys.argv[1:]])
