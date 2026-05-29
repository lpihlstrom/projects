##Project: GP2_EUR_YOPD_metaGWAS##

##Author: MV##

#!/usr/bin/env python3
"""
Filter a VEP-annotated VCF for:
at least N HOM-ALT cases, at most M HOM-ALT controls (default 0)
and likely damaging variant judged by:
always keep clear LoF / essential splice / IMPACT=HIGH (regardless of scores)
otherwise keep if any score threshold is met (CADD / SpliceAI / missense predictors)
"""

import argparse
import re
import sys
from typing import Dict, List, Optional, Tuple

from cyvcf2 import VCF


def read_names(path: str) -> List[str]:
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            out.append(line)
    return out


def parse_csq_format_from_header(raw_header: str) -> List[str]:
    """
    Extract CSQ field format list from VCF header line
    """
    for line in raw_header.splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            m = re.search(r"Format:\s*([^\">]+)", line)
            if not m:
                break
            fmt = m.group(1).strip()
            return fmt.split("|")
    raise RuntimeError("Could not find CSQ format in VCF header")


def safe_float_max(x: Optional[str]) -> Optional[float]:
    """
    Safely get the max number from numeric fields that may be single values or &-separated lists
    """
    if x is None:
        return None
    x = str(x).strip()
    if x in ("", ".", "NA"):
        return None
    parts = re.split(r"[&]", x)
    vals: List[float] = []
    for p in parts:
        p = p.strip()
        if p in ("", ".", "NA", "N"):
            continue
        try:
            vals.append(float(p))
        except ValueError:
            continue
    return max(vals) if vals else None


def any_pred_is(field: str, target: str) -> bool:
    """
    Check dnNSFP pred field for match in any subfield
    """
    if not field:
        return False
    for token in re.split(r"[&|,]", field):
        token = token.strip()
        if token == target:
            return True
    return False


def get_first_csq_for_alt(csq: str, csq_fields: List[str], alt: str) -> Dict[str, str]:
    """
    Return a dict of CSQ subfields for the entry corresponding to ALT allele.
    """
    entries = csq.split(",")
    if not entries:
        return {}
    allele_idx = 0
    chosen_parts: Optional[List[str]] = None

    for e in entries:
        parts = e.split("|")
        if len(parts) > allele_idx and parts[allele_idx] == alt:
            chosen_parts = parts
            break

    if chosen_parts is None:
        chosen_parts = entries[0].split("|")

    d: Dict[str, str] = {}
    for i, k in enumerate(csq_fields):
        d[k] = chosen_parts[i] if i < len(chosen_parts) else ""
    return d


def spliceai_max(csq_d: Dict[str, str]) -> Optional[float]:
    """
    Get max value from different SpliceAI fields
    """
    vals: List[float] = []
    for k, v in csq_d.items():
        if "SpliceAI_pred_DS_" in k or "SpliceAI_DS_" in k:
            fv = safe_float_max(v)
            if fv is not None:
                vals.append(fv)
    return max(vals) if vals else None


def is_lof_or_essential_splice(csq_d: Dict[str, str]) -> bool:
    """
    Keep clear LoF / essential splice regardless of scores (VEP impact/consequence)
    """
    impact = (csq_d.get("IMPACT") or "").strip()
    consequence = (csq_d.get("Consequence") or "").strip()

    loftee = (csq_d.get("LoF") or "").strip()
    if loftee == "HC":
        return True
    if impact == "HIGH":
        return True

    cons_terms = set(consequence.split("&")) if consequence else set()

    lof_terms = {
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "start_lost",
    }
    return len(cons_terms.intersection(lof_terms)) > 0


def is_damaging_by_scores(
    csq_d: Dict[str, str],
    min_cadd: float,
    min_spliceai: float,
    min_revel: float,
    min_clinpred: float,
    min_primateai: float,
) -> Tuple[bool, Dict[str, Optional[float]]]:
    """
    For non-LoF variants, pass if any of the score criteria hit
    """
    consequence = (csq_d.get("Consequence") or "").strip()

    cadd = safe_float_max(csq_d.get("CADD_PHRED"))
    spx = spliceai_max(csq_d)

    revel = safe_float_max(csq_d.get("REVEL_score"))
    clinpred = safe_float_max(csq_d.get("ClinPred_score"))
    primateai = safe_float_max(csq_d.get("PrimateAI_score"))

    metasvm = (csq_d.get("MetaSVM_pred") or "").strip()
    metalr = (csq_d.get("MetaLR_pred") or "").strip()

    scores_out = {
        "CADD_PHRED": cadd,
        "SpliceAI_MAX": spx,
        "REVEL": revel,
        "ClinPred": clinpred,
        "PrimateAI": primateai,
    }

    # Any strong signal passes
    if cadd is not None and cadd >= min_cadd:
        return True, scores_out
    if spx is not None and spx >= min_spliceai:
        return True, scores_out

    # Missense-specific predictors when applicable (dbNSFP)
    if "missense_variant" in consequence:
        if revel is not None and revel >= min_revel:
            return True, scores_out
        if clinpred is not None and clinpred >= min_clinpred:
            return True, scores_out
        if primateai is not None and primateai >= min_primateai:
            return True, scores_out
        if any_pred_is(metasvm, "D") or any_pred_is(metalr, "D"):
            return True, scores_out

    return False, scores_out


def count_hom_alt_for_alt_index(
    genotypes, sample_idxs: List[int], alt_index: int
) -> int:
    """
    Count samples homozygous for the given ALT allele index
    """
    n = 0
    for i in sample_idxs:
        a1, a2, _ph = genotypes[i]
        if a1 == alt_index and a2 == alt_index:
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="VEP-annotated vcf.gz")
    ap.add_argument(
        "--cases", required=True, help="File listing one case sample name per line"
    )
    ap.add_argument(
        "--controls",
        required=True,
        help="File listing one control sample name per line",
    )
    ap.add_argument(
        "--out-tsv", required=True, help="Output file with passing variants"
    )
    ap.add_argument("--min-hom-cases", type=int, default=1)
    ap.add_argument("--max-hom-controls", type=int, default=0)
    ap.add_argument("--min-cadd", type=float, default=20.0)
    ap.add_argument("--min-spliceai", type=float, default=0.2)
    ap.add_argument("--min-revel", type=float, default=0.5)
    ap.add_argument("--min-clinpred", type=float, default=0.5)
    ap.add_argument("--min-primateai", type=float, default=0.8)
    args = ap.parse_args()

    vcf = VCF(args.vcf)
    csq_fields = parse_csq_format_from_header(vcf.raw_header)

    cases = set(read_names(args.cases))
    controls = set(read_names(args.controls))

    sample_to_idx = {s: i for i, s in enumerate(vcf.samples)}
    case_idxs = [sample_to_idx[s] for s in cases if s in sample_to_idx]
    ctrl_idxs = [sample_to_idx[s] for s in controls if s in sample_to_idx]

    # Output variants
    sys.stdout.write(vcf.raw_header)

    with open(args.out_tsv, "w") as tsv:
        tsv.write(
            "\t".join(
                [
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "ID",
                    "SYMBOL",
                    "Consequence",
                    "IMPACT",
                    "CADD_PHRED",
                    "SpliceAI_MAX",
                    "REVEL",
                    "ClinPred",
                    "PrimateAI",
                    "HOM_CASES",
                    "HOM_CONTROLS",
                ]
            )
            + "\n"
        )

        for var in vcf:
            csq_str = var.INFO.get("CSQ")
            if not csq_str:
                continue

            genotypes = var.genotypes
            alts = var.ALT

            keep_this_record = False
            best_row = None  # choose best ALT by HOM_CASES

            for alt_i, alt in enumerate(alts, start=1):
                hom_cases = count_hom_alt_for_alt_index(genotypes, case_idxs, alt_i)
                if hom_cases < args.min_hom_cases:
                    continue

                hom_ctrls = (
                    count_hom_alt_for_alt_index(genotypes, ctrl_idxs, alt_i)
                    if ctrl_idxs
                    else 0
                )
                if hom_ctrls > args.max_hom_controls:
                    continue

                csq_d = get_first_csq_for_alt(csq_str, csq_fields, alt)

                # Keep clear LoF / essential splice / HIGH impact
                if is_lof_or_essential_splice(csq_d):
                    damaging = True
                    scores = {
                        "CADD_PHRED": safe_float_max(csq_d.get("CADD_PHRED")),
                        "SpliceAI_MAX": spliceai_max(csq_d),
                        "REVEL": safe_float_max(csq_d.get("REVEL_score")),
                        "ClinPred": safe_float_max(csq_d.get("ClinPred_score")),
                        "PrimateAI": safe_float_max(csq_d.get("PrimateAI_score")),
                    }
                else:
                    damaging, scores = is_damaging_by_scores(
                        csq_d,
                        min_cadd=args.min_cadd,
                        min_spliceai=args.min_spliceai,
                        min_revel=args.min_revel,
                        min_clinpred=args.min_clinpred,
                        min_primateai=args.min_primateai,
                    )

                if not damaging:
                    continue

                keep_this_record = True

                symbol = csq_d.get("SYMBOL", "")
                consequence = csq_d.get("Consequence", "")
                impact = csq_d.get("IMPACT", "")

                row = (
                    var.CHROM,
                    str(var.POS),
                    var.REF,
                    alt,
                    (var.ID or "."),
                    symbol,
                    consequence,
                    impact,
                    scores.get("CADD_PHRED"),
                    scores.get("SpliceAI_MAX"),
                    scores.get("REVEL"),
                    scores.get("ClinPred"),
                    scores.get("PrimateAI"),
                    hom_cases,
                    hom_ctrls,
                )

                # Prefer ALT with more HOM cases
                if best_row is None or hom_cases > best_row[-2]:
                    best_row = row

            if keep_this_record and best_row is not None:
                sys.stdout.write(str(var))

                def fmt(x):
                    if x is None:
                        return "."
                    if isinstance(x, float):
                        return f"{x:.3f}"
                    return str(x)

                tsv.write("\t".join(fmt(x) for x in best_row) + "\n")


if __name__ == "__main__":
    main()
