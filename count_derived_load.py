#!/usr/bin/env python3

import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path


# Change this if your ancestral allele INFO field has another name.
AA_TAG = "AA"

# Priority prevents one SNP being counted in several load categories.
CATEGORY_PRIORITY = {
    "LoF_HIGH": 3,
    "missense_MODERATE": 2,
    "synonymous_LOW": 1,
}


def open_text(filename):
    """Open plain-text or gzip-compressed files."""
    if str(filename).endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, "rt")


def parse_info(info_string):
    """Convert the VCF INFO column into a dictionary."""
    info = {}

    if info_string in {".", ""}:
        return info

    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True

    return info


def parse_csq_header(line):
    """
    Extract VEP CSQ field names from:
    Format: Allele|Consequence|IMPACT|SYMBOL|...
    """
    match = re.search(r"Format: ([^\">]+)", line)

    if not match:
        raise ValueError("Could not identify CSQ field structure.")

    return match.group(1).strip().split("|")


def classify_variant(csq_value, csq_fields):
    """
    Classify one SNP once.

    Priority:
      1. Any HIGH-impact transcript annotation
      2. missense_variant
      3. synonymous_variant
    """
    field_index = {name: index for index, name in enumerate(csq_fields)}

    required = {"Consequence", "IMPACT"}
    missing = required.difference(field_index)

    if missing:
        raise ValueError(
            f"Required CSQ fields missing: {', '.join(sorted(missing))}"
        )

    observed_categories = set()

    for annotation in csq_value.split(","):
        values = annotation.split("|")

        # Pad incomplete annotation records.
        if len(values) < len(csq_fields):
            values.extend([""] * (len(csq_fields) - len(values)))

        consequence = values[field_index["Consequence"]]
        impact = values[field_index["IMPACT"]]

        consequence_terms = set(consequence.split("&"))

        if impact == "HIGH":
            observed_categories.add("LoF_HIGH")

        if (
            impact == "MODERATE"
            and "missense_variant" in consequence_terms
        ):
            observed_categories.add("missense_MODERATE")

        if (
            impact == "LOW"
            and "synonymous_variant" in consequence_terms
        ):
            observed_categories.add("synonymous_LOW")

    if not observed_categories:
        return None

    return max(
        observed_categories,
        key=lambda category: CATEGORY_PRIORITY[category]
    )


def derived_genotype_count(gt, ref, alt, ancestral):
    """
    Return:
      derived allele copies,
      called allele copies,
      homozygous-derived status.

    Supports phased and unphased diploid genotypes.
    """
    if not gt or gt in {".", "./.", ".|."}:
        return None

    alleles = re.split(r"[\/|]", gt)

    # Ignore partially missing genotypes.
    if not alleles or any(a == "." for a in alleles):
        return None

    try:
        allele_numbers = [int(a) for a in alleles]
    except ValueError:
        return None

    # This analysis is restricted to biallelic SNPs.
    if any(a not in {0, 1} for a in allele_numbers):
        return None

    called_alleles = len(allele_numbers)
    alt_count = sum(a == 1 for a in allele_numbers)

    ancestral = ancestral.upper()
    ref = ref.upper()
    alt = alt.upper()

    if ancestral == ref:
        # ALT is derived.
        derived_count = alt_count
    elif ancestral == alt:
        # REF is derived.
        derived_count = called_alleles - alt_count
    else:
        # AA does not match REF or ALT.
        return None

    homozygous_derived = int(derived_count == called_alleles)

    return derived_count, called_alleles, homozygous_derived


def main():
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: python count_derived_load.py "
            "caracal.vep.vcf caracal.derived_load.tsv"
        )

    input_vcf = Path(sys.argv[1])
    output_tsv = Path(sys.argv[2])

    samples = []
    csq_fields = None

    # counts[sample][category][metric]
    counts = defaultdict(
        lambda: defaultdict(
            lambda: {
                "derived_alleles": 0,
                "homozygous_derived_sites": 0,
                "heterozygous_sites": 0,
                "called_alleles": 0,
                "called_sites": 0,
            }
        )
    )

    site_counts = defaultdict(int)

    skipped_no_csq = 0
    skipped_no_aa = 0
    skipped_bad_aa = 0

    with open_text(input_vcf) as handle:
        for line in handle:

            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = parse_csq_header(line)
                continue

            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                samples = columns[9:]
                continue

            if line.startswith("#"):
                continue

            if csq_fields is None:
                raise RuntimeError("No VEP CSQ header found in the VCF.")

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            chrom, pos, variant_id, ref, alt = fields[:5]
            info_string = fields[7]
            format_string = fields[8]
            sample_values = fields[9:]

            # Skip multiallelic sites.
            if "," in alt:
                continue

            info = parse_info(info_string)

            if "CSQ" not in info:
                skipped_no_csq += 1
                continue

            if AA_TAG not in info:
                skipped_no_aa += 1
                continue

            ancestral = str(info[AA_TAG]).split("|")[0].upper()

            # Remove common ancestral-allele suffixes such as A, C, G or T.
            ancestral = ancestral.split(",")[0]

            if ancestral not in {ref.upper(), alt.upper()}:
                skipped_bad_aa += 1
                continue

            category = classify_variant(info["CSQ"], csq_fields)

            if category is None:
                continue

            site_counts[category] += 1

            format_fields = format_string.split(":")

            if "GT" not in format_fields:
                raise RuntimeError(
                    f"No GT field at {chrom}:{pos}"
                )

            gt_index = format_fields.index("GT")

            for sample, sample_field in zip(samples, sample_values):
                values = sample_field.split(":")

                if gt_index >= len(values):
                    continue

                gt = values[gt_index]

                result = derived_genotype_count(
                    gt=gt,
                    ref=ref,
                    alt=alt,
                    ancestral=ancestral,
                )

                if result is None:
                    continue

                derived, called_alleles, hom_derived = result

                record = counts[sample][category]

                record["derived_alleles"] += derived
                record["called_alleles"] += called_alleles
                record["called_sites"] += 1
                record["homozygous_derived_sites"] += hom_derived

                if derived == 1 and called_alleles == 2:
                    record["heterozygous_sites"] += 1

    categories = [
        "LoF_HIGH",
        "missense_MODERATE",
        "synonymous_LOW",
    ]

    with open(output_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")

        writer.writerow([
            "sample",
            "category",
            "available_variant_sites",
            "called_sites",
            "called_alleles",
            "derived_alleles",
            "homozygous_derived_sites",
            "heterozygous_sites",
            "derived_alleles_per_1000_called_alleles",
        ])

        for sample in samples:
            for category in categories:
                record = counts[sample][category]
                called = record["called_alleles"]

                if called > 0:
                    normalized = (
                        1000 * record["derived_alleles"] / called
                    )
                    normalized_text = f"{normalized:.6f}"
                else:
                    normalized_text = "NA"

                writer.writerow([
                    sample,
                    category,
                    site_counts[category],
                    record["called_sites"],
                    called,
                    record["derived_alleles"],
                    record["homozygous_derived_sites"],
                    record["heterozygous_sites"],
                    normalized_text,
                ])

    print(f"Output written to: {output_tsv}")
    print("\nUnique classified sites:")

    for category in categories:
        print(f"  {category}: {site_counts[category]}")

    print("\nSkipped variants:")
    print(f"  No CSQ annotation: {skipped_no_csq}")
    print(f"  No {AA_TAG} annotation: {skipped_no_aa}")
    print(f"  Ancestral allele did not match REF/ALT: {skipped_bad_aa}")


if __name__ == "__main__":
    main()
