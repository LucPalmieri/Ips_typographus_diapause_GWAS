#!/bin/bash

###############################################################################
# vcftools VCF Filtering Script
#
# Purpose: Remove variants with excessive missing data
# 
# Reference: Danecek et al. (2011) Bioinformatics
# Tool: vcftools v0.1.16+
#
# Usage: bash vcftools_filtering.sh input.vcf output_prefix
###############################################################################

# Input & output files
INPUT_VCF="${1}"
OUTPUT_PREFIX="${2}"

if [ -z "$INPUT_VCF" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Usage: bash vcftools_filtering.sh input.vcf output_prefix"
    echo ""
    echo "Example:"
    echo "  bash vcftools_filtering.sh Ips_ipyrad.vcf filtered75_Ips_ipyrad"
    exit 1
fi

echo "================================================"
echo "vcftools VCF Filtering"
echo "================================================"
echo "Input VCF:        $INPUT_VCF"
echo "Output prefix:    $OUTPUT_PREFIX"
echo "Date:             $(date)"
echo ""

# Filter parameters
MAX_MISSING_RATE=0.25  # Allow up to 25% missing (--max-missing 0.75)

echo "Filtering parameters:"
echo "  Max missing rate: $MAX_MISSING_RATE (i.e., --max-missing 0.75)"
echo ""

# Run vcftools filter
# --max-missing 0.75  : Keep only loci with ≤25% missing data
# --recode            : Output in VCF format
# --out               : Output prefix

echo "Running vcftools..."

vcftools --vcf "$INPUT_VCF" \
    --max-missing 0.75 \
    --recode \
    --out "$OUTPUT_PREFIX"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ vcftools filtering completed successfully"
    echo ""
    echo "Output files:"
    echo "  - ${OUTPUT_PREFIX}.recode.vcf    (filtered VCF)"
    echo "  - ${OUTPUT_PREFIX}.log            (vcftools log)"
    echo ""
    
    # Print statistics
    echo "Filtering statistics:"
    if [ -f "${OUTPUT_PREFIX}.log" ]; then
        echo "  (see ${OUTPUT_PREFIX}.log for details)"
        grep "After filtering" "${OUTPUT_PREFIX}.log"
    fi
    
    echo ""
    echo "Next step: PLINK LD pruning"
    echo "  plink --vcf ${OUTPUT_PREFIX}.recode.vcf --indep-pairwise 50 10 0.2"
    
else
    echo ""
    echo "ERROR: vcftools filtering failed"
    exit 1
fi

echo ""
echo "================================================"
