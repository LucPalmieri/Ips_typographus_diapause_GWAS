#!/bin/bash

###############################################################################
# Liftoff Annotation Transfer Script
# 
# Purpose: Transfer genome annotations from old assembly to new assembly
# 
# Reference: Shumate & Salzberg (2021) Bioinformatics
# https://github.com/agshumate/Liftoff
#
# Usage: bash liftoff_annotation.sh
###############################################################################

# Set input/output paths (modify as needed)
OLD_ANNOTATION="old_annotation.gff3"      # Original genome annotation
OLD_GENOME="old_genome.fasta"              # Original reference genome
NEW_GENOME="new_genome.fasta"              # New reference genome
AGP_FILE="agp_file.agp"                    # AGP mapping file (optional)
OUTPUT_ANNOTATION="new_annotation.gff3"    # Output annotation file

###############################################################################
# Liftoff Annotation Transfer
###############################################################################

echo "================================================"
echo "Liftoff Annotation Transfer"
echo "================================================"
echo "Date: $(date)"
echo ""
echo "Input annotation:    $OLD_ANNOTATION"
echo "Old genome:          $OLD_GENOME"
echo "New genome:          $NEW_GENOME"
echo "AGP file:            $AGP_FILE"
echo "Output annotation:   $OUTPUT_ANNOTATION"
echo ""

# Check if input files exist
if [ ! -f "$OLD_ANNOTATION" ]; then
    echo "ERROR: Old annotation file not found: $OLD_ANNOTATION"
    exit 1
fi

if [ ! -f "$OLD_GENOME" ]; then
    echo "ERROR: Old genome file not found: $OLD_GENOME"
    exit 1
fi

if [ ! -f "$NEW_GENOME" ]; then
    echo "ERROR: New genome file not found: $NEW_GENOME"
    exit 1
fi

# Run Liftoff
# Parameters:
#   -g: Input annotation (GFF3 format)
#   -o: Output annotation file
#   -a: AGP file for scaffolding (optional, improves accuracy)
#   -chrOrder: Chromosome order file (optional)

echo "Running Liftoff..."
echo ""

if [ -f "$AGP_FILE" ]; then
    echo "Using AGP file for additional mapping information..."
    liftoff -g "$OLD_ANNOTATION" \
        -o "$OUTPUT_ANNOTATION" \
        -a "$AGP_FILE" \
        "$NEW_GENOME" \
        "$OLD_GENOME"
else
    echo "No AGP file provided, proceeding without additional mapping..."
    liftoff -g "$OLD_ANNOTATION" \
        -o "$OUTPUT_ANNOTATION" \
        "$NEW_GENOME" \
        "$OLD_GENOME"
fi

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "================================================"
    echo "Liftoff completed successfully!"
    echo "================================================"
    echo ""
    echo "Output files created:"
    echo "  - $OUTPUT_ANNOTATION"
    
    # Count features
    echo ""
    echo "Annotation statistics:"
    TOTAL_FEATURES=$(grep -v '^#' "$OUTPUT_ANNOTATION" | wc -l)
    echo "  Total features: $TOTAL_FEATURES"
    
    GENES=$(grep -w 'gene' "$OUTPUT_ANNOTATION" | wc -l)
    echo "  Genes: $GENES"
    
    MRNAS=$(grep -w 'mRNA' "$OUTPUT_ANNOTATION" | wc -l)
    echo "  mRNAs: $MRNAS"
    
    EXONS=$(grep -w 'exon' "$OUTPUT_ANNOTATION" | wc -l)
    echo "  Exons: $EXONS"
    
    # Additional output files that Liftoff creates
    if [ -f "${OUTPUT_ANNOTATION%.gff3}.unlifted.features" ]; then
        UNLIFTED=$(wc -l < "${OUTPUT_ANNOTATION%.gff3}.unlifted.features")
        echo "  Features NOT lifted: $UNLIFTED"
    fi
    
    echo ""
    echo "Next step: Use $OUTPUT_ANNOTATION for SNP-to-gene mapping"
    echo "           (Gene_feature_sorting.rmd)"
    
else
    echo ""
    echo "ERROR: Liftoff failed with exit code $?"
    exit 1
fi

###############################################################################
# Optional: Generate feature-only GFF3 for SNP mapping
###############################################################################

# Extract only CDS, exon, and UTR features for efficient SNP overlap detection
OUTPUT_FEATURES="${OUTPUT_ANNOTATION%.gff3}_features.txt"

echo ""
echo "Extracting genomic features for SNP mapping..."
echo ""

# Create tab-separated file with: chrom, start, end, feature, gene_id
grep -v '^#' "$OUTPUT_ANNOTATION" | \
    awk -F'\t' '$3 ~ /^(CDS|exon|five_prime_UTR|three_prime_UTR)$/ {
        chrom = $1
        start = $4
        end = $5
        feature = $3
        
        # Extract gene_id and parent info from attributes (column 9)
        attrs = $9
        if (match(attrs, /gene_id=([^;]+)/)) {
            gene_id = substr(attrs, RSTART+8, RLENGTH-9)
        } else if (match(attrs, /Parent=([^;]+)/)) {
            gene_id = substr(attrs, RSTART+7, RLENGTH-8)
        } else {
            gene_id = "NA"
        }
        
        print chrom "\t" start "\t" end "\t" feature "\t" gene_id
    }' | sort -k1,1 -k2,2n > "$OUTPUT_FEATURES"

if [ $? -eq 0 ] && [ -s "$OUTPUT_FEATURES" ]; then
    echo "✓ Features extracted: $OUTPUT_FEATURES"
    echo "  This file is ready for SNP-to-gene mapping"
    echo "  (Use as input for Gene_feature_sorting.rmd)"
    echo ""
    echo "  File format:"
    echo "  chrom    start       end       feature      gene_id"
    head -5 "$OUTPUT_FEATURES"
else
    echo "⚠ Warning: Feature extraction produced empty file or failed"
fi

echo ""
echo "================================================"
echo "Liftoff annotation transfer complete!"
echo "================================================"
