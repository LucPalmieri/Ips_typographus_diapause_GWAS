import pandas as pd
import cyvcf2

# Define input and output file names
vcf_file = "PHENOTYPED.vcf"  # Change this to your VCF filename
output_csv = "GENOTYPE_PHENOTYPED_INTEGER.csv"

# Open VCF file
vcf = cyvcf2.VCF(vcf_file)

# Extract sample names
samples = vcf.samples
columns = ["CHROM", "POS", "ID", "REF", "ALT"] + samples  # Column names for output

# Process each variant in the VCF file
data = []
for variant in vcf:
    chrom, pos, var_id, ref, alt = variant.CHROM, variant.POS, variant.ID, variant.REF, variant.ALT[0]
    row = [chrom, pos, var_id, ref, alt]
    
    # Convert genotype calls to INTEGER representation (0, 2, 3)
    # Following the manuscript encoding scheme:
    # 0/0 (homozygous reference) â†’ 0
    # 0/1 or 1/0 (heterozygous) â†’ 2
    # 1/1 (homozygous alternate) â†’ 3
    # ./. (missing) â†’ NA
    for gt in variant.genotypes:
        if gt[0] == -1 or gt[1] == -1:  # Handle missing genotypes
            row.append("NA")
        elif gt[0] == 0 and gt[1] == 0:  # Homozygous reference
            row.append(0)
        elif gt[0] == 1 and gt[1] == 1:  # Homozygous alternate
            row.append(3)
        elif (gt[0] == 0 and gt[1] == 1) or (gt[0] == 1 and gt[1] == 0):  # Heterozygous
            row.append(2)
        else:
            # Catch any unexpected genotype patterns
            row.append("NA")
    
    data.append(row)

# Convert to Pandas DataFrame
df = pd.DataFrame(data, columns=columns)

# Save as CSV file
df.to_csv(output_csv, index=False)

print(f"Genotypes converted to integer encoding and saved to {output_csv}")
print("\nEncoding scheme used:")
print("  0/0 (homozygous reference) â†’ 0")
print("  0/1 or 1/0 (heterozygous)  â†’ 2")
print("  1/1 (homozygous alternate) â†’ 3")
print("  ./. (missing data)         â†’ NA")
