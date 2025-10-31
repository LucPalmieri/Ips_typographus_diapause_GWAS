# Data Dictionary

Complete specification of input and output files throughout the pipeline.

## Input Data Files (User-Provided)

### Raw Sequencing Data

| File | Format | Description | Example |
|------|--------|-------------|---------|
| `raw_reads_*.fastq` | FASTQ | Demultiplexed or raw sequencing reads | Sample_001_R1.fastq |
| `sorted_reads/*.fastq` | FASTQ | Pre-sorted reads by barcode | Sample_001.fastq |

**Location**: Input path specified in ipyrad params file

---

### Reference Genome & Annotations

| File | Format | Description | Example |
|------|--------|-------------|---------|
| `reference_genome.fasta` | FASTA | Reference genome sequence | Ips_typographus_LG16.fasta |
| `old_annotation.gff3` | GFF3 | Genome annotation (for Liftoff transfer) | reference_annotation.gff3 |
| `old_genome.fasta` | FASTA | Old genome assembly (for Liftoff mapping) | old_genome_v1.fasta |
| `agp_file.agp` | AGP | Assembly AGP file (optional, for Liftoff) | genome_scaffolds.agp |
| `pop_assign.txt` | TXT | Population assignment for individuals | See format below |

**Location**: Specified in analysis parameters

---

### Population Assignment File

Format: Two-column, tab-separated, no header

```
sample_001    LOW
sample_002    LOW
sample_003    HIGH
sample_004    HIGH
sample_005    NORTH
sample_006    NORTH
```

**Column 1**: Individual sample ID (must match VCF sample names)  
**Column 2**: Population assignment (LOW, HIGH, or NORTH)

---

### Phenotype Data

| File | Format | Description |
|------|--------|-------------|
| `phenotypes.txt` | TXT | Diapause phenotype labels for phenotyped individuals |

Format: Two-column, tab-separated, no header

```
sample_001    0
sample_002    1
sample_003    0
sample_004    1
```

**Column 1**: Sample ID (from phenotyped dataset)  
**Column 2**: Phenotype code:
- `0` = Obligate diapause
- `1` = Facultative diapause

---

## Intermediate Output Files by Stage

### Stage 1: Assembly (ipyrad)

| File | Format | Description | Purpose |
|------|--------|-------------|---------|
| `Ips_ipyrad.vcf` | VCF | Raw VCF from ipyrad | Contains all called SNPs & indels |
| `Ips_ipyrad.stats.txt` | TXT | Assembly statistics | QC summary |

**Location**: `Ips_ipyrad/Ips_ipyrad_outfiles/`

---

### Stage 2: Variant Filtering

| File | Format | Description | Purpose |
|------|--------|-------------|---------|
| `filtered75_Ips_ipyrad.recode.vcf` | VCF | VCF after missing data filter | vcftools output |
| `pruned_filtered75_Ips_ipyrad.bed` | BED | PLINK binary genotype format | Machine learning input |
| `pruned_filtered75_Ips_ipyrad.bim` | BIM | PLINK SNP information | SNP map |
| `pruned_filtered75_Ips_ipyrad.fam` | FAM | PLINK family information | Sample info |
| `pruned.prune.in` | TXT | SNP list after LD pruning | Variants to keep |
| `pruned_filtered75_Ips_ipyrad.recode.vcf` | VCF | Final pruned VCF | Downstream analyses |

**Important**: `pruned_filtered75_Ips_ipyrad.recode.vcf` is the main output used for all downstream analyses.

**Key Statistics** (from example):
- Input variants: 16,587 (from ipyrad)
- After missing data filter: ~14,000
- After LD pruning: 11,867 variants
- Samples: 574 individuals (297 phenotyped, 286 wild)
- Genotyping rate: ~79.8%

---

### Stage 3: Population Genetics

| File | Format | Description | Purpose |
|------|--------|-------------|---------|
| `VCA_results.scores` | TXT | Variance component analysis scores | GWAS covariate |
| `VCA_results.kinship` | TXT | Genomic relationship matrix | Population structure |
| `STRUCTURE_run_1.q` | TXT | STRUCTURE membership coefficients | Admixture proportions |
| `STRUCTURE_barplot.pdf` | PDF | Visualization of admixture | Publication figure |

**Format - VCA scores** (tab-separated):

```
sample_id     PC1        PC2        PC3
sample_001    0.1234     -0.0456    0.0078
sample_002    0.1245     -0.0423    0.0089
```

**Format - STRUCTURE output** (tab-separated):

```
sample_id    K1_LOW    K2_HIGH   K3_NORTH
sample_001   0.85      0.10      0.05
sample_002   0.80      0.15      0.05
```

---

### Stage 4: Genome Annotation

| File | Format | Description | Purpose |
|------|--------|-------------|---------|
| `new_annotation.gff3` | GFF3 | Transferred genome annotation | Feature locations |
| `new_annotation_features.txt` | TXT | Extracted genomic features | SNP mapping input |
| `SNP_feature_overlaps_GWAS_moderateBF.tsv` | TSV | SNP-to-gene mapping results | Functional annotation |

**Format - Feature file** (tab-separated):

```
chrom     start      end        feature      gene_id
LG1       1000       1500       CDS          gene_001
LG1       1500       1800       exon         gene_001
LG1       2000       2100       5_UTR        gene_001
```

**Format - SNP-feature overlaps**:

```
seqnames    ranges        MRK         gene_id       feature      gene_id.1
LG1         1234-1234     snp_001     gene_001      exon         gene_001
LG2         5678-5678     snp_002     gene_002      CDS          gene_002
```

---

### Stage 5: GWAS (BayPass)

| File | Format | Description | Purpose |
|------|--------|-------------|---------|
| `phenotyped.geno` | BayPass format | Allele count format for BayPass | GWAS input |
| `phenotyped_core_run1_summary_pi_xtx.txt` | TXT | XtX matrix from core model | Population covariance |
| `phenotyped_importance_sampling_cov_run1_summary_betai_fst.txt` | TXT | Bayes Factors with covariate | Association test results |
| `GWAS_moderateBF.txt` | TXT | SNPs with BF > 6 | Candidate loci |
| `gwas_summary.txt` | TXT | GWAS statistics summary | QC metrics |

**Format - BayPass allele counts** (BayPass-specific binary format)

**Format - Bayes Factor output** (tab-separated):

```
SNP_ID          CHROM    POS       REF   ALT    AF_OBL    AF_FAC    BF_median
snp_000001      LG1      1234      A     T      0.45      0.55      3.2
snp_000002      LG1      5678      G     C      0.30      0.70      8.9
snp_000003      LG2      9012      T     G      0.48      0.52      1.5
```

**Format - Moderate BF SNPs** (tab-separated):

```
CHROM    POS       SNP_ID          BF_median    AF_obligate    AF_facultative    gene_id
LG1      1234      snp_000001      8.9          0.30           0.70              gene_001
LG2      5678      snp_000002      7.2          0.40           0.65              gene_002
```

**Key Output**: `GWAS_moderateBF.txt` contains SNPs with BF > 6 (moderate to strong associations)

---

### Stage 6: Machine Learning

| File | Format | Description | Purpose |
|------|--------|-------------|---------|
| `GENOTYPE_PHENOTYPED_INTEGER.csv` | CSV | Integer-encoded genotypes | ML model input |
| `optuna_results/best_params.json` | JSON | Optimal hyperparameters | Model configuration |
| `SNP_importance_ranking.csv` | CSV | Ranked SNP feature importance | Candidate SNP ranking |
| `cross_validation_metrics.txt` | TXT | k-fold CV performance | Model quality assessment |
| `model_predictions.csv` | CSV | Predicted phenotype probabilities | Prediction scores |

**Format - Integer-encoded genotypes** (CSV):

```
sample_id,snp_001,snp_002,snp_003,...
sample_001,0,2,3
sample_002,2,3,NA
sample_003,3,0,2
```

**Genotype encoding**:
- `0` = Homozygous reference (0/0)
- `2` = Heterozygous (0/1 or 1/0)
- `3` = Homozygous alternate (1/1)
- `NA` = Missing (./.)

**Format - SNP importance ranking** (CSV):

```
SNP_ID,CHROM,POS,importance_score,rank,p_value
snp_001,LG1,1234,0.0456,1,0.0001
snp_002,LG2,5678,0.0423,2,0.0002
snp_003,LG1,9012,0.0389,3,0.0003
```

**Format - Cross-validation metrics** (TXT):

```
Metric              Mean ± SD
AUROC              0.89 ± 0.03
Accuracy           0.82 ± 0.04
F1 Score           0.80 ± 0.05
Sensitivity        0.85 ± 0.06
Specificity        0.80 ± 0.05
Precision          0.82 ± 0.04
```

---

### Stage 7: Functional Annotation

| File | Format | Description | Purpose |
|------|--------|-------------|---------|
| `foreground_mod_BF6.txt` | TXT | Foreground genes (SNP hits, BF>6) | GO enrichment input |
| `gene_length.txt` | TXT | Background gene lengths | Bias correction for GO |
| `Gene2GO.txt` | TXT | Gene-to-GO term mapping | GO enrichment mapping |
| `GOseq_all_terms_modBF6.tsv` | TSV | All GO terms tested | Complete GO results |
| `GOseq_significant_terms_modBF6.tsv` | TSV | Significant terms (p<0.05) | Enriched pathways |

**Format - Foreground genes**:

```
gene_id    isDE
gene_001   1
gene_002   1
gene_003   0
```

**Format - Gene lengths**:

```
gene_id    length
gene_001   1500
gene_002   2300
gene_003   1800
```

**Format - Gene2GO mapping**:

```
gene_id    GO_term
gene_001   GO:0006629
gene_001   GO:0008150
gene_002   GO:0003677
```

**Format - GO enrichment results** (TSV):

```
GO_ID           ontology    term                          p_value    FDR
GO:0006629      BP          lipid metabolism              0.001      0.025
GO:0008150      BP          biological_process           0.002      0.030
GO:0003677      MF          DNA binding                   0.005      0.045
```

**Key Output**: `GOseq_significant_terms_modBF6.tsv` contains significantly enriched GO categories

---

### Stage 8: Wild Population Predictions

| File | Format | Description | Purpose |
|------|--------|-------------|---------|
| `GENOTYPE_WILD_INTEGER.csv` | CSV | Integer-encoded wild genotypes | Prediction input |
| `VCA_wild.scores` | TXT | VCA scores for wild individuals | Covariate input |
| `wild_predictions.csv` | CSV | Predicted phenotypes for wild individuals | Results |
| `population_frequencies.txt` | TXT | Population-level phenotype frequencies | Summary stats |
| `outbreak_risk_assessment.txt` | TXT | Outbreak risk classification | Risk metrics |

**Format - Wild predictions** (CSV):

```
sample_id,population,predicted_phenotype,probability_facultative
wild_001,LOW,1,0.78
wild_002,LOW,0,0.32
wild_003,HIGH,1,0.85
wild_004,NORTH,0,0.15
```

**Predicted phenotype**: 0 = Obligate, 1 = Facultative

**Format - Population frequencies** (TXT):

```
Population    N_samples    Obligate_%    Facultative_%
LOW           96           45            55
HIGH          96           38            62
NORTH         94           72            28
```

**Format - Outbreak risk** (TXT):

```
Population    Risk_Category    Reasoning
LOW           MODERATE         55% facultative predicted
HIGH          HIGH             62% facultative predicted
NORTH         LOW              72% obligate predicted
```

---

## Output Directory Structure

```
results/
├── 01_assembly/
│   ├── Ips_ipyrad.vcf
│   └── stats/
├── 02_filtering/
│   ├── filtered75_Ips_ipyrad.recode.vcf
│   ├── pruned_filtered75_Ips_ipyrad.*
│   └── plink_report.log
├── 03_population_genetics/
│   ├── VCA_results.scores
│   ├── STRUCTURE_run_*.q
│   └── plots/
│       └── STRUCTURE_barplot.pdf
├── 04_annotation/
│   ├── new_annotation.gff3
│   ├── new_annotation_features.txt
│   └── SNP_feature_overlaps_GWAS_moderateBF.tsv
├── 05_gwas/
│   ├── phenotyped.geno
│   ├── baypass_results/
│   │   ├── *_core_run*.txt
│   │   └── *_importance_sampling*.txt
│   ├── GWAS_moderateBF.txt
│   └── gwas_summary.txt
├── 06_machine_learning/
│   ├── optuna_results/
│   │   ├── best_params.json
│   │   └── optimization_plots.html
│   ├── GENOTYPE_PHENOTYPED_INTEGER.csv
│   ├── SNP_importance_ranking.csv
│   ├── cross_validation_metrics.txt
│   └── model_predictions.csv
├── 07_functional_annotation/
│   ├── GO_foreground_genes.txt
│   ├── GOseq_all_terms_modBF6.tsv
│   ├── GOseq_significant_terms_modBF6.tsv
│   └── enrichment_plots/
├── 08_visualization/
│   ├── manhattan_plot.pdf
│   ├── covariance_heatmap.pdf
│   └── structure_barplot.pdf
└── 09_wild_predictions/
    ├── wild_predictions.csv
    ├── population_frequencies.txt
    └── outbreak_risk_assessment.txt
```

---

## File Size Estimates

| Type | Typical Size |
|------|--------------|
| Raw VCF (11,867 SNPs × 574 samples) | ~50-80 MB |
| Pruned VCF | ~30-50 MB |
| Encoded genotypes (CSV) | ~10-20 MB |
| BayPass output (all runs) | ~100-200 MB |
| GO enrichment results | <1 MB |
| Final predictions | <1 MB |
| **Total output** | ~200-500 MB |

---

## Important Notes

1. **VCF format**: Uses `./. ` for missing genotypes (standard VCF)
2. **SNP IDs**: Generated as `CHROM_POS_REF_ALT`
3. **Allele frequency**: Calculated from allele counts (may include missing data)
4. **P-values**: When multiple testing, apply FDR correction (see goseq output)
5. **Phenotype coding**: Consistent throughout (0=obligate, 1=facultative)

---

For file-specific questions, see the corresponding stage README in `scripts/[stage]/README.md`
