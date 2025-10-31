# Set the working directory
setwd("C:/Users/LPalmieriRocha/Documents/RADips_results/Baypass/phenotyped")

```{r}
#source the baypass R functions (check PATH)
source("C:/Users/LPalmieriRocha/Documents/RADips_results/Baypass/baypass_utils.R")
#upload the estimated Omega matrix
omega=as.matrix(read.table("noWild_core_mat_omega.out"))
pop.names=c("Ndia","Hdia","Hnod","Ldia","Lnod")
dimnames(omega)=list(pop.names,pop.names)
# Visualization of the matrix
# Using SVD decomposition
plot.omega(omega=omega,pop.names=pop.names)

#heatmap and hierarchical clustering tree (using the average agglomeration method)
cor.mat=cov2cor(omega)
hclust.ave <- function(x) hclust(x, method="average")
heatmap(1-cor.mat,hclustfun = hclust.ave,
        main=expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

--------------------------------------------------------------------------------
  # Load necessary libraries
  library(gplots)

# Close all previous graphic devices
while (!is.null(dev.list())) dev.off()

# Open a new graphics device (use x11 for Windows/Linux, quartz for macOS)
png("heatmap_output.png", width = 800, height = 800)

# Define the omega matrix
omega <- matrix(c(0.1167109, -0.00281574, -0.02369615, -0.01686128, -0.02605021,
                  -0.00281574, 0.08377436, 0.04001008, -0.02469262, -0.032499,
                  -0.02369615, 0.04001008, 0.03726688, -0.01275034, -0.0179783,
                  -0.01686128, -0.02469262, -0.01275034, 0.07345348, 0.05334868,
                  -0.02605021, -0.032499, -0.0179783, 0.05334868, 0.05781831),
                nrow = 5, byrow = TRUE)

# Define population names
pop.names <- c("Ldia", "Lnod", "Ndia", "Hdia", "Hnod")
rownames(omega) <- pop.names
colnames(omega) <- pop.names

# Calculate the correlation matrix
cor.mat <- cor(omega)

# Create the heatmap with hierarchical clustering and annotations
heatmap.2(1 - cor.mat,
          dendrogram = "both", 
          Rowv = TRUE, 
          Colv = TRUE, 
          trace = "none", 
          col = colorRampPalette(c("darkgreen", "darkorange"))(50), 
          key = TRUE, 
          key.title = "1 - Correlation",
          density.info = "none", 
          scale = "none", 
          labRow = pop.names, 
          labCol = pop.names,
          margins = c(7,7),
          cellnote = round((1 - cor.mat), 3), # Convert to percentage and round #round(1 - cor.mat, 2), # Annotate with rounded values
          notecol = "white", 
          notecex = 2, 
          hclustfun = function(x) hclust(x, method = "average"))

# Add a title
title(main = expression("Heatmap of" ~ hat(Omega) ~ "(" * d[ij] * "= 1 - " * rho[ij] * ")"))

# Close the graphics device if needed
dev.off()

--------------------------------------------------------------------------------
#Estimates of the XtX differentiation measures (using the calibrated XtXst estimator)
noWild.snp=read.table("noWild_core_summary_pi_xtx.out",h=T)
#check behavior of the p-values associated to the XtXst estimator
hist(10**(-1*noWild.snp$log10.1.pval.),freq=F,breaks=50)
abline(h=1,col = "red")
layout(matrix(1:2,2,1))
plot(noWild.snp$XtXst)
plot(noWild.snp$log10.1.pval.,ylab="XtX P-value (-log10 scale)")
abline(h=3,lty=2,col = "red") #0.001 p--value theshold
```

```{r}
#Calibrating statistics with the simulation and analysis of 
# PODs (pseudoobserved data sets)

#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution
pi.beta.coef=read.table("noWild_core_summary_beta_params.out",h=T)$Mean
#upload the original data to obtain total allele count
noWild.data<-geno2YN("noWild.geno")
#Create the POD
simu.bta<-simulate.baypass(omega.mat=omega,nsnp=1000,sample.size=noWild.data$NN,
                           beta.pi=pi.beta.coef,pi.maf=0,suffix="btapods")
```

```{python}
## run baypass with the pseudoobserved data set
i_baypass -gfile G.btapods -outprefix noWild_pseudo 
```


```{r}
#######################################################
#Sanity Check: Compare POD and original data estimates
#######################################################
#get estimate of omega from the POD analysis
pod.omega=as.matrix(read.table("noWild_pseudo_mat_omega.out"))
plot(pod.omega,omega) ; abline(a=0,b=1,col = "red")
fmd.dist(pod.omega,omega)
#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table("noWild_pseudo_summary_beta_params.out",h=T)$Mean
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1,col = "red")

#######################################################
#XtX calibration
#######################################################
#get the pod XtX
pod.xtx=read.table("noWild_pseudo_summary_pi_xtx.out",h=T)$M_XtX
#compute the 1% threshold
pod.thresh=quantile(pod.xtx,probs=0.99)
#add the thresh to the actual XtX plot
plot(noWild.snp$M_XtX)
abline(h=pod.thresh,lty=2,col = "red")
```
#Similarly, when considering analysis of association with population-specific covariable under
#the core model, one may also calibrate of the different measures (BFis, regression coefficients,
#correlation coefficients, etc.) by analyzing a POD together with the covariables.
#More generally, the PODs distribution may also be used to compute empirical Pâ€“values and
#to derive from them q-values to control for multiple testing (see the qvalue package, Storey and Tibshirani, 2003).
```
###################################
# Used the on the covariate file
Ndia Hdia Hnod Ldia Lnod
300 1500 1500 350 350 # altitude
1 1 0 1 0 # binary diapause no diapause status
-2.9336 0.9231 0.9231 0.3796 0.3796 #PCA1 latitude longitude coordinates
#################################
# Analysis under the Importance Sampling (IS) covariate mode (MCMC is run under the core model before, see above)

# following the manual, let's test one variant at time 1-phenotype, 2-altitude, 3-latitude-longitude
# check the results

```{python}
i_baypass -gfile noWild.geno -efile noWild_phenotype.covariate -omegafile noWild_core_mat_omega.out -outprefix noWild_phenotype

i_baypass -gfile noWild.geno -efile noWild_altitude.covariate -omegafile noWild_core_mat_omega.out -outprefix noWild_altitude

i_baypass -gfile noWild.geno -efile noWild_coordinates.covariate -omegafile noWild_core_mat_omega.out -outprefix noWild_coordinates
```

# plot the results Diapause vs no diapause covariate A.K.A phenotype

```{r}
pheno.snp.res=read.table("noWild_phenotype_summary_betai_reg.out",h=T)
significant_markers <- subset(pheno.snp.res, BF.dB. >= 10) #select only markers with very strong bayes factor association
significant_markers <- subset(significant_markers, eBPis >= 0.99) #from the strong associated markes, select the ones with a high prob of inclusion 
print(significant_markers)


graphics.off()
layout(matrix(1:3,3,1))
plot(pheno.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
abline(h=10,lty=2,col = "red")
plot(pheno.snp.res$eBPis,xlab="SNP",ylab="eBPis")
abline(h =1, lty = 2, col = "red")
plot(pheno.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
abline(h=0,lty=2,col = "white")
```

# plot the results altitude covariate

```{r}
alt.snp.res=read.table("noWild_altitude_summary_betai_reg.out",h=T)
alt.loci <- read.table("loci_list.txt", header = FALSE, col.names = "Loci")
alt.snp.res <- as.data.frame(cbind(alt.snp.res, alt.loci))
significant_markers_alt <- subset(alt.snp.res, BF.dB. >= 10) #select only markers with very strong bayes factor association
significant_markers_alt <- subset(significant_markers_alt, eBPis >= 0.99) #from the strong associated markes, select the ones with a high prob of inclusion 
print(significant_markers_alt)

graphics.off()
layout(matrix(1:3,3,1))
plot(alt.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
abline(h=10,lty=2,col = "red")
plot(alt.snp.res$eBPis,xlab="SNP",ylab="eBPis")
abline(h = 3, lty = 2, col = "red")
plot(alt.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
abline(h=0,lty=2,col = "white")
```
significant_loci_alt <- significant_markers_alt$Loci
print(significant_loci_alt)
# plot the results latitude, longitude covariate

```{r}
graphics.off()
layout(matrix(1:3,3,1))
plot(coord.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
abline(h=10,lty=2,col = "red")
plot(coord.snp.res$eBPis,xlab="SNP",ylab="eBPis")
abline(h =3, lty = 2, col = "red")
plot(coord.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
abline(h=0,lty=2,col = "white")
```

# Now extract the significant snps list  for each setup
###############################################################################
# ALTITUDE
# lod baypass result and filter only SNPs over the 1% threshold
alt.snp.res=read.table("noWild_altitude_summary_betai_reg.out",h=T)
alt.loci <- read.table("loci_list.txt", header = FALSE, col.names = "Loci")
alt.snp.res <- as.data.frame(cbind(alt.snp.res, alt.loci))
significant_markers_alt <- subset(pheno.snp.res, BF.dB. >= 10) #select only markers with very strong bayes factor association
significant_markers_alt <- subset(significant_markers_alt, eBPis >= 0.99) #from the strong associated markes, select the ones with a high prob of inclusion 
print(significant_markers_alt)

significant_loci_alt <- significant_markers_alt$Loci
print(significant_loci_alt)
write.table(significant_loci_alt, file = "significant_loci_altitudeBF.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Extract the significant loci from original vcf file
# First need to edit the vcf file and remove all the comment lines before the header
# The first line should look like below:
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DP001	DP002	DP004	DP005	DP006	DP007	DP008	DP009	DP010	DP011	DP012	DP013	DP014	DP015	DP017	DP018	DP019	DP020	DP021	DP022	DP023	DP024	DP025	DP026	DP027	DP028	DP029	DP030	DP031	DP032	DP033	DP034	DP035	DP036	DP037	DP038	DP039	DP040	DP041	DP042	DP043	DP044	DP045	DP046	DP047	DP048	DP049	DP050	DP051	DP052	DP053	DP054	DP055	DP056	DP057	DP058	DP059	DP060	DP061	DP062	DP063	DP064	DP065	DP066	DP067	DP068	DP069	DP070	DP071	DP072	DP073	DP074	DP075	DP076	DP077	DP078	DP079	DP080	DP082	DP083	DP084	DP085	DP086	DP087	DP088	DP089	DP090	DP091	DP092	DP093	DP094	DP095	DP096	DP097	DP098	DP099	DP100	DP101	DP102	DP103	DP104	DP105	DP106	DP107	DP108	DP109	DP110	DP111	DP112	DP113	DP114	DP115	DP116	DP117	DP118	DP119	DP120	DP121	DP122	DP123	DP124	DP125	DP126	DP127	DP128	DP129	DP130	DP131	DP132	DP133	DP134	DP135	DP136	DP137	DP138	DP139	DP140	DP141	DP142	DP143	DP144	DP145	DP146	DP147	DP148	DP149	DP150	DP151	DP152	DP153	DP154	DP155	DP156	DP157	DP158	DP159	DP160	DP161	DP162	DP163	DP164	DP165	DP166	DP167	DP168	DP169	DP170	DP171	DP172	DP173	DP174	DP175	DP176	DP177	DP178	DP179	DP180	DP181	DP182	DP183	DP184	DP185	DP186	DP187	DP188	DP189	DP190	DP191	DP192	DP193	DP194	DP195	DP196	DP197	DP198	DP199	DP200	DP201	DP202	DP203	DP204	DP205	DP206	DP207	DP208	DP209	DP210	DP211	DP212	DP213	DP214	DP215	DP216	DP217	DP218	DP219	DP220	DP221	DP222	DP223	DP224	DP225	DP226	DP227	DP228	DP230	DP231	DP232	DP233	DP234	DP235	DP236	DP237	DP238	DP239	DP240	DP241	DP242	DP243	DP244	DP245	DP246	DP247	DP248	DP249	DP250	DP251	DP252	DP253	DP254	DP255	DP256	DP257	DP258	DP259	DP260	DP261	DP262	DP263	DP264	DP265	DP266	DP267	DP268	DP269	DP270	DP271	DP272	DP273	DP274	DP275	DP276	DP277	DP278	DP279	DP280	DP281	DP282	DP283	DP284	DP285	DP286	DP287	DP288	DP289	DP290	DP291	DP292	DP293	DP294	DP295	DP296	DP297	DP298	DP299	DP300

# Load necessary libraries
library(data.table)

# Load the VCF file as a data.table for efficient processing
vcf_file <- fread("noWild2.vcf", header = TRUE)

# Load the significant loci IDs from the text file
significant_loci_alt <- readLines("significant_loci_altitudeBF.txt")

# Filter the VCF data to include only rows with matching IDs from the text file
filtered_vcf <- vcf_file[ID %in% significant_loci_alt]

# Display the filtered VCF rows
print(filtered_vcf)

# Save the filtered results to a new text file
fwrite(filtered_vcf, "filtered_noWild2_significant_altitudeBF.vcf", sep = "\t")

# mapping loci of vcf file on the updated genome annotation
# The genome annotation from D Powell et al. 2021 was updated to the Krystyna genome version with 65 contigs.
# the file Ityp65.gff3 contain the updated annotation done with liftoff on the hpc cluster using the following commands.
```{bash}
 cd /data/users/lpalmierirocha/genome/Ips_typographus
 module load anaconda3/
 conda activate liftoff

 liftoff -g Ityp.gff3 -o Ityp65.gff3 Krystyna_ips.fasta ncbi_genome_renamed.fasta
```
# now to extracted the mapped loci from the updated annotation
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")
BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)

# Load vcf file
#vcf_data <- read.table("filtered_noWild2_significant_altitude.txt", header = TRUE)
vcf_data <- filtered_vcf
# Check the data types of the key columns
str(vcf_data$CHROM)
str(vcf_data$POS)
str(vcf_data$ID)
# Ensure no missing or unusual values are in these columns
summary(vcf_data$CHROM)
summary(vcf_data$POS)
summary(vcf_data$ID)

# obtain the ranges
vcf_gr <- GRanges(seqnames = vcf_data$CHROM,
                  ranges = IRanges(start = vcf_data$POS, end = vcf_data$POS),
                  names = vcf_data$ID)

# Check the length and contents of vcf_gr
length(vcf_gr)
vcf_gr

# Load gene annotations
gff_data <- import("Ityp65.gff3", format = "gff3")

# find the overlaps between vcf loci and gff3 annotations
overlaps <- findOverlaps(vcf_gr, gff_data)

# Extract information about the matches
vcf_hits <- vcf_gr[queryHits(overlaps)]
gff_hits <- gff_data[subjectHits(overlaps)]

# Combine information for inspection
hits_data <- data.frame(
  Locus = mcols(vcf_hits)$names,
  Chromosome = seqnames(vcf_hits),
  Position = start(vcf_hits),
  FeatureType = mcols(gff_hits)$type,
  FeatureID = mcols(gff_hits)$ID,
  FeatureName = mcols(gff_hits)$Name
)

# Display the results
head(hits_data)
# output the results to a file
write.table(hits_data, file = "noWild_ALTITUDE_gff65_overlaps.txt", sep = "\t", row.names = FALSE, quote = FALSE)
##############################################################################
# PHENOTYPE

pheno.loci <- read.table("loci_list.txt", header = FALSE, col.names = "Loci")
pheno.snp.res <- as.data.frame(cbind(pheno.snp.res, pheno.loci))

# lod baypass result and filter only SNPs over the 1% threshold
pheno.snp.res=read.table("noWild_phenotype_summary_betai_reg.out",h=T)
pheno.loci <- read.table("loci_list.txt", header = FALSE, col.names = "Loci")
pheno.snp.res <- as.data.frame(cbind(pheno.snp.res, pheno.loci))
significant_markers_pheno <- subset(pheno.snp.res, BF.dB. >= 10) #select only markers with very strong bayes factor association
significant_markers_pheno <- subset(significant_markers_pheno, eBPis >= 0.99) #from the strong associated markes, select the ones with a high prob of inclusion 
print(significant_markers_pheno)

significant_loci_pheno <- significant_markers_pheno$Loci
print(significant_loci_pheno)
write.table(significant_loci_pheno, file = "significant_loci_phenotypeBF.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Extract the significant loci from original vcf file
# First need to edit the vcf file and remove all the comment lines before the header
# The first line should look like below:
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DP001	DP002	DP004	DP005	DP006	DP007	DP008	DP009	DP010	DP011	DP012	DP013	DP014	DP015	DP017	DP018	DP019	DP020	DP021	DP022	DP023	DP024	DP025	DP026	DP027	DP028	DP029	DP030	DP031	DP032	DP033	DP034	DP035	DP036	DP037	DP038	DP039	DP040	DP041	DP042	DP043	DP044	DP045	DP046	DP047	DP048	DP049	DP050	DP051	DP052	DP053	DP054	DP055	DP056	DP057	DP058	DP059	DP060	DP061	DP062	DP063	DP064	DP065	DP066	DP067	DP068	DP069	DP070	DP071	DP072	DP073	DP074	DP075	DP076	DP077	DP078	DP079	DP080	DP082	DP083	DP084	DP085	DP086	DP087	DP088	DP089	DP090	DP091	DP092	DP093	DP094	DP095	DP096	DP097	DP098	DP099	DP100	DP101	DP102	DP103	DP104	DP105	DP106	DP107	DP108	DP109	DP110	DP111	DP112	DP113	DP114	DP115	DP116	DP117	DP118	DP119	DP120	DP121	DP122	DP123	DP124	DP125	DP126	DP127	DP128	DP129	DP130	DP131	DP132	DP133	DP134	DP135	DP136	DP137	DP138	DP139	DP140	DP141	DP142	DP143	DP144	DP145	DP146	DP147	DP148	DP149	DP150	DP151	DP152	DP153	DP154	DP155	DP156	DP157	DP158	DP159	DP160	DP161	DP162	DP163	DP164	DP165	DP166	DP167	DP168	DP169	DP170	DP171	DP172	DP173	DP174	DP175	DP176	DP177	DP178	DP179	DP180	DP181	DP182	DP183	DP184	DP185	DP186	DP187	DP188	DP189	DP190	DP191	DP192	DP193	DP194	DP195	DP196	DP197	DP198	DP199	DP200	DP201	DP202	DP203	DP204	DP205	DP206	DP207	DP208	DP209	DP210	DP211	DP212	DP213	DP214	DP215	DP216	DP217	DP218	DP219	DP220	DP221	DP222	DP223	DP224	DP225	DP226	DP227	DP228	DP230	DP231	DP232	DP233	DP234	DP235	DP236	DP237	DP238	DP239	DP240	DP241	DP242	DP243	DP244	DP245	DP246	DP247	DP248	DP249	DP250	DP251	DP252	DP253	DP254	DP255	DP256	DP257	DP258	DP259	DP260	DP261	DP262	DP263	DP264	DP265	DP266	DP267	DP268	DP269	DP270	DP271	DP272	DP273	DP274	DP275	DP276	DP277	DP278	DP279	DP280	DP281	DP282	DP283	DP284	DP285	DP286	DP287	DP288	DP289	DP290	DP291	DP292	DP293	DP294	DP295	DP296	DP297	DP298	DP299	DP300

# Load necessary libraries
library(data.table)

# Load the VCF file as a data.table for efficient processing
vcf_file <- fread("noWild2.vcf", header = TRUE)

# Load the significant loci IDs from the text file
significant_loci_pheno <- readLines("significant_loci_phenotypeBF.txt")

# Filter the VCF data to include only rows with matching IDs from the text file
filtered_vcf <- vcf_file[ID %in% significant_loci_pheno]

# Display the filtered VCF rows
print(filtered_vcf)

# Save the filtered results to a new text file
fwrite(filtered_vcf, "filtered_noWild2_significant_phenotypeBF.txt", sep = "\t")

# now to extracted the mapped loci from the updated annotation
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")
BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)

# Load vcf file
#vcf_data <- read.table("filtered_noWild2_significant_phenotype.txt", header = TRUE)
vcf_data <- filtered_vcf
# Check the data types of the key columns
str(vcf_data$CHROM)
str(vcf_data$POS)
str(vcf_data$ID)
# Ensure no missing or unusual values are in these columns
summary(vcf_data$CHROM)
summary(vcf_data$POS)
summary(vcf_data$ID)

# obtain the ranges
vcf_gr <- GRanges(seqnames = vcf_data$CHROM,
                  ranges = IRanges(start = vcf_data$POS, end = vcf_data$POS),
                  names = vcf_data$ID)

# Check the length and contents of vcf_gr
length(vcf_gr)
vcf_gr

# Load gene annotations
gff_data <- import("Ityp65.gff3", format = "gff3")

# find the overlaps between vcf loci and gff3 annotations
overlaps <- findOverlaps(vcf_gr, gff_data)

# Extract information about the matches
vcf_hits <- vcf_gr[queryHits(overlaps)]
gff_hits <- gff_data[subjectHits(overlaps)]

# Combine information for inspection
hits_data <- data.frame(
  Locus = mcols(vcf_hits)$names,
  Chromosome = seqnames(vcf_hits),
  Position = start(vcf_hits),
  FeatureType = mcols(gff_hits)$type,
  FeatureID = mcols(gff_hits)$ID,
  FeatureName = mcols(gff_hits)$Name
)

# Display the results
head(hits_data)
# output the results to a file
write.table(hits_data, file = "noWild_PHENOTYPE_gff65_overlaps.txt", sep = "\t", row.names = FALSE, quote = FALSE)
##################################################################
# LATITUDE

coord.loci <- read.table("loci_list.txt", header = FALSE, col.names = "Loci")
coord.snp.res <- as.data.frame(cbind(coord.snp.res, coord.loci))

# lod baypass result and filter only SNPs over the 1% threshold
coord.snp.res=read.table("noWild_coordinates_summary_betai_reg.out",h=T)
coord.loci <- read.table("loci_list.txt", header = FALSE, col.names = "Loci")
coord.snp.res <- as.data.frame(cbind(coord.snp.res, coord.loci))
significant_markers_coord <- subset(coord.snp.res, BF.dB. >= 10) #select only markers with very strong bayes factor association
significant_markers_coord <- subset(significant_markers_coord, eBPis >= 0.99) #from the strong associated markes, select the ones with a high prob of inclusion 
print(significant_markers_coord)
significant_loci_coord <- significant_markers_coord$Loci
print(significant_loci_coord)
write.table(significant_loci_coord, file = "significant_loci_latitudeBF.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Extract the significant loci from original vcf file
# First need to edit the vcf file and remove all the comment lines before the header
# The first line should look like below:
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DP001	DP002	DP004	DP005	DP006	DP007	DP008	DP009	DP010	DP011	DP012	DP013	DP014	DP015	DP017	DP018	DP019	DP020	DP021	DP022	DP023	DP024	DP025	DP026	DP027	DP028	DP029	DP030	DP031	DP032	DP033	DP034	DP035	DP036	DP037	DP038	DP039	DP040	DP041	DP042	DP043	DP044	DP045	DP046	DP047	DP048	DP049	DP050	DP051	DP052	DP053	DP054	DP055	DP056	DP057	DP058	DP059	DP060	DP061	DP062	DP063	DP064	DP065	DP066	DP067	DP068	DP069	DP070	DP071	DP072	DP073	DP074	DP075	DP076	DP077	DP078	DP079	DP080	DP082	DP083	DP084	DP085	DP086	DP087	DP088	DP089	DP090	DP091	DP092	DP093	DP094	DP095	DP096	DP097	DP098	DP099	DP100	DP101	DP102	DP103	DP104	DP105	DP106	DP107	DP108	DP109	DP110	DP111	DP112	DP113	DP114	DP115	DP116	DP117	DP118	DP119	DP120	DP121	DP122	DP123	DP124	DP125	DP126	DP127	DP128	DP129	DP130	DP131	DP132	DP133	DP134	DP135	DP136	DP137	DP138	DP139	DP140	DP141	DP142	DP143	DP144	DP145	DP146	DP147	DP148	DP149	DP150	DP151	DP152	DP153	DP154	DP155	DP156	DP157	DP158	DP159	DP160	DP161	DP162	DP163	DP164	DP165	DP166	DP167	DP168	DP169	DP170	DP171	DP172	DP173	DP174	DP175	DP176	DP177	DP178	DP179	DP180	DP181	DP182	DP183	DP184	DP185	DP186	DP187	DP188	DP189	DP190	DP191	DP192	DP193	DP194	DP195	DP196	DP197	DP198	DP199	DP200	DP201	DP202	DP203	DP204	DP205	DP206	DP207	DP208	DP209	DP210	DP211	DP212	DP213	DP214	DP215	DP216	DP217	DP218	DP219	DP220	DP221	DP222	DP223	DP224	DP225	DP226	DP227	DP228	DP230	DP231	DP232	DP233	DP234	DP235	DP236	DP237	DP238	DP239	DP240	DP241	DP242	DP243	DP244	DP245	DP246	DP247	DP248	DP249	DP250	DP251	DP252	DP253	DP254	DP255	DP256	DP257	DP258	DP259	DP260	DP261	DP262	DP263	DP264	DP265	DP266	DP267	DP268	DP269	DP270	DP271	DP272	DP273	DP274	DP275	DP276	DP277	DP278	DP279	DP280	DP281	DP282	DP283	DP284	DP285	DP286	DP287	DP288	DP289	DP290	DP291	DP292	DP293	DP294	DP295	DP296	DP297	DP298	DP299	DP300

# Load necessary libraries
library(data.table)

# Load the VCF file as a data.table for efficient processing
vcf_file <- fread("noWild2.vcf", header = TRUE)

# Load the significant loci IDs from the text file
significant_loci_coord <- readLines("significant_loci_latitudeBF.txt")

# Filter the VCF data to include only rows with matching IDs from the text file
filtered_vcf <- vcf_file[ID %in% significant_loci_coord]

# Display the filtered VCF rows
print(filtered_vcf)

# Save the filtered results to a new text file
fwrite(filtered_vcf, "filtered_noWild2_significant_latitudeBF.txt", sep = "\t")

# now to extracted the mapped loci from the updated annotation
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")
BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)

# Load vcf file
#vcf_data <- read.table("filtered_noWild2_significant_latitude.txt", header = TRUE)
vcf_data <- filtered_vcf
# Check the data types of the key columns
str(vcf_data$CHROM)
str(vcf_data$POS)
str(vcf_data$ID)
# Ensure no missing or unusual values are in these columns
summary(vcf_data$CHROM)
summary(vcf_data$POS)
summary(vcf_data$ID)

# obtain the ranges
vcf_gr <- GRanges(seqnames = vcf_data$CHROM,
                  ranges = IRanges(start = vcf_data$POS, end = vcf_data$POS),
                  names = vcf_data$ID)

# Check the length and contents of vcf_gr
length(vcf_gr)
vcf_gr

# Load gene annotations
gff_data <- import("Ityp65.gff3", format = "gff3")

# find the overlaps between vcf loci and gff3 annotations
overlaps <- findOverlaps(vcf_gr, gff_data)

# Extract information about the matches
vcf_hits <- vcf_gr[queryHits(overlaps)]
gff_hits <- gff_data[subjectHits(overlaps)]

# Combine information for inspection
hits_data <- data.frame(
  Locus = mcols(vcf_hits)$names,
  Chromosome = seqnames(vcf_hits),
  Position = start(vcf_hits),
  FeatureType = mcols(gff_hits)$type,
  FeatureID = mcols(gff_hits)$ID,
  FeatureName = mcols(gff_hits)$Name
)

# Display the results
head(hits_data)
# output the results to a file
write.table(hits_data, file = "noWild_LATITUDE_gff65_overlaps.txt", sep = "\t", row.names = FALSE, quote = FALSE)
##################################################################
```{r}
## we will also make a Venn diagram to check the common SNPs between the three comparisons
library(ggvenn)
library(grid)

phenotype = significant_markers_pheno
altitude = significant_markers_alt
coordinates = significant_markers_coord

x <- list ('Phenotype' = phenotype$VLoci,
           'Altitude' = altitude$Loci,
           'Coordinates' = coordinates$Loci)

venn_plot <- ggvenn(x, c("Phenotype", "Altitude", "Coordinates"), fill_color = c("#b35806", "#f1a340", "#fee0b6", "#d8daeb", "#998ec3", "#542788"), fill_alpha = 0.7, stroke_size = 0.2, set_name_size = 4, stroke_color = "black", show_percentage = FALSE)

# Save the plot as a PNG file
#ggsave("venn_diagram_noWild_signif_markers.png", plot = venn_plot, width = 10, height = 10)
ggsave("venn_diagram_noWild_signif_markers.pdf", plot = venn_plot, width = 10, height = 10)
```
=================================================
library(ggvenn)
library(grid)

# Preprocess data
clean_column <- function(column) {
  # Ensure all values are characters, trim whitespace, and convert to lowercase
  unique(tolower(trimws(as.character(column))))
}

phenotype_clean <- clean_column(significant_markers_pheno$V11)
altitude_clean <- clean_column(significant_markers_alt$V11)
coordinates_clean <- clean_column(significant_markers_coord$V11)

x <- list(
  'Phenotype' = phenotype_clean,
  'Altitude' = altitude_clean,
  'Coordinates' = coordinates_clean
)

# Create the Venn diagram
venn_plot <- ggvenn(
  x, 
  c("Phenotype", "Altitude", "Coordinates"), 
  fill_color = c("#b35806", "#f1a340", "#fee0b6", "#d8daeb", "#998ec3", "#542788"), 
  fill_alpha = 0.7, 
  stroke_size = 0.2, 
  set_name_size = 4, 
  stroke_color = "black", 
  show_percentage = FALSE
)

# Save the plot as a PNG file
ggsave("venn_diagram_noWild_signif_markersclean.png", plot = venn_plot, width = 10, height = 10)



