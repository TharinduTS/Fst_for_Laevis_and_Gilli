# Fst_for_Laevis_and_Gilli

# subset variant sites 
```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL


module load nixpkgs/16.09 
module load intel/2018.3
module load vcftools/0.1.16

vcftools --gzvcf Xlaevis_and_gilli_all_samples_merged_sorted.bam.vcf.gz --keep all_sample_list --non-ref-ac-any 1 --recode --recode-INFO-all --stdout > all_sample_varient_sites.vcf
```

# filter vcf keeping just variant sites in chromosomes(removing scaffolds) from all samples

```bash

vcftools --gzvcf ../Xlaevis_and_gilli_all_samples_merged_sorted.bam.vcf.gz --keep ../all_sample_list --chr chr1L --chr chr1S --chr chr2L --chr chr2S --chr chr3L --chr chr3S --chr chr4L --chr chr4S --chr chr5L --chr chr5S --chr chr6L --chr chr6S --chr chr7L --chr chr7S --chr chr8L --chr chr8S --chr chr9_10L --chr chr9_10S --non-ref-ac-any 1 --recode --recode-INFO-all --stdout > all_samples_varient_chrs_only
```
# remove scaffold info from header
```bash
grep -v '^##contig=<ID=Scaffold' all_samples_varient_chrs_only > all_samples_varient_chrs_only_header_scaffolds_cleared
```

# filter L only from all
```bash
vcftools --vcf ./all_samples_varient_chrs_only_header_scaffolds_cleared.vcf --keep ../populations_and_species/all_sample_list --chr chr1L  --chr chr2L --chr chr3L --chr chr4L --chr chr5L  --chr chr6L  --chr chr7L  --chr chr8L  --chr chr9_10L  --recode --recode-INFO-all --stdout > L_only_varient_chrs_only_header_scaffolds_cleared.vcf
```
# filter S only from all
```bash
vcftools --vcf ./all_samples_varient_chrs_only_header_scaffolds_cleared.vcf --keep ../populations_and_species/all_sample_list --chr chr1S  --chr chr2S --chr chr3S --chr chr4S --chr chr5S  --chr chr6S  --chr chr7S  --chr chr8S  --chr chr9_10S  --recode --recode-INFO-all --stdout > S_only_varient_chrs_only_header_scaffolds_cleared.vcf
```
# for laevis only, filtered laevis samples only from all of three above genomes

# all chrs
```bash
vcftools --vcf all_samples_varient_chrs_only_header_scaffolds_cleared.vcf --keep ../populations_and_species/laevis_sample_list --recode --recode-INFO-all --out ./laevis_varient_chrs_only_header_scaffolds_cleared.vcf
```
# L only
```bash
vcftools --vcf L_only_varient_chrs_only_header_scaffolds_cleared.vcf --keep ../populations_and_species/laevis_sample_list --recode --recode-INFO-all --out ./laevis_L_only_varient_chrs_only_header_scaffolds_cleared.vcf
```
#S only
```bash
vcftools --vcf S_only_varient_chrs_only_header_scaffolds_cleared.vcf --keep ../populations_and_species/laevis_sample_list --recode --recode-INFO-all --out ./laevis_S_only_varient_chrs_only_header_scaffolds_cleared.vcf
```


# keep different subgenomes in different folders when proceeding so you can use same codes in seperate folders

# convert data to hierfstat input format
``` bash
module load nixpkgs/16.09
module load gcc/7.3.0
module load r
R

library(radiator)

test_all_pop <- radiator::genomic_converter(
                                              data = "../../vcfs/laevis_varient_chrs_only_header_scaffolds_cleared.vcf.recode.vcf", 
                                                strata = "../../strata/strata_pops_L_only.tsv", 
                                                output = "hierfstat"
                                                )

```

# installing R packages in computecanada

# install packages first  IN A PERSONAL LIBRARY

# Change libPath (try this if package installation gives you errors)
```bash
export R_LIBS="~/R/lib"
```

```bash
module load nixpkgs/16.09 gcc/7.3.0 r/4.0.0 gdal/3.0.1
module load proj/6.3.0 udunits/2.2.26

R
install.packages("devtools")
library(devtools)

install_github("jgx65/hierfstat”)
library(hierfstat)

```

# Rscript for calculating basic stats(including Fst without CI)

# load R and run
```bash
library(hierfstat)
library(data.table)



my_samples<-read.fstat("./01_radiator_genomic_converter_20201115@1142/radiator_data_20201115@1142_hierfstat.dat")

fst<-basic.stats(my_samples,diploid=TRUE)
print("done calculating pairwise fst for whole genome")


final<-fst$overall
write.csv(final,"global_fst_without_ci.csv")
```
# Calculating ppfst/ creating usable data table and saving as CSV


