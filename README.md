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

```bash
library(hierfstat)
library(data.table)



my_samples<-read.fstat("01_radiator_genomic_converter_20201115@1303/radiator_data_20201115@1303_hierfstat.dat")

pairwise_all_chrs<-boot.ppfst(dat = my_samples,nboot = 1000,quant=c(0.025,0.975))
print("done calculating pairwise fst for whole genome")

#making a proper/writable dataframe from those data / pairwise_all_chrs is not writable

x<-pairwise_all_chrs$ll
y<-pairwise_all_chrs$ul
no_of_pops<-sqrt(length(x))

#lower limit values
df1 <- data.frame(matrix(unlist(x), nrow=sqrt(length(x)), byrow=T),stringsAsFactors=FALSE)

# Create full table for lower limit by duplicating values(to use values again as 1:2 and 2:1)
# this fills df1 with lower limit values
a<-1
load
for (j in 1:no_of_pops) {
          for (i in a:no_of_pops) {
                      df1[a,i]<-df1[i,a]
  }
  a<-a+1
}

#upper limit values
df2 <- data.frame(matrix(unlist(y), nrow=sqrt(length(y)), byrow=T),stringsAsFactors=FALSE)
# this fills df2 with upper limit values
a<-1
load
for (j in 1:no_of_pops) {
          for (i in a:no_of_pops) {
                      df2[a,i]<-df2[i,a]
  }
  a<-a+1
}

# output filled tables as CSV
write.csv(df1,file = "all_chrs_lower_limit_ppfst.csv")

write.csv(df2,file = "all_chrs_upper_limit_ppfst.csv")
```

# plot_ppfst_values.r

```r
#set current path as wd
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(data.table)

# load data all chrs
all_chrs_lower_limit<-read.csv("./data/all_chrs/all_chrs_lower_limit_ppfst.csv")
all_chrs_upper_limit<-read.csv("./data/all_chrs/all_chrs_upper_limit_ppfst.csv")

# load data L only(disregard the file name "all_chrs", didnt change file name in script. But they are filtered for L and S)
l_only_lower_limit<-read.csv("./data/l_only/all_chrs_lower_limit_ppfst.csv")
l_only_upper_limit<-read.csv("./data/l_only/all_chrs_upper_limit_ppfst.csv")

# load data S only(disregard the file name "all_chrs", didnt change file name in script. But they are filtered for L and S)
s_only_lower_limit<-read.csv("./data/s_only/all_chrs_lower_limit_ppfst.csv")
s_only_upper_limit<-read.csv("./data/s_only/all_chrs_upper_limit_ppfst.csv")


plot_list<-list()

p<-2


#colors for all,L,S in that order
plot_colours<-c("black","red","lightblue")

plot_no<-1


for (j in 1:10) {
  
  n<-1
  
  emp_df_x<-1:10
  emp_df_y<-c(0,0,0,0,0,0,0,0,0,0)
  emp_df<-data.frame(emp_df_x,emp_df_y)
  
  emp_plot<-ggplot(data = emp_df,aes(x = emp_df_x,y = emp_df_y))+
    scale_x_continuous("Population", labels = as.character(emp_df_x), breaks = emp_df_x)+
    theme_bw()
  
  
  for (i in 1:10) {
    x<-c(n,n,n)
    
    # order here from all data is all,L,S, therefore with x_for errorbars values assigned in the following line, n-0.25/S comes left , n/all comes middles, n+0.25/L comes last
    x_for_errorbars<-c(n,n+0.25,n-0.25)
    
    # get mean fst values for plot in the middle
    y<-c((all_chrs_lower_limit[n,p]+all_chrs_upper_limit[n,p])/2,(l_only_lower_limit[n,p]+l_only_upper_limit[n,p])/2,(s_only_lower_limit[n,p]+s_only_upper_limit[n,p])/2)
    
    id<-1:3
    
    # get lower limit fst for y
    ci_l<-c(all_chrs_lower_limit[n,p],l_only_lower_limit[n,p],s_only_lower_limit[n,p])
    
    # get upper limit for y
    ci_u<-c(all_chrs_upper_limit[n,p],l_only_upper_limit[n,p],s_only_upper_limit[n,p])
    
    data<-data.frame(x,x_for_errorbars,y,id,ci_l,ci_u)
    
    
    if (n==1) {
      plot_set<-emp_plot+
        geom_point(data = data,aes(x = x_for_errorbars,y=y),colour=plot_colours,shape=4)+
        geom_errorbar(data = data, mapping = aes(x = x_for_errorbars, y = y, ymin =ci_l, ymax = ci_u), size=1, color=plot_colours, width=.1) +
        geom_boxplot(data = data,aes(x = x,y = y))+
        ylab(paste("pop",plot_no))+
        ylim(0,1)+
        theme_bw()
    
    }
  
    if (n!=1) {
      plot_set<-plot_set+
        geom_point(data = data,aes(x = x_for_errorbars,y=y),colour=plot_colours,shape=4)+
        geom_errorbar(data = data, mapping = aes(x = x_for_errorbars, y = y, ymin =ci_l, ymax = ci_u), size=1, color=plot_colours, width=.1) +
        geom_boxplot(data = data,aes(x = x,y = y))+
        ylab(paste("pop",plot_no))+
        ylim(0,1)+
        theme_bw()
    }
    
    
    n<-n+1
  }
  
  plot_list[[plot_no]]=plot_set
  
  
  p<-p+1
  
  
  plot_no<-plot_no+1
  
}

library(gridExtra)
Final_plot_grid<-grid.arrange(
  grobs = plot_list,
  ncol=3
  
  #widths = c(2, 2),
  #heights=c(2,1)
)
  
ggsave(filename = "ppfst.pdf",plot = Final_plot_grid,height = 45,width =25 )


print("DONE")
print("all=black,L=red,S=Green")
```
# Run cal_ppfst.sh

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
module load gcc/7.3.0
module load r
R

Rscript cal_ppfst.r
```
# cal_fst_without_ci.r

```bash
library(hierfstat)
library(data.table)



my_samples<-read.fstat("01_radiator_genomic_converter_20201115@1303/radiator_data_20201115@1303_hierfstat.dat")

fst<-basic.stats(my_samples,diploid=TRUE)
print("done calculating pairwise fst for whole genome")


final<-fst$overall
write.csv(final,"global_fst_without_ci.csv")
```
# run_cal_fst_without_ci.r

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
module load gcc/7.3.0
module load r
R

Rscript cal_fst_without_ci.r
```
# cal_cal_golbal_fst_with_ci.r

```bash
library(hierfstat)
dat<-read.fstat("01_radiator_genomic_converter_20201115@1303/radiator_data_20201115@1303_hierfstat.dat")
 
#get FST FIS and Variance components
x<-wc(dat)
 
#function to bootstrap the variance components
get.boot<-function(x){
        a<-sample(nrow(x),replace=TRUE)
         
        tmp<-colSums(x[a,])
        FST<-tmp[1]/sum(tmp)
        FIS<-tmp[2]/sum(tmp[2:3])
        c(FIS,FST)
}
 
#1000 bootstraps
boot.fst<-replicate(1000,get.boot(x$sigma.loc))
row.names(boot.fst)<-c("FIS","FST")
 
#get empirical quantiles
fst<-apply(boot.fst,1,quantile,c(0.025,0.975))
 
write.csv(fst,file = "fst_global_with_ci.csv"
```




