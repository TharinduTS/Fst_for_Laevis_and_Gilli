# Fst_for_Laevis_and_Gilli

# Filter Gilli samples only from the VCF file
```bash
bcftools view -s  XGL713_179_sorted.bam,XG12_07_sorted.bam,XGUAE_43_sorted.bam,XGUAE_42_sorted.bam,XGUAE_36_sorted.bam,XG92_sorted.bam,XGL713_123_sorted.bam,XGUAE_44_sorted.bam,XG153_sorted.bam,XGL713_177_sorted.bam,XGL713_180_sorted.bam,XGL713_181_sorted.bam Xlaevis_and_gilli_all_samples_merged_sorted.bam.vcf.gz > Gilli_only.vcf
```
# List into comma seperated values
```text
Open Word
"Paste special" as text only
Select the data in Word (the one that you need to convert to text separated with ,), press Ctrl-H (Find & replace)
In "Find what" box type ^p
In "Replace with" box type ,
Select "Replace all"
```
