# Fst_for_Laevis_and_Gilli

# List into comma seperated values
```text
Open Word
"Paste special" as text only
Select the data in Word (the one that you need to convert to text separated with ,), press Ctrl-H (Find & replace)
In "Find what" box type ^p
In "Replace with" box type ,
Select "Replace all"
```

# Filter Gilli samples only from the VCF file
```bash
bcftools view -s  XGL713_179_sorted.bam,XG12_07_sorted.bam,XGUAE_43_sorted.bam,XGUAE_42_sorted.bam,XGUAE_36_sorted.bam,XG92_sorted.bam,XGL713_123_sorted.bam,XGUAE_44_sorted.bam,XG153_sorted.bam,XGL713_177_sorted.bam,XGL713_180_sorted.bam,XGL713_181_sorted.bam Xlaevis_and_gilli_all_samples_merged_sorted.bam.vcf.gz > Gilli_only.vcf
```
# Filter Laevis samples only
```bash
bcftools view -s XGUAE_71_sorted.bam,XGL713_179_sorted.bam,XG12_07_sorted.bam,XGUAE_93_sorted.bam,XGUAE_92_sorted.bam,BJE267_sorted.bam,BJE263_sorted.bam,XGL713_232_sorted.bam,BJE3608.fqsorted.bam,BJE3639_sorted.bam,jonk_02.fqsorted.bam,XGUAE_43_sorted.bam,BJE3545_sorted.bam,XGUAE_65_sorted.bam,BJE266_sorted.bam,BJE261_sorted.bam,XGUAE_72_sorted.bam,XGUAE_42_sorted.bam,XGUAE_97_sorted.bam,XGUAE_36_sorted.bam,XLJONK_14_sorted.bam,CoGH105.fqsorted.bam,BJE264_sorted.bam,XGUAE_59_sorted.bam,XG92_sorted.bam,XL_CPT1_sorted.bam,BJE265_sorted.bam,XGUAE_124_sorted.bam,BJE3574.fqsorted.bam,XGL713_123_sorted.bam,XGUAE_44_sorted.bam,JMEC006.fqsorted.bam,BJE3581.fqsorted.bam,XG153_sorted.bam,BJE3536.fqsorted.bam,XGL713_177_sorted.bam,XGL713_180_sorted.bam,BJE1488_sorted.bam,BJE1489_sorted.bam,XGL713_181_sorted.bam,XL_CPT4_sorted.bam,XL_CPT2_sorted.bam,JMEC003.fqsorted.bam,XGUAE_70_sorted.bam,BJE3540.fqsorted.bam,XL_CPT3_sorted.bam Xlaevis_and_gilli_all_samples_merged_sorted.bam.vcf.gz > Laevis_only.vcf
```


