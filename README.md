# PanGenome_VCF_PostProcess
A simple pipeline to process the VCF deconstructed from the VG pangenome graph 
1. Remove large variants with `vcfbub -l 0 -r 10000000`;
2. Trim ALT alleles not seen in the genotype fields with `bcftools view -a`;
3. Group the alternative alleles in multiallelic SV loci by their length (see `allele-group.pdf`);
4. Split multiallelics with `bcftools norm -m -any`;
5. Classify the variants into small variants and SVs (abs(length(ALT_allele)-length(REF_allele))>=50). 

# Usage
### 1. Installation
conda or mamba is required
```
# clone the repo
git clone https://github.com/Shuhua-Group/PanGenome_VCF_PostProcess.git
# create the enviroment
cd PanGenome_VCF_PostProcess
mamba env create -f environment.yaml
# activate the enviroment
mamba activate PanVCF
```
### 2. Modify the configfile
Add the working directory `work_dir:` and PanGenome VCF file `pan_vcf` to `config.yaml` 
```
# where you want the pipeline work
work_dir: ${wd}
# the vcf file deconstructed from the vg pangenome 
pan_vcf: ${PanGenome_VCF}
# prefix of output files 
out_prefix: CPC.HPRC.Phase1.processed
# allele grouping threads
threads: 16
```
### 3. Run the pipeline
Make sure that `PanGenome_vcf.post_process.smk`, `config.yaml`, and `module.py` are under the same directory:
```
snakemake -c 128 -s PanGenome_vcf.post_process.smk
```
It cost about two hours in a 128 core sever to process the [CPC+HPRC pangenome](https://pog.fudan.edu.cn/cpc/files/CPC.HPRC.Phase1.CHM13v2/CPC.HPRC.Phase1.CHM13v2.vcf.gz), and the processed VCF files are:
```
# biallelics small variants (abs(length(ALT_allele)-length(REF_allele))<50)
CPC.HPRC.Phase1.processed.small_variants.normed.vcf.gz
# biallelics SVs (abs(length(ALT_allele)-length(REF_allele))>=50)
CPC.HPRC.Phase1.processed.SVs.normed.vcf.gz
```
The allele correspondence before and after the merging is output as the following file:
```
temp/{chr}/temp.SVs.merge_map.csv
```

**[Note]** We have updated the process strategy about the complex loci with both small variants and SVs alleles, so that the number of variants is slightly different than in the paper.
