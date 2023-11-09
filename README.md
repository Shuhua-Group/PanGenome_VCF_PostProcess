# PanGenome_VCF_PostProcess
A simple pipeline to process the VCF deconstructed from the VG pangenome graph 
1. Remove large variants with `vcfbub -l 0 -r 10000000`;
2. Trim ALT alleles not seen in the genotype fields with `bcftools view -a`;
3. Group the alternative alleles in multiallelic SV loci by their length (see `allele_group.pdf`);
4. Split multiallelics with `bcftools norm -m -any`;
5. Classify the variants into small variants and SVs (abs(length(ALT_allele)-length(REF_allele))>=50).
