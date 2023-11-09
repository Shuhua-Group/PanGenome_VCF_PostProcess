import sys,os
import gzip
from collections import Counter
import pandas as pd
from multiprocessing import Pool

from src import *

configfile: "config.yaml"
pan_vcf = config['pan_vcf']
out_prefix = config['out_prefix']
wd = config['work_dir'] + '/'
threads = int(config['threads'])

CHROMOSOMES = ['chr'+i for i in list(map(str,range(1,23)))+["X","Y","M"]]


rule all:
    input:
        SVs_vcf = wd + out_prefix + '.SVs.normed.vcf.gz',
        SVs_vcf_tbi = wd + out_prefix + '.SVs.normed.vcf.gz.tbi',        
        small_variants_vcf = wd + out_prefix + '.small_variants.normed.vcf.gz',
        small_variants_vcf_tbi = wd + out_prefix + '.small_variants.normed.vcf.gz.tbi',

rule pan_vcf_to_chromosomes:
    input:
        pan_vcf = pan_vcf
    output:
        chr_vcf = wd + 'temp/{chr}.vcf.gz',
        chr_vcf_tbi = wd + 'temp/{chr}.vcf.gz.tbi'
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        "bcftools view -a -r {params.chr} {input.pan_vcf} | bgzip > {output.chr_vcf} && tabix -f {output.chr_vcf}"

rule chr_vcf_bub:
    input:
        chr_vcf = wd + 'temp/{chr}.vcf.gz'
    output:
        chr_bub_vcf = wd + 'temp/{chr}.bub_10mb.vcf.gz'
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        "vcfbub -i {input.chr_vcf} -l 0 -r 10000000 | bgzip > {output.chr_bub_vcf} && tabix -f {output.chr_bub_vcf}"

ruleorder: pan_vcf_to_chromosomes > chr_vcf_bub

rule alt_alleles_grouping:
    input:
        chr_bub_vcf = wd + 'temp/{chr}.bub_10mb.vcf.gz'
    output:
        chr_small_variants_raw_vcf = wd + 'temp/{chr}.small_variants.raw.vcf.gz',
        chr_SVs_raw_vcf = wd + 'temp/{chr}.SVs.raw.vcf.gz',
    params:
        chr = lambda wildcards: wildcards.chr
    threads: threads
    run:
        headers,vcf_main = vcf_file_read(input.chr_bub_vcf)
        vcf_small_variants, vcf_SVs = multi_process(vcf_main,variant_batch_process,threads)
        with open(wd +'temp/{chr}.small_variants.raw.vcf','w') as f:
            f.writelines(headers[:-1])
            f.write(headers[-1])
            f.writelines(vcf_small_variants)
        with open(wd +'temp/{chr}.SVs.raw.vcf','w') as f:
            f.writelines(headers[:-1])
            f.write(headers[-1])
            f.writelines(vcf_SVs)
        shell('bgzip '+ wd +'temp/{chr}.small_variants.raw.vcf && tabix -f '+ output.chr_small_variants_raw_vcf)    
        shell('bgzip '+ wd +'temp/{chr}.SVs.raw.vcf && tabix -f '+ output.chr_SVs_raw_vcf)    

ruleorder: pan_vcf_to_chromosomes > chr_vcf_bub > alt_alleles_grouping

rule allele_norm:
    input:
        chr_small_variants_raw_vcf = wd + 'temp/{chr}.small_variants.raw.vcf.gz',
        chr_SVs_raw_vcf = wd + 'temp/{chr}.SVs.raw.vcf.gz',
    output:
        chr_small_variants_raw_normed_vcf = wd + 'temp/{chr}.small_variants.raw.normed.vcf.gz',
        chr_SVs_raw_normed_vcf = wd + 'temp/{chr}.SVs.raw.normed.vcf.gz',
        chr_SVs_normed_vcf = wd + 'temp/{chr}.SVs.normed.vcf.gz',
        chr_small_variants_normed_vcf = wd + 'temp/{chr}.small_variants.normed.vcf.gz',
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        bcftools norm -m -any {input.chr_small_variants_raw_vcf} | bgzip > output.chr_small_variants_raw_normed_vcf && tabix -f output.chr_small_variants_raw_normed_vcf
        bcftools norm -m -any {input.chr_SVs_raw_vcf} | bgzip > output.chr_SVs_raw_normed_vcf && tabix -f output.chr_SVs_raw_normed_vcf
        bcftools view -h output.chr_SVs_raw_normed_vcf > {wd}/temp/{params.chr}.header
        cp {wd}/temp/{params.chr}.header {wd}/temp/{params.chr}.SV_loci.small_variants.vcf
        mv {wd}/temp/{params.chr}.header {wd}/temp/{params.chr}.SVs.normed.vcf
        bcftools view -H output.chr_SVs_raw_normed_vcf | \
            awk '{{
                if (length($4) - length($5) >= 50) {{
                    print $0 >> {wd}/temp/{params.chr}.SVs.normed.vcf;
                }} else {{
                    print $0 >> {wd}/temp/{params.chr}.SV_loci.small_variants.vcf;
                }}
            }}'
        bgzip {wd}/temp/{params.chr}.SVs.normed.vcf && tabix -f {output.chr_SVs_normed_vcf}
        bgzip {wd}/temp/{params.chr}.SV_loci.small_variants.vcf && tabix -f {wd}/temp/{params.chr}.SV_loci.small_variants.vcf.gz
        bcftools concat output.chr_small_variants_raw_normed_vcf {wd}/temp/{params.chr}.SV_loci.small_variants.vcf.gz | bgzip > {output.chr_small_variants_normed_vcf} && tabix -f {output.chr_small_variants_normed_vcf}
        """

ruleorder: pan_vcf_to_chromosomes > chr_vcf_bub > alt_alleles_grouping > allele_norm

rule chromosomes_concat:
    input:
        chr_SVs_normed_vcfs = expand(wd + 'temp/{chr}.SVs.normed.vcf.gz',chr = CHROMOSOMES),
        chr_small_variants_normed_vcfs = expand(wd + 'temp/{chr}.small_variants.normed.vcf.gz',chr = CHROMOSOMES),
    output:
        SVs_vcf = wd + out_prefix + '.SVs.normed.vcf.gz',
        SVs_vcf_tbi = wd + out_prefix + '.SVs.normed.vcf.gz.tbi',        
        small_variants_vcf = wd + out_prefix + '.small_variants.normed.vcf.gz',
        small_variants_vcf_tbi = wd + out_prefix + '.small_variants.normed.vcf.gz.tbi',
    shell:
        """
        bcftools concat {input.chr_SVs_normed_vcfs} | bgzip > {output.SVs_vcf} && tabix -f {output.SVs_vcf}
        bcftools concat {input.chr_small_variants_normed_vcfs} | bgzip > {output.small_variants_vcf} && tabix -f {output.small_variants_vcf}
        """