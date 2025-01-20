# import grouping fuctions
from module import *

# Load configuration from config.yaml
configfile: "config.yaml"
pan_vcf = config['pan_vcf']
out_prefix = config['out_prefix']
wd = config['work_dir'] + '/'
threads = int(config['threads'])

# Get the reference chromosomes
CHROMOSOMES = ['CHM13v2.chr'+i for i in list(map(str,range(1,23)))+["X","Y","M"]]

rule all:
    input:
        SVs_vcf = wd + out_prefix + '.SVs.normed.vcf.gz',
        SVs_vcf_tbi = wd + out_prefix + '.SVs.normed.vcf.gz.tbi',        
        small_variants_vcf = wd + out_prefix + '.small_variants.normed.vcf.gz',
        small_variants_vcf_tbi = wd + out_prefix + '.small_variants.normed.vcf.gz.tbi',
        SVs_map = expand(wd + 'temp/{chr}/temp.SVs.merge_map.csv',chr = CHROMOSOMES)


# Remove too large variants
rule vcfbub:
    input:
        pan_vcf = pan_vcf,
    output:
        bub_vcf = wd + 'temp/temp.bub_10mb.vcf.gz',
        bub_vcf_tbi = wd + 'temp/temp.bub_10mb.vcf.gz.tbi',
    shell:
        "vcfbub -i {input.pan_vcf} -l 0 -r 10000000 | bgzip -f > {output.bub_vcf} && tabix -f {output.bub_vcf}"

# Split to chromosomes
rule pan_vcf_to_chromosomes:
    input:
        bub_vcf = wd + 'temp/temp.bub_10mb.vcf.gz',
        bub_vcf_tbi = wd + 'temp/temp.bub_10mb.vcf.gz.tbi',
    output:
        chr_vcf = wd + 'temp/{chr}/temp.vcf.gz',
        chr_vcf_tbi = wd + 'temp/{chr}/temp.vcf.gz.tbi'
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        "bcftools view -a -r {params.chr} {input.bub_vcf} -Oz -o {output.chr_vcf} && tabix -f {output.chr_vcf}"

ruleorder: vcfbub > pan_vcf_to_chromosomes

# Grouping SV alleles
rule alt_alleles_grouping:
    input:
        chr_vcf = wd + 'temp/{chr}/temp.vcf.gz',
    output:
        chr_small_variants_raw_vcf = wd + 'temp/{chr}/temp.small_variants.raw.vcf.gz',
        chr_SVs_raw_vcf = wd + 'temp/{chr}/temp.SVs.raw.vcf.gz',
    params:
        chr = lambda wildcards: wildcards.chr
    threads: threads
    run:
        headers,vcf_main = vcf_file_read(input.chr_vcf)
        vcf_small_variants, vcf_SVs df_map = multi_process(vcf_main,variant_batch_process,threads)
        with open(wd +'temp/'+params.chr+'/temp.small_variants.raw.vcf','w') as f:
            f.writelines(headers)
            f.writelines(vcf_small_variants)
        with open(wd +'temp/'+params.chr+'/temp.SVs.raw.vcf','w') as f:
            f.writelines(headers)
            f.writelines(vcf_SVs)
        df_map.to_csv(wd +'temp/'+params.chr+'/temp.SVs.merge_map.csv',index=False)
        shell('bgzip -f '+ wd +'temp/'+params.chr+'/temp.small_variants.raw.vcf && tabix -f '+ output.chr_small_variants_raw_vcf)    
        shell('bgzip -f '+ wd +'temp/'+params.chr+'/temp.SVs.raw.vcf && tabix -f '+ output.chr_SVs_raw_vcf)    

ruleorder: vcfbub > pan_vcf_to_chromosomes > alt_alleles_grouping

# Split to biallelics and classify to SVs and small variants
rule allele_norm:
    input:
        chr_small_variants_raw_vcf = wd + 'temp/{chr}/temp.small_variants.raw.vcf.gz',
        chr_SVs_raw_vcf = wd + 'temp/{chr}/temp.SVs.raw.vcf.gz',
    output:
        chr_small_variants_raw_normed_vcf = wd + 'temp/{chr}/temp.small_variants.raw.normed.vcf.gz',
        chr_SVs_raw_normed_vcf = wd + 'temp/{chr}/temp.SVs.raw.normed.vcf.gz',
        chr_SVs_normed_vcf = wd + 'temp/{chr}/temp.SVs.normed.vcf.gz',
        chr_small_variants_normed_vcf = wd + 'temp/{chr}/temp.small_variants.normed.vcf.gz',
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        bcftools norm -m -any {input.chr_small_variants_raw_vcf} -Oz -o {wd}/temp/{params.chr}/temp.small_variants.raw.normed.vcf.gz && tabix -f {wd}/temp/{params.chr}/temp.small_variants.raw.normed.vcf.gz
        bcftools norm -m -any {input.chr_SVs_raw_vcf} -Oz -o {wd}/temp/{params.chr}/temp.SVs.raw.normed.vcf.gz && tabix -f {wd}/temp/{params.chr}/temp.SVs.raw.normed.vcf.gz
        bcftools view -h {wd}/temp/{params.chr}/temp.SVs.raw.normed.vcf.gz > {wd}/temp/{params.chr}/temp.header
        cp {wd}/temp/{params.chr}/temp.header {wd}/temp/{params.chr}/temp.SV_loci.small_variants.vcf
        mv {wd}/temp/{params.chr}/temp.header {wd}/temp/{params.chr}/temp.SVs.normed.vcf
        bcftools view -H {wd}/temp/{params.chr}/temp.SVs.raw.normed.vcf.gz | \
            awk '{{
                if (length($4) - length($5) >= 50 || length($5) - length($4) >= 50) {{
                    print $0 >> "{wd}/temp/{params.chr}/temp.SVs.normed.vcf";
                }} else {{
                    print $0 >> "{wd}/temp/{params.chr}/temp.SV_loci.small_variants.vcf";
                }}
            }}'
        bgzip -f {wd}/temp/{params.chr}/temp.SVs.normed.vcf && tabix -f {output.chr_SVs_normed_vcf}
        bgzip -f {wd}/temp/{params.chr}/temp.SV_loci.small_variants.vcf && tabix -f {wd}/temp/{params.chr}/temp.SV_loci.small_variants.vcf.gz
        bcftools concat {wd}/temp/{params.chr}/temp.small_variants.raw.normed.vcf.gz {wd}/temp/{params.chr}/temp.SV_loci.small_variants.vcf.gz | bcftools sort -Oz -o {output.chr_small_variants_normed_vcf} && tabix -f {output.chr_small_variants_normed_vcf}
        """

ruleorder: vcfbub > pan_vcf_to_chromosomes > alt_alleles_grouping > allele_norm

# Concat the chromosomes
rule chromosomes_concat:
    input:
        chr_SVs_normed_vcfs = expand(wd + 'temp/{chr}/temp.SVs.normed.vcf.gz',chr = CHROMOSOMES),
        chr_small_variants_normed_vcfs = expand(wd + 'temp/{chr}/temp.small_variants.normed.vcf.gz',chr = CHROMOSOMES),
    output:
        SVs_vcf = wd + out_prefix + '.SVs.normed.vcf.gz',
        SVs_vcf_tbi = wd + out_prefix + '.SVs.normed.vcf.gz.tbi',        
        small_variants_vcf = wd + out_prefix + '.small_variants.normed.vcf.gz',
        small_variants_vcf_tbi = wd + out_prefix + '.small_variants.normed.vcf.gz.tbi',
    run:
        SVs_vcfs = [wd + 'temp/' + chrom + '/temp.SVs.normed.vcf.gz' for chrom in CHROMOSOMES]
        small_variants_vcfs = [wd + 'temp/' + chrom + '/temp.small_variants.normed.vcf.gz' for chrom in CHROMOSOMES]
        cmd = 'bcftools concat ' + ' '.join(SVs_vcfs) + ' -Oz -o {output.SVs_vcf} && tabix -f {output.SVs_vcf}'
        shell(cmd)
        cmd = 'bcftools concat ' + ' '.join(small_variants_vcfs) + ' -Oz -o {output.small_variants_vcf} && tabix -f {output.small_variants_vcf}'
        shell(cmd)
