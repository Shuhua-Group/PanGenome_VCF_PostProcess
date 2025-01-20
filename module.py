import gzip
from collections import Counter
from multiprocessing import Pool
import sys,re
import pandas as pd
import argparse


def vcf_file_read(filename:str):
    headers = []
    vcf_main = []
    if filename.endswith('.gz'):
        with gzip.open(filename,'r') as f:
            for line in f:
                if line.decode().startswith('#'):
                    headers.append(line.decode())
                else:
                    vcf_main.append(line.decode().rstrip('\n').split('\t'))
                    break
            vcf_main+=[i.decode().rstrip('\n').split('\t') for i in f.readlines()]
    else:
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('#'):
                    headers.append(line)
                else:
                    vcf_main.append(line.rstrip('\n').split('\t'))
                    break
            vcf_main = [i.rstrip('\n').split('\t') for i in f.readlines()]
    return headers,vcf_main



def alt_combine(alts, ref):
    alts_sort = sorted(enumerate([len(i) for i in  alts]),key = lambda x:x[1],reverse=True)
    # init alt_dict
    alt_dict = {'.':'.',0:0}
    alt_dict.update({i+1:i+1 for i in range(len(alts))})
    # only group alts longer than 50 bp
    #alts_sort = [alt for alt in alts_sort if alt[1] > 50]
    if len(alts_sort) < 2:
        return ','.join(alts), {str(key): str(value) for key, value in alt_dict.items()}
    else:
        combined = [[alts_sort[0][0]]]
        alt_now = 0
        alt_drop = []
        len_now = alts_sort[0][1]
        for i in range(1,len(alts_sort)):
            if len_now - alts_sort[i][1] >= 50:
                alt_now+=1
                combined.append([alts_sort[i][0]])
                len_now = alts_sort[i][1]
            else:
                alt_drop.append(alts_sort[i][0])
                combined[alt_now].append(alts_sort[i][0])
        for i in range(len(combined)):
            alt_dict.update({alt_now+1:i+1 for alt_now in combined[i]})
        alts_new = [alts[i[0]] for i in combined]

        alts_out = []
        alt_dict_out = {'.':'.',0:0}
        alt_dict.pop('.')
        for i, alt in enumerate(alts_new):
            if abs(len(alt) - len(ref)) >= 50:
                alts_out.append(alt)
                alt_dict_out.update({key:value for key, value in alt_dict.items() if value == i+1})
            else:
                alt_dict_out.update({key:0 for key, value in alt_dict.items() if value == i+1})
                for j in range(i+1, len(alts_new)):
                    alt_dict_out.update({key:value-1 for key, value in alt_dict.items() if value > i})
                break
        return ','.join(alts_out), {str(key): str(value) for key, value in alt_dict_out.items()}





def geno_update(geno,alt_dict):  
    geno = ['|'.join([alt_dict[j] for j in re.split(r'[|/]', i)]) for i in geno]
    return geno


def info_cal(geno):
    alleles = sum([re.split(r'[|/]', i) for i in geno],[])
    allele_count = Counter(alleles)
    alts = set(allele_count.keys()) - {'0','.'}
    alleles_avail = [i for i in alleles if i!='.']
    AN = len(alleles_avail)
    AC = []
    AF = []
    for i in range(1,len(alts)+1):
        AC.append(str(allele_count[str(i)]))
        AF.append(str(round(allele_count[str(i)]/AN,6)))
    INFO = ['AC='+','.join(AC),'AF='+','.join(AF),'AN='+str(AN),'NS='+str(len([i for i in geno if i not in {'.|.','.'}]))]
    return ';'.join(INFO)

def merge_map(row, alt_dict):
    alts = row[4].split(',')
    merge_map = []
    for i in range(len(alts)):
        merge_map.append(row[:2]+[len(alts[i]),i+1,alt_dict[str(i+1)], False])
    merge_map.append(row[:2]+[len(row[3]),0,0, True])
    df_merge = pd.DataFrame(merge_map,columns=['CHROM','POS','Len','Before','After','is_REF'])
    return df_merge


def variant_process(row):
    geno = row[9:]
    ref = row[3]
    alts = row[4].split(',')
    # whether SV loci
    alleles = [ref]+alts
    len_max = max([len(i)-len(ref) if isinstance(i,str) else 0 for i in alleles])
    if len_max > 50:
        alts_new, alt_dict = alt_combine(alts,ref)
        geno_new = geno_update(geno,alt_dict)
        info_new = info_cal(geno_new)
        row_out = row[:4] + [alts_new,'60', '.',info_new,'GT'] + geno_new
        df_map = merge_map(row, alt_dict)
        return True, '\t'.join(row_out)+'\n', df_map
    else:
        return False, '\t'.join(row)+'\n'



def variant_batch_process(vcf_rows):
    grouped_vcf =  [variant_process(i) for i in vcf_rows]
    raw_small_variants = [i[1] for i in grouped_vcf if i[0] == False]
    raw_SVs = [i[1] for i in grouped_vcf if i[0] == True]
    df_map = pd.concat([i[2] for i in grouped_vcf if i[0] == True])
    return raw_small_variants, raw_SVs, df_map


def multi_process(vcf_main,threads=16):
    len_sub = len(vcf_main)//threads
    list_multi_thread = [[]]
    sub_now = 0
    for i,e in enumerate(vcf_main):
        if i//(len_sub) > sub_now:
            sub_now+=1
            list_multi_thread.append([])
        list_multi_thread[sub_now].append(e)
    p=Pool(threads)
    res_dict = {}
    for i in range(threads):
        res_dict[i] = p.apply_async(variant_batch_process,(list_multi_thread[i],),)
    p.close()
    p.join()
    out = [res_dict[i].get() for i in range(threads)]
    vcf_small_variants = sum([i[0] for i in out],[])
    vcf_SVs = sum([i[1] for i in out],[])
    df_map = pd.concat([i[2] for i in out])
    return vcf_small_variants, vcf_SVs, df_map

def main(input_vcf, out_prefix, threads=16):
    headers, vcf_main = vcf_file_read(input_vcf)
    vcf_small_variants, vcf_SVs, df_map = multi_process(vcf_main, threads)
    with open(f'{out_prefix}.small.vcf', 'w') as f:
        f.writelines(headers+vcf_small_variants)
    with open(f'{out_prefix}.SV.vcf', 'w') as f:
        f.writelines(headers+vcf_SVs)
    df_map.to_csv(f'{out_prefix}.map.csv', index=False)
