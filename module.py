import gzip
from collections import Counter
from multiprocessing import Pool

def vcf_file_read(filename:str):
    headers = []
    if filename.endswith('.gz'):
        with gzip.open(filename,'r') as f:
            for line in f:
                if line.decode().startswith('#'):
                    headers.append(line.decode())
                else:
                    break
            vcf_main = [i.decode().rstrip('\n').split('\t') for i in f.readlines()]
    else:
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('#'):
                    headers.append(line)
                else:
                    break
            vcf_main = [i.rstrip('\n').split('\t') for i in f.readlines()]
    return headers,vcf_main



def alt_combine(alts):
    alts_sort = sorted(enumerate([len(i) for i in  alts]),key = lambda x:x[1],reverse=True)
    # init alt_dict
    alt_dict = {'.':'.','0':'0'}
    alt_dict.update({str(i+1):str(i+1) for i in range(len(alts))})
    # only group alts longer than 50 bp
    alts_process = [alt for alt in alts_sort if alt[1] > 50]
    if len(alts_process) < 2:
        return ','.join(alts), alt_dict
    else:
        combined = [[alts_process[0][0]]]
        alt_now = 0
        alt_drop = []
        len_now = alts_process[0][1]
        for i in range(1,len(alts_process)):
            if len_now - alts_process[i][1] >= 50:
                alt_now+=1
                combined.append([])
                len_now = alts_process[i][1]
            else:
                alt_drop.append(alts_process[i][0])
            combined[alt_now].append(alts_process[i][0])
            alt_dict[str(alts_process[i][0]+1)] = str(combined[alt_now][0]+1)
        alts_new = [alts[i] for i in range(len(alts)) if i not in alt_drop]
        return ','.join(alts_new), alt_dict


def geno_update(geno,alt_dict):  
    geno = ['|'.join([alt_dict[j] for j in i.split('|')]) for i in geno]
    return geno


def info_cal(geno):
    alleles = sum([i.split('|') for i in geno],[])
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



def variant_process(row):
    geno = row[9:]
    ref = row[3]
    alts = row[4].split(',')
    # whether SV loci
    alleles = [ref]+alts
    len_max = max([len(i) if isinstance(i,str) else 0 for i in alleles])
    if len_max > 50:
        alts_new, alt_dict = alt_combine(alts)
        geno_new = geno_update(geno,alt_dict)
        info_new = info_cal(geno_new)
        row_out = row[:4] + [alts_new,'60', '.',info_new,'GT'] + geno_new
        return True, '\t'.join(row_out)+'\n'
    else:
        return False, '\t'.join(row)+'\n'



def variant_batch_process(vcf_rows):
    grouped_vcf =  [variant_process(i) for i in vcf_rows]
    raw_small_variants = [i[1] for i in grouped_vcf if i[0] == False]
    raw_SVs = [i[1] for i in grouped_vcf if i[0] == True]
    return raw_small_variants, raw_SVs


def multi_process(vcf_main,variant_batch_process,thread_num=16):
    len_sub = len(vcf_main)//thread_num
    list_multi_thread = [[]]
    sub_now = 0
    for i,e in enumerate(vcf_main):
        if i//(len_sub) > sub_now:
            sub_now+=1
            list_multi_thread.append([])
        list_multi_thread[sub_now].append(e)
    p=Pool(thread_num)
    res_dict = {}
    for i in range(thread_num):
        res_dict[i] = p.apply_async(variant_batch_process,(list_multi_thread[i],),)
    p.close()
    p.join()
    out = [res_dict[i].get() for i in range(thread_num)]
    vcf_small_variants = sum([i[0] for i in out],[])
    vcf_SVs = sum([i[1] for i in out],[])
    return vcf_small_variants, vcf_SVs