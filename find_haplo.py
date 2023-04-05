import argparse
import subprocess
import re
import os
from os import path
from datetime import datetime
import numpy as np
import pandas as pd
from Bio import SeqIO
from natsort import natsorted, ns
from itertools import combinations as combi
from functools import reduce


def parse_args():
    parser = argparse.ArgumentParser(description='Find reads contaning SNPs of interest')
    parser.add_argument('-t', '--threads', help='Number of CPU threads to use. Default is 4.', default='4')
    parser.add_argument('-r', '--ref-dir', help='Reference directory which contains all relevant reference files for pipeline.', required=True)
    parser.add_argument('-i','--input', help='Bam file to search for reads containing SNPs of interest', required=True)
    parser.add_argument('-v','--vcf', help='VCF file containing the SNPs of interest', required=True)
    parser.add_argument('-o', '--outputbase', help='Output basename. Default uses fastq name')
    parser.add_argument('-d', '--output-dir', help='Output directory Default is current fastq path/outputbase')
    parser.add_argument('-f','--freq-only', help='Only find frequencies of SNPs, and not look at haplotypes', default=False, action='store_true')

   
    # 
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    args = parser.parse_args()
    return args

def runProcess(cmd, pipe=None):
    print('\n',' '.join(cmd))
    if pipe is None:
        process = subprocess.run(cmd)
    else:
        process = subprocess.run(cmd, stdout=pipe)
    if process.returncode != 0:
        print('\n',cmd,'FAILED\n')
        print()
        quit()

def validateInput(input1, ref_dir):
    if not path.exists(input1):
        print(input1, 'cannot be found')
        quit()
   
    if not path.exists(ref_dir):
        print(ref_dir, 'cannot be found')
        quit()
    else:
        ref_list = ('.fasta', '.dict')
        files_list = os.listdir(ref_dir)
        for ref in ref_list:
            if not any(file.endswith(ref) for file in files_list):
                print('[WARNING] Cannot find reference file ending in', ref)
                quit()
            elif sum([file.endswith(ref) for file in files_list]) >1:
                print('[WARNING] More than one file ending in', ref, 'is found in the provided reference directory')
                quit()


def getAbsPath(input1, ref_dir):
    validateInput(input1,ref_dir)
    input1 = path.abspath(input1)
    ref_dir = path.abspath(ref_dir)
    return input1, ref_dir

def validateThreads(threads):
    cpu_no = os.cpu_count()
    if int(threads) > cpu_no:
        print('\nThere are only', cpu_no, 'threads available.')
        threads = str(cpu_no)
    return threads

def getOutputbaseName(args):
    if args.outputbase is None:
        args.outputbase =  args.input[:-4]
    
    if args.outputbase is not None:
        if args.outputbase.endswith('.bam') or args.outputbase.endswith('.sam'):
            args.outputbase = args.outputbase[:-4]
    return args.outputbase

def getReferenceFasta(ref_dir):
    for file in os.listdir(ref_dir):
        if file.endswith('.fasta'):
            fasta_path = ref_dir+'/'+file
            return fasta_path

def createOutputDir(args):
    if args.output_dir is None:
        file_path = path.abspath(args.input)
        file_dir = path.dirname(file_path)
    elif args.output_dir.endswith('\/'):
        file_dir = args.output_dir[:-1]
    else:
        file_dir = args.output_dir

    if not path.isdir(file_dir):
        os.makedirs(file_dir)
        print('\nOutput directory:', file_dir)
        os.chdir(file_dir)

    return file_dir

def checkArgs(args):
    args.input, args.ref_dir = getAbsPath(args.input, args.ref_dir)
    args.outputbase = getOutputbaseName(args)
    args.output_dir = createOutputDir(args)
    args.threads = validateThreads(args.threads)

    return args

def createConfig(args, print_only):  
    args_tuple = []
    for arg in vars(args):
        argument = arg
        value = getattr(args, arg) 
        if value is None:
            value = 'None'
        args_tuple.append((argument, value))
    
    tab = pd.DataFrame(args_tuple)
    tab.columns = ['Argument', 'Value']
    if print_only:
        print(tab)
    else: 
        output_config = args.outputbase+'_snp.config'
        tab.to_csv(output_config, sep='\t', index=False, header=True)
        return output_config

def Sam2Tsv(args):
    ref_fasta = getReferenceFasta(args.ref_dir)
    tsv_out = args.input[:-4]+'_snp.tsv'
    sub_sam, variant_dict = subSampleHaplo(args)
    cmd = ['sam2tsv', '-r',  ref_fasta, '-o', tsv_out, sub_sam]
    runProcess(cmd)
    return tsv_out, variant_dict

def parseVCF(args):
    df = pd.read_table(args.vcf, usecols=['CHROM', 'POS', 'REF', 'ALT', 'FILTER'])
    pass_df = df[df['FILTER']=='PASS']
    variant_tuple = list(pass_df.itertuples(index=False, name=None))
    return variant_tuple


def haploTSV(args):
    tsv_out, variant_dict = Sam2Tsv(args)

    tsv_file = open(tsv_out, 'r')
    tsv_dict = dict()
    for line in tsv_file:
        line_split = line.strip().split('\t')
        read_name = line_split[0]
        read_base = line_split[4]
        ref_pos = line_split[6]
        ref_base = line_split[7]
        if ref_pos in variant_dict:
            vcf_ref, vcf_alt,_ = variant_dict.get(ref_pos)
            if read_base == vcf_alt and ref_base == vcf_ref:
                if ref_pos not in tsv_dict:
                    tsv_dict[ref_pos] = {read_name}
                else:
                    dict_set = tsv_dict.get(ref_pos)
                    dict_set.update({read_name})
                    tsv_dict[ref_pos] = dict_set

    return tsv_dict, variant_dict

def getVariantHaploDict(variant_tuple, common_set):
    variant_dict = dict()
    for _, pos,ref,alt,_ in variant_tuple:
        variant_dict[str(pos)] = (ref, alt, len(common_set))
    return variant_dict

def getCommonSet(args, variant_tuple):
    common_set = set()
    i=0
    prev_positions = []
    #get read names that are in common for all positions
    for chrom, pos,_,_,_ in variant_tuple:
        i +=1
        pos = str(pos)
        region = chrom+':'+pos+'-'+pos
        names_set = regionNamesSet(args, region)
        repeat_names = common_set.intersection(names_set)
        if i ==1:
            common_set = names_set
        else:
            common_set = repeat_names
        
        if len(common_set) == 0:
            print('\n[WARNING]',pos, 'is not found on the same reads as',', '.join(prev_positions))
            print('Use -f/--freq-only to get their individual frequencies. \nTerminating\n')
            quit()
        prev_positions.append(pos)


    variant_dict = getVariantHaploDict(variant_tuple, common_set)
    return common_set, variant_dict

def getSNPSet(args, variant_tuple):
    merged_set = set()
    #get read names that are in common for all positions
    variant_dict = dict()
    for chrom, pos,ref,alt,_ in variant_tuple:
        pos = str(pos)
        region = chrom+':'+pos+'-'+pos
        names_set = regionNamesSet(args, region)
        repeat_names = merged_set.intersection(names_set)
        uniq_names = names_set-repeat_names
        merged_set.update(uniq_names)
        variant_dict[pos] = (ref, alt, len(names_set))

    return merged_set, variant_dict

def subSampleHaplo(args):
    variant_tuple = parseVCF(args)
    if args.freq_only:
        reads_set, variant_dict = getSNPSet(args, variant_tuple)
    else:
        reads_set, variant_dict = getCommonSet(args, variant_tuple)

    cmd = ['samtools', 'view','-@', args.threads,'-h', args.input]
    samtools_view = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdout, _ = samtools_view.communicate()
    stdout = stdout.decode('utf-8').split('\n') 
    sub_sam = args.outputbase+'.sam'

    fout = open(sub_sam, 'w')
    for line in stdout:
        if line.startswith('@'):
            fout.write(line+'\n')
        elif len(line) > 1:
            line_split = line.split()
            if line_split[0] in reads_set:
                fout.write(line+'\n')
    
    fout.close()

    
    return sub_sam, variant_dict

def regionNamesSet(args, region):
    samtools_view = subprocess.Popen(['samtools', 'view', '-@', args.threads, args.input, region], stdout = subprocess.PIPE)
    stdout, _ = samtools_view.communicate()
    stdout = stdout.decode('utf-8').split('\n')  
    set_out = set()
    for line in stdout:
        line_split = line.split()
        if len(line_split) > 0:
            read_name = line_split[0]
            if read_name not in set_out:
                set_out.add(read_name)
    return set_out

def findFreqs(args, tsv_dict, variant_dict):
    snps = natsorted(list(tsv_dict.keys()))

    snp_tuple = []
    for snp in snps:
        _,_,total_reads = variant_dict.get(snp)
        snp_no = len(tsv_dict.get(snp))
        snp_perc = snp_no/total_reads
        snp_tuple.append((snp, snp_no, snp_perc, total_reads))

    tab = pd.DataFrame(snp_tuple)
    tab.columns = ['SNP','Reads', 'Percent', 'Total_reads']
    print('\n##### SNP FREQUENCIES #####')
    print(tab)
    tab.to_csv(args.outputbase+'_snp-freq.csv', sep='\t', index=False, header=False)


def findHaplo(args, tsv_dict, variant_dict):
    dict_keys = list(tsv_dict.keys())
    _,_,total_reads = variant_dict.get(dict_keys[0])

    results = []
    for r in range(len(tsv_dict)):
        for combo in combi(natsorted(list(tsv_dict.keys())), r+1):
            name = ";".join(combo)
            intersection = reduce(set.intersection, (tsv_dict[n] for n in combo))
            unique = reduce(set.difference, (tsv_dict[n] for n in tsv_dict if n not in combo), intersection)
            uniq_no = len(unique)
            results.append((name, uniq_no, uniq_no/total_reads, intersection))

    tab = pd.DataFrame(results)
    tab.columns = ['SNP','Reads', 'Percent', 'Names']

    tab = tab.drop(labels='Names', axis = 1)
    tab.loc[len(tab)] = ['No_SNPs', total_reads-tab['Reads'].sum(), 1-tab['Percent'].sum()]
    print('\n########## SNP HAPLOTYPES ##########')
    print(tab)
    tab.to_csv(args.outputbase+'_haplotype.csv', sep='\t', index=False, header=False)
    
    for item in results:
        snp_pos, _,_,read_names = item
        fname = snp_pos+'_haplo.names'
        writeHaplos(read_names, fname)

def writeHaplos(region_list, fname):
    original_dir = os.getcwd()
    haplo_dir = './haplotypes'
    if not os.path.exists(haplo_dir):
        os.makedirs(haplo_dir)
    os.chdir(haplo_dir)

    fout = open(fname, 'w')
    for readname in region_list:
        fout.write(readname+'\n')
    fout.close()

    os.chdir(original_dir)

def main():
    args = parse_args()
    args = checkArgs(args)
    createConfig(args, True)    
    tsv_dict, variant_dict = haploTSV(args)
    findFreqs(args, tsv_dict, variant_dict)
    if not args.freq_only:
        findHaplo(args, tsv_dict, variant_dict)

    output_config = createConfig(args, False)


if __name__ == "__main__":
   main()