# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:58:12 2020

@author: lijiang
"""

from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
from itertools import product 
import pandas as pd
import datetime
import gzip
import os
import re

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-f", "--inputFile", help = "input file name", required=True)
parser.add_argument("-m", "--Motif", help = "which format, fa=fasta, m=MEME, s=sequence", type=str, choices=["fas", "mem", "se","gi","all"], required=True)
parser.add_argument("-fm", "--inputmotifFile", help = "input motif file", required=True)
parser.add_argument("-MF", "--MEMEformat", help = "whether to output the MEME format file",type=str, choices=["Y", "N",""], default="N")
parser.add_argument("-CW", "--Combinationw", help = "whether to combination mutations",type=str, choices=["Y", "N",""], default="N")
parser.add_argument("-MT", "--mutationtype", help = "gain or loss",type=str, choices=["Gain", "Loss","All",""],  default="All")
parser.add_argument("-spe", "--species", help = "The name of species", required=True)
parser.add_argument("-out", "--outFile", help = "output ID of file name", required=True)
parser.add_argument("-fan", "--outname", help = "output file name", required=True)
args = parser.parse_args()



inputfile_mutation = args.inputFile
Motif_lag = args.Motif
Motif_file = args.inputmotifFile
MEME_lag = args.MEMEformat
if MEME_lag == "":
    MEME_lag = "N"    
Combination_lag = args.Combinationw
if Combination_lag == "":
    Combination_lag = "N"
mutationtype_lag = args.mutationtype
if mutationtype_lag == "":
    mutationtype_lag = "All"
species = args.species
outnames = args.outFile
out_file_name = args.outname
main_dic = "/var/www/html/Annotation/SBSA/"
os.system("rm %sresult/%s*"%(main_dic,outnames))

###outnames = "test_1"
 
log = open("%sresult/%s.log"%(main_dic,outnames), 'a+')

mapping_name = { "fas": "Multiple Sequences (FASTA)","mem":"Motif (MEME)", "se":"Single Sequence","gi":"Genomic interval file", "all": "Identify all somatic motifs","Y":"Yes","N":"No","All":"Gain and Loss","Gain":"Gain","Loss":"Loss" }   

log.write("Options selected:" +"\n")
log.write("--The input Variant File:   " + inputfile_mutation.split("/")[-1] +"\n")
log.write("--The input Binding Targets:   " + mapping_name[Motif_lag] +"\n")
if len(Motif_file.split("/")) > 1:
    log.write("--More about Binding Targets:   " + Motif_file.split("/")[-1] +"\n")
else:    
    log.write("--More about Binding Targets:   " + Motif_file +"\n")
#log.write("--Whether to output the MEME format file:   " +  mapping_name[MEME_lag] +"\n")
#log.write("--Whether to analysis the combination of multiple mutations:   " +  mapping_name[Combination_lag] +"\n")

if Motif_lag in ["fas", "mem", "se"]:
    if mutationtype_lag != "All":
        log.write("--Somatic Motif effect:   " +  mutationtype_lag +"\n")
    else:
        log.write("--Somatic Motif effect:  Gain and Loss \n")


species_list = pd.read_csv("%sfiles_check/species_list.csv"%main_dic)
species_1 = species_list[species_list["Abb"]  == species]["Name"]
species_1 = ''.join(species_1)
log.write("--The name of species:   " + species_1 +"\n")
log.write("--The name of output file:   " + outnames +"\n\n")




def find_possible_variant_sequences(refSequence, variants):
    resSeq = refSequence.split("$")
    variantsBases = [v[1:] for v in variants]
    res = list(product(*variantsBases))    
    result = []
    for resValues in res:
        for i in range(len(variants)):
            variant = variants[i]
            resSeq[variant[0]] = resValues[i]
        result.append("".join(resSeq))
    return(result,res)


###to find continute location 
def get_blocks(values,dis_value):
    mi, ma = 0, 0
    result = []
    temp = []
    for v in sorted(values):
        if not temp:
            mi = ma = v
            temp.append(v)
        else:
            if abs(v - mi) < dis_value or abs(v - ma) < dis_value:
                temp.append(v)
                if v < mi:
                    mi = v
                elif v > ma:
                    ma = v
            else:
                if len(temp) > 1:
                    result.append(temp)
                mi = ma = v
                temp = [v]
    if len(temp) > 1:
        result.append(temp)    
    return result




###to change the number of chromosome
#species = "aae"
#chr_ID = "supercont1.1"
    
def get_chr(species,chr_ID,main_dic): 
    arr = os.listdir('%sdata_fa/%s'%(main_dic,species))     
    new_strings = []
    for string in arr:
        new_string = string.replace(".fa.gz", "")
        new_strings.append(new_string)    
    
    #print(new_strings)
    if chr_ID=="chrMT":
        chr_ID = "chrM"
   
    if chr_ID in new_strings:
        #print("OK")
        return chr_ID
    else:
        return "Wrong"




###to judge whether the motif is in a sequence
###seq = "AACCTTTTCCAA"        
###ref_seq = "AACTTTTTCCAA" 
##Motifs = ["TTCC","AACCTT","TTCCAA","CCTTTT"]
    
def get_Judge(seq,ref_seq,Motifs):
    seq1 = seq
    ref_seq1 = ref_seq
    seq = seq.upper()
    ref_seq = ref_seq.upper()
    ##print(seq)
    ##print(ref_seq)
    ##print(Motifs)
    re_motif = []
    lags = "delete"
    for motif in Motifs:
        if motif in seq and motif not in ref_seq:
            idx = [i.start() for i in re.finditer(motif,seq)]
            for jj in idx:
                motif1 = seq1[jj:(jj+len(motif))] 
                re_motif.append(["Gain",motif1])
            lags = "Gain"
            
        if motif not in seq and motif in ref_seq:
            idx = [i.start() for i in re.finditer(motif,ref_seq)]
            for jj in idx:
                motif1 = ref_seq1[jj:(jj+len(motif))] 
                re_motif.append(["Loss",motif1])        
            lags = "Loss"      
            
    return lags,re_motif



###to find the complement sequence
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 

def reverse_complement(seq):    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases



def Convert_matrix(sequence):
    sequence = sequence.replace("-","").replace(".","").replace("_","").upper()
    sequence = list(sequence)
    matrix = np.zeros((len(sequence),4))
    base = ["A","C","G","T"]
    num = 0
    for ii in sequence:
        matrix[num,base.index(ii)] = 1
        num = num + 1
    return matrix





def Produce_Motif_sequences(variants):
    variantsBases = [v[1:] for v in variants]
    res = list(product(*variantsBases)) 
    result = []
    for resValues in res:
        result.append("".join(resValues))        
    return(result)


def Produce_variants(random_matrix_array):
    variants = []
    base = ["A","C","G","T"]
    for index in range(np.shape(random_matrix_array)[0]):
        subvariants = [str(index)]
        vectors = random_matrix_array[index,]
        for ii in range(len(vectors)):
            if vectors[ii] == 1:
                subvariants.append(base[ii])
        variants.append(subvariants)
    return(variants)
    



#file_names = "test.txt"    
def Read_MEME_format(file_names,main_dic,outnames):
    meme_o = open("%sresult/%s_meme_fa.txt"%(main_dic,outnames), 'w')
    all_motif = []
    lags1 = 0
    lag2 = 0
    f = open(file_names, 'r')
    sub_matrix = []
    for line in f:
        if line.startswith("MOTIF"):
            motif_id = line.replace("MOTIF","").replace(" ","").rstrip()
    
        if line.startswith("letter"):
            lags1 = 1
            sub_matrix = []
        if not line.startswith("letter") and lags1 == 1 and line not in ['\n', '\r\n'] and not line.startswith("URL"):
            sub_matrix.append(line.split())
            lag2 = 1
            
            
        if lag2 == 1 and lags1 == 1 and (line in ['\n', '\r\n',' \n'] or  line.startswith("URL")):
            lags1 = 0
            lag2 = 0
            sub_matrix = np.asarray(sub_matrix).astype(float)
            sub_matrix[sub_matrix<0.25] = 0
            sub_matrix[sub_matrix>=0.25] = 1            
            seqs_result = Produce_Motif_sequences(Produce_variants(sub_matrix))
            for ii in seqs_result:
                meme_o.write(motif_id +"\t"+ ii +"\n")
                all_motif.append(ii)
                sub_matrix = []
                
    if sub_matrix != []:
        sub_matrix = np.asarray(sub_matrix).astype(float)
        sub_matrix[sub_matrix<0.25] = 0
        sub_matrix[sub_matrix>=0.25] = 1            
        seqs_result = Produce_Motif_sequences(Produce_variants(sub_matrix))
        for ii in seqs_result:
            meme_o.write(motif_id +" "+ ii +"\n")
            all_motif.append(ii)



    meme_o.close()
    return all_motif       
            
           
def Convert_format(result_all,outnames,main_dic,Motif_lag): 
    d = open("%sresult/%s_MEME_format.txt"%(main_dic,outnames), 'w')
    d.write("MEME version 4 \n \nALPHABET= ACGT \n \nstrands: + - \n \nBackground letter frequencies \nA 0.25 C 0.25 G 0.25 T 0.25 \n \n")    
    
    if Motif_lag in ["fas", "mem", "se"]:
        for index, row in result_all.iterrows():
            d.write("MOTIF "+ row['Chromosome'] +":" + str(row['Location']) +":"+row['Mutation']+"||"+ row["Seq.alt"]+"||"+ row['Strand']+"\n" )        
            mat = Convert_matrix(row['Seq.alt'])
            [m,n] = np.shape(mat)
            d.write("letter-probability matrix: alength = %s  w = %s \n"%(n,m))
            np.savetxt(d, mat, fmt='%s',newline = "\n")                   
            d.write( "\n\n\n" ) 
    if Motif_lag == "gi":
        for index, row in result_all.iterrows():
            d.write("MOTIF "+ row['Chromosome'] +":" + str(row['Location']) +":"+row['Mutation']+"||"+ row["New_Seed"]+"\n" )        
            mat = Convert_matrix(row['New_Seed'])
            [m,n] = np.shape(mat)
            d.write("letter-probability matrix: alength = %s  w = %s \n"%(n,m))
            np.savetxt(d, mat, fmt='%s',newline = "\n")                   
            d.write( "\n\n\n" ) 
            
    if Motif_lag == "all":
        for index, row in result_all.iterrows():
            d.write("MOTIF "+ row['Chromosome'] +":" + str(row['Location']) +":"+row['Mutation']+"||"+ row["Seq.alt"]+"||"+ row['Strand']+"\n" )        
            mat = Convert_matrix(row['Seq.alt'])
            [m,n] = np.shape(mat)
            d.write("letter-probability matrix: alength = %s  w = %s \n"%(n,m))
            np.savetxt(d, mat, fmt='%s',newline = "\n")                   
            d.write( "\n\n\n" ) 
                   
                
               
    d.close()
    
    
    

def Class_motif(Motifs):
    Motifs_list = {}
    length_all = []
    for ii in range(len(Motifs)):        
        length_all.append(len(Motifs[ii]))
    for len_ID in np.unique(np.array(length_all)):
        aa = [index for index, value in enumerate(length_all) if (value == len_ID) ] 
        sub_len = [Motifs[i] for i in aa] 
        Motifs_list[len_ID] = sub_len
    return Motifs_list 
    
    


#########################To get the abb of species############################################



def Get_motif_combination(offset,Motifs,result_all,data_chrom,data_position,data_ref,data_A1,data_A2,main_dic,species): 
    offset = offset-1
    for chrom_ID in np.unique(np.array(data_chrom)):
        #print(np.unique(np.array(data_chrom)))
        aa = [index for index, value in enumerate(data_chrom) if (value == chrom_ID) ] 
        sub_position = [data_position[i] for i in aa] 
        sub_ref = [data_ref[i] for i in aa] 
        sub_A1 = [data_A1[i] for i in aa] 
        sub_A2 = [data_A2[i] for i in aa] 
        chrom_ID_lag = get_chr(species,chrom_ID,main_dic)
        
        if chrom_ID_lag ==  "Wrong":
            log.write("%s is not a standard chromosome. Please check the SBSA file all_chromosomes.csv or user-uploaded reference genome file for standard chromosome names presumed for the designated reference genome. \n\n" %(chrom_ID))                
        else:
            chrom_ID = chrom_ID_lag                            
            with gzip.open("%sdata_fa/%s/%s.fa.gz"%(main_dic,species,chrom_ID), "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        seqs = [record]    

            aa = [index for index, value in enumerate(sub_ref) if (len(value) < offset) ] 
            sub_position1 = [sub_position[i] for i in aa]                 
            pos_list = get_blocks(sub_position1,offset)            
            if len(pos_list) != 0:        
                for pos in pos_list: 
                    various_lag = "Yes" 
                    ref_seq = str(seqs[0].seq[(min(pos)-offset-1):(max(pos)+offset+len(sub_ref[sub_position.index(max(pos))])-1)]).lower()
                    
                    for ii in range(len(pos)-1,-1,-1):
                        
                        index_value = sub_position.index(pos[ii])                            
                        ref_p = sub_ref[index_value].replace(" ", "")                           
                        if ii != len(pos)-1:
                            if len(ref_p) <= (pos[ii+1] - pos[ii]):
                                if ref_p != "-":
                                   ref_seq = ref_seq[0:(pos[ii]-min(pos)+offset)].lower() + "$"+ref_p.upper()+ "$"+ref_seq[(pos[ii]-min(pos)+offset)+len(ref_p):].lower()
                                if ref_p == "-":
                                    ref_seq = ref_seq[0:(pos[ii]-min(pos)+offset)].lower() + "$-$"+ref_seq[(pos[ii]-min(pos)+offset):].lower()                                  
                            else:
                                various_lag = "No"
                                break
                        else:
                            if ref_p != "-":
                                ref_seq = ref_seq[0:(pos[ii]-min(pos)+offset)].lower() + "$"+ref_seq[(pos[ii]-min(pos)+offset):(pos[ii]-min(pos)+offset+len(ref_p))].upper()+ "$"+ref_seq[(pos[ii]-min(pos)+offset)+len(ref_p):].lower()  
                            if ref_p == "-": 
                                ref_seq = ref_seq[0:(pos[ii]-min(pos)+offset)].lower() + "$-$"+ref_seq[(pos[ii]-min(pos)+offset):].lower()
                            
                    
                    if various_lag == "Yes":
                        variants = []
                        for ii in range(len(pos)-1,-1,-1):                            
                            index_value = sub_position.index(pos[ii])                            
                            ref_p = sub_ref[index_value].replace(" ", "").upper()
                            list_sub = [(2*ii+1),ref_p]
                                
                            ll = [sub_A1[index_value],sub_A2[index_value]]
                            if list_sub[1] in ll:
                                ll.remove(list_sub[1])
                            list_sub.extend(np.unique(np.array(ll)))
                            variants.append(list_sub)
                            sub_position.pop(index_value)
                            sub_ref.pop(index_value)
                            sub_A1.pop(index_value)
                            sub_A2.pop(index_value)
                            
                        variants.sort()                        
                        results,res = find_possible_variant_sequences(ref_seq, variants)

                        ref_seq_out = results[0]
                        results.remove(results[0])
                        title_o = res[0]
                        res.remove(res[0])
                        
                        
                        for ii in range(len(res)):
                            results[ii] = results[ii].replace("-", "")
                                                                 
                            ref_seq =  ref_seq.replace("$", "").replace("-", "").upper()
                            
                            lags,re_motif =  get_Judge(results[ii],ref_seq,Motifs)
                            if lags != "delete" : 
                                for temp in re_motif:                                  
                                    result_all.append([temp[1],chrom_ID,str(pos),ref_seq_out,results[ii],str(title_o) +">"+ str(res[ii]),temp[0],"+"])
                            
                            
                            lags,re_motif =  get_Judge(reverse_complement(results[ii]),reverse_complement(ref_seq),Motifs)
                            if lags != "delete" :
                                for temp in re_motif:                                    
                                    result_all.append([temp[1],chrom_ID,str(pos),reverse_complement(ref_seq_out),reverse_complement(results[ii]),"("+reverse_complement("-".join(list(title_o))) +")>("+ reverse_complement("-".join(list(res[ii]))) + ")",temp[0],"-"])
           
            
            for nn in range(len(sub_position)):  
                p=int(sub_position[nn])
                ref_p = sub_ref[nn].replace(" ", "").upper()  
                n_ref = len(ref_p)-1
                
                s_L=str(seqs[0].seq[(p-offset-1):(p-1)]).lower()
                if ref_p != "-":
                    s_R=str(seqs[0].seq[(p+n_ref):(p+offset+n_ref)]).lower()
                    ref_seq = (s_L+ref_p+s_R)
                    #ref_seq = str(seqs[0].seq[(p-offset-1):(p+offset+n_ref)]).upper()
                if ref_p == "-":
                    s_R=str(seqs[0].seq[(p-1):(p+offset+n_ref-1)]).lower()
                    ref_seq = (s_L+ref_p+s_R)
                    #ref_seq = str(seqs[0].seq[(p-offset-1):(p+offset+n_ref-1)]).upper()                       
                
                                    
                s = sub_A1[nn].upper()
                if s != ref_p:
                    seq = (s_L+s+s_R)
                    lags,re_motif =  get_Judge(seq.replace("-", ""),ref_seq,Motifs)
                    if lags != "delete" : 
                        for temp in re_motif:
                            result_all.append([temp[1], chrom_ID,str(p),(s_L+ref_p+s_R),seq,ref_p  +">"+ str(s),temp[0],"+"])   
                      
                    lags,re_motif =  get_Judge(reverse_complement(seq.replace("-", "")),reverse_complement(ref_seq),Motifs)
                    if lags != "delete" :   
                        for temp in re_motif:
                            result_all.append([temp[1], chrom_ID,str(p),reverse_complement((s_L+ref_p+s_R)),reverse_complement(seq.replace("-", "")), reverse_complement(ref_p)  +">"+ reverse_complement(str(s)) , temp[0],"-"])
                       
                        
                s = sub_A2[nn].upper()       
                if s != ref_p and s != sub_A1[nn].upper() :            
                    seq = (s_L+s+s_R)
                    lags,re_motif =  get_Judge(seq.replace("-", ""),ref_seq,Motifs)
                    if lags != "delete" :
                        for temp in re_motif:
                            result_all.append([temp[1],chrom_ID,str(p),(s_L+ref_p+s_R),seq,ref_p  +">"+ str(s),temp[0],"+"])
                        
                    lags,re_motif =  get_Judge(reverse_complement(seq.replace("-", "")),reverse_complement(ref_seq),Motifs)
                    if lags != "delete" : 
                        for temp in re_motif:
                            result_all.append([temp[1],chrom_ID,str(p),reverse_complement((s_L+ref_p+s_R)),reverse_complement(seq), reverse_complement(ref_p)  +">"+ reverse_complement(str(s)), temp[0],"-"])   
    #print(result_all)
    return result_all
    
               

def Get_motif_single(offset,Motifs,result_all,data_chrom,data_position,data_ref,data_A1,data_A2,main_dic,species): 
    offset = offset-1
    for chrom_ID in np.unique(np.array(data_chrom)):
        #print(chrom_ID)
        aa = [index for index, value in enumerate(data_chrom) if (value == chrom_ID) ] 
        sub_position = [data_position[i] for i in aa] 
        sub_ref = [data_ref[i] for i in aa] 
        sub_A1 = [data_A1[i] for i in aa] 
        sub_A2 = [data_A2[i] for i in aa] 
        chrom_ID_lag = get_chr(species,chrom_ID,main_dic)
        
        if chrom_ID_lag ==  "Wrong":
            log.write("%s is not a standard chromosome. Please check the SBSA file all_chromosomes.csv or user-uploaded reference genome file for standard chromosome names presumed for the designated reference genome.\n\n" %(chrom_ID))                
        else:
            chrom_ID = chrom_ID_lag                            
            with gzip.open("%sdata_fa/%s/%s.fa.gz"%(main_dic,species,chrom_ID), "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        seqs = [record]    
           
            for nn in range(len(sub_position)):  
                p=int(sub_position[nn])
                ref_p = sub_ref[nn].replace(" ", "").upper()  
                n_ref = len(ref_p)-1
                
                s_L=str(seqs[0].seq[(p-offset-1):(p-1)]).lower()
                if ref_p != "-":
                    s_R=str(seqs[0].seq[(p+n_ref):(p+offset+n_ref)]).lower()
                    #ref_seq = str(seqs[0].seq[(p-offset-1):(p+offset+n_ref)]).upper()
                    ref_seq = (s_L+ref_p+s_R)
                if ref_p == "-":
                    s_R=str(seqs[0].seq[(p-1):(p+offset+n_ref-1)]).lower()
                    #ref_seq = str(seqs[0].seq[(p-offset-1):(p+offset+n_ref-1)]).upper()  
                    ref_seq = (s_L+s_R)                    
                                   
                s = sub_A1[nn].upper()
                if s != ref_p:
                    seq = (s_L+s+s_R)
                    lags,re_motif =  get_Judge(seq.replace("-", ""),ref_seq,Motifs)
                    if lags != "delete" : 
                        for temp in re_motif:
                            result_all.append([temp[1], chrom_ID,str(p),(s_L+ref_p+s_R),seq,ref_p  +">"+ str(s),temp[0],"+"])   
                      
                    lags,re_motif =  get_Judge(reverse_complement(seq.replace("-", "")),reverse_complement(ref_seq),Motifs)
                    if lags != "delete" :   
                        for temp in re_motif:
                            result_all.append([temp[1], chrom_ID,str(p),reverse_complement((s_L+ref_p+s_R)),reverse_complement(seq.replace("-", "")), reverse_complement(ref_p)  +">"+ reverse_complement(str(s)) , temp[0],"-"])
                       
                        
                s = sub_A2[nn].upper()       
                if s != ref_p and s != sub_A1[nn].upper() :            
                    seq = (s_L+s+s_R)
                    lags,re_motif =  get_Judge(seq.replace("-", ""),ref_seq,Motifs)
                    if lags != "delete" :
                        for temp in re_motif:
                            result_all.append([temp[1],chrom_ID,str(p),(s_L+ref_p+s_R),seq,ref_p  +">"+ str(s),temp[0],"+"])
                        
                    lags,re_motif =  get_Judge(reverse_complement(seq.replace("-", "")),reverse_complement(ref_seq),Motifs)
                    if lags != "delete" : 
                        for temp in re_motif:
                            result_all.append([temp[1],chrom_ID,str(p),reverse_complement((s_L+ref_p+s_R)),reverse_complement(seq), reverse_complement(ref_p)  +">"+ reverse_complement(str(s)), temp[0],"-"])   
    
    return result_all                


##inputfile_mutation = r"D:\Jiang_work\Projects\MOTIF\supplement_experiment\Example.csv"

####to read the input file
data_chrom = []
data_position = []
data_ref = []
data_A1 = []
data_A2 = []


f = open(inputfile_mutation, 'r')
words = f.readline().split(",") 
if words[1].isdigit() and "#" not in words:
    words[-1] = words[-1].strip() 
    words[0] = words[0].replace('"', "")    
    data_chrom.append(words[0])
    data_position.append(int(words[1]))
    data_ref.append(words[2].upper().replace(".", "-"))
    data_A1.append(words[3].upper().replace(".", "-")  )
    data_A2.append(words[4].upper().replace(".", "-"))

for line in f:
    if line not in [' \n','\n','\t\n','\r\n']:
        line = line.replace("\"", "").replace("\'", "")
        words = line.split(",")    
        words[-1] = words[-1].strip() 
        words[0] = words[0].replace('"', "")        
        data_chrom.append(words[0])        
        data_position.append(int(words[1]))
        data_ref.append(words[2].upper().replace(".", "-"))
        data_A1.append(words[3].upper().replace(".", "-")  )
        data_A2.append(words[4].upper().replace(".", "-"))       
        
f.close()


#Motif_lag = "fas"
#####################to find motifs###########################################

if Motif_lag in ["fas", "mem", "se"]:
    if Motif_lag == "fas":
        Motifs = []
        # Motif_file =  "seq_unique.fasta"
        seqs = list(SeqIO.parse(Motif_file, "fasta")) 
        for ii in seqs:
            Motifs.append(str(ii.seq))
        #print(Motifs)
        Motifs = np.unique(np.array(Motifs))   

        
    if Motif_lag == "mem":
        Motifs = Read_MEME_format(Motif_file,main_dic,outnames)  
        #print(Motifs)
        Motifs = np.unique(np.array(Motifs)) 
        
        
        
    if Motif_lag == "se":    
        Motifs = [Motif_file.upper().replace(" ", "")]
    

    Motifs = Class_motif(Motifs)
    
    result_all = [] 
    for key_word in Motifs:
        print(key_word)
        sub_Motif = Motifs[key_word]
        log.write("--User supplied %s binding sequences of length %s nt.\n\n"%(len(sub_Motif), key_word))
        print(len(sub_Motif))
        if Combination_lag == "Y":
            result_all = Get_motif_combination(key_word,sub_Motif,result_all,data_chrom,data_position,data_ref,data_A1,data_A2,main_dic,species)
        
        if Combination_lag == "N":
            result_all = Get_motif_single(key_word,sub_Motif,result_all,data_chrom,data_position,data_ref,data_A1,data_A2,main_dic,species)
        
        
    if not result_all == []: 
        
        df = pd.DataFrame(result_all)
        title = ["Bind_seq","Chromosome","Location","Seq.ref","Seq.alt","Mutation","Types","Strand"]
        
        if mutationtype_lag == "Loss":
            aa = (df[5] == "Loss")
            df1 = df[aa]   
            
        if mutationtype_lag == "Gain": 
            aa = (df[5] == "Gain")
            df1 = df[aa]
            
        if mutationtype_lag == "All": 
            df1 = df
        
        if len(df1) > 0:
            log.write("Some somatic motifs were found and saved in %s%s."%(outnames,out_file_name) +"\n\n")  
            df1 = pd.DataFrame({title[0]: df1[0], title[1]: df1[1],title[2]: df1[2],title[3]: df1[3],title[4]: df1[4],title[5]: df1[5],title[6]: df1[6],title[7]: df1[7]})
            column_order = ['Chromosome','Location','Mutation','Strand',"Seq.ref","Seq.alt",'Types','Bind_seq']           
            df1[column_order].to_csv('%sresult/%s%s'%(main_dic,outnames,out_file_name), header=True, index=False)
            if MEME_lag == "Y":
                Convert_format(df1,outnames,main_dic,Motif_lag)

        else:
            log.write("No somatic motifs were found." +"\n\n")             
                
    else:
        log.write("No somatic motifs were found." +"\n\n")                        






#####################to find miRNA###########################################

##Motif_file = "hsa.bed.csv"
if Motif_lag == "gi":

    d_miR = pd.read_csv("%s"%(Motif_file))
    title_original = d_miR.columns.values
    title_original = np.array(title_original)
    title_original[0:4] = ['V1','V2','V3','V4']
    d_miR.columns.values[0:4] = ['V1','V2','V3','V4']
    miRNA_result = pd.DataFrame()
    
    
    for chrom_ID in np.unique(np.array(data_chrom)):
        #print(chrom_ID)
        aa = [index for index, value in enumerate(data_chrom) if (value == chrom_ID) ] 
        sub_position = [data_position[i] for i in aa] 
        sub_ref = [data_ref[i] for i in aa] 
        sub_A1 = [data_A1[i] for i in aa] 
        sub_A2 = [data_A2[i] for i in aa] 
        chrom_ID_lag = get_chr(species,chrom_ID,main_dic)
        
        if chrom_ID_lag ==  "Wrong":
            log.write("%s is not a standard chromosome. Please check the SBSA file all_chromosomes.csv or user-uploaded reference genome file for standard chromosome names presumed for the designated reference genome.\n\n" %(chrom_ID))                
        else:
            chrom_ID = chrom_ID_lag
            with gzip.open("%sdata_fa/%s/%s.fa.gz"%(main_dic,species,chrom_ID), "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        seqs = [record]
           

            sub_d = d_miR[d_miR["V1"] == chrom_ID]
            for nn in range(len(sub_position)):
                pos = int(sub_position[nn])   
                set1 = sub_d["V2"] <= pos
                set1 = [i for i, x in enumerate(set1) if x]
                set2 = sub_d["V3"] >= pos
                set2 = [i for i, x in enumerate(set2) if x]
                index = list(set(set1) & set(set2))
                if len(index) != 0:  
                    aa_all = sub_d.iloc[index]
                    for idex_aa in range(aa_all.shape[0]):
                        title_o_t = title_original
                        aa = aa_all.iloc[idex_aa]
                        length_beed = int(aa["V3"]) - int(aa["V2"])+1
                        aa['Chromosome'] =  chrom_ID       
                        title_o_t = np.append(title_o_t,'Chromosome') 
                        aa['Location'] = pos
                        title_o_t = np.append(title_o_t,'Location')                        
                        ref_p = sub_ref[nn].replace(" ", "").upper() 
                        ref_n = len(ref_p) - 1             
                    
                           
                        s_L=str(seqs[0].seq[int(aa["V2"])-1:(pos-1)]).lower()                        
                        if ref_p != "-":
                            s_R=str(seqs[0].seq[(pos+ref_n):(int(aa["V3"])+ref_n)]).lower()
                        if ref_p == "-":
                            s_R=str(seqs[0].seq[(pos-1):int(aa["V3"])]).lower()  
                        #seq = str(seqs[0].upper().seq[int(aa["V2"]) -1:int(aa["V3"])])  
                        seq = str(s_L+ref_p+s_R)                         
                        
                        
                        if ''.join(aa["V4"]) == "-":
                            seq = reverse_complement(seq)  
                        aa['Original_Seed'] = seq  
                        title_o_t = np.append(title_o_t,'Original_Seed')     
        
                        s = sub_A1[nn].upper().replace(" ", "")
                        if s != ref_p:
                            seq =str(s_L+s+s_R) 
                            #seq1 = seq
                            #if ''.join(aa["V4"]) == "-":
                            #   seq1 = reverse_complement(seq1)                        
                            #aa['Mutation_Sequence'] = seq1  
                            #title_o_t = np.append(title_o_t,'Mutation_Sequence')    
                            seq = seq.replace("-", "")
                            if len(seq)>length_beed:
                                seq = seq[0:length_beed]
                            if len(seq)<length_beed:
                                temp_len = length_beed-len(seq)
                                if ref_p != "-":
                                    seq = seq + str(seqs[0].upper().seq[(int(aa["V3"])+ref_n):(int(aa["V3"])+ref_n+temp_len)])
                                if ref_p == "-":
                                    seq = seq + str(seqs[0].upper().seq[(int(aa["V3"])):(int(aa["V3"])+temp_len)])
                            
                            if ''.join(aa["V4"]) == "-":
                                seq = reverse_complement(seq) 
                            aa['New_Seed'] = seq
                            title_o_t = np.append(title_o_t,'New_Seed')    
                            if ''.join(aa["V4"]) == "-":
                                aa['Mutation'] = reverse_complement(ref_p)  +">"+ reverse_complement(s)
                            else:
                                aa['Mutation'] = ref_p  +">"+ s
                            title_o_t = np.append(title_o_t,'Mutation')    
                            aa.index =  title_o_t 
                            
                            miRNA_result = miRNA_result.append(aa)
                        
                        
                        s = sub_A2[nn].upper().replace(" ", "")
                        if s != ref_p and s != sub_A1[nn].upper():
                            seq =str(s_L+s+s_R)    
                            #seq1 = seq
                            #if ''.join(aa["V4"]) == "-":
                            #    seq1 = reverse_complement(seq1)                        
                            #aa['Mutation_Sequence'] = seq1 
                            #title_o_t = np.append(title_o_t,'Mutation_Sequence') 
                            seq = seq.replace("-", "")
                            if len(seq)>length_beed:
                                seq = seq[0:length_beed]
                            if len(seq)<length_beed:
                                temp_len = length_beed-len(seq)
                                if ref_p != "-":
                                    seq = seq + str(seqs[0].upper().seq[(int(aa["V3"])+ref_n):(int(aa["V3"])+ref_n+temp_len)])
                                if ref_p == "-":
                                    seq = seq + str(seqs[0].upper().seq[(int(aa["V3"])):(int(aa["V3"])+temp_len)])
                            
                            if ''.join(aa["V4"]) == "-":
                                seq = reverse_complement(seq) 
                            aa['New_Seed'] = seq
                            if 'New_Seed' not in title_o_t:
                                title_o_t = np.append(title_o_t,'New_Seed')  
                            if ''.join(aa["V4"]) == "-":
                                aa['Mutation'] = reverse_complement(ref_p)  +">"+ reverse_complement(s)
                            else:
                                aa['Mutation'] = ref_p  +">"+ s
                            if 'Mutation' not in title_o_t:
                                title_o_t = np.append(title_o_t,'Mutation')
                            aa.index =  title_o_t 
                            miRNA_result = miRNA_result.append(aa)
    
    if miRNA_result.shape[0] > 0:
        log.write("Some affected miRNA binding sequences were found and saved in %s%s"%(outnames,out_file_name)+".\n\n")        
        miRNA_result[title_o_t].to_csv('%sresult/%s%s'%(main_dic,outnames,out_file_name), header=True, index=False)
        
        
        if MEME_lag == "Y":
            Convert_format(miRNA_result,outnames,main_dic,Motif_lag)
        
        
    else:
        log.write("No affected miRNA binding sequences were found." +"\n\n")
        
        
#####################to extract fragments for all mutation###########################################

if Motif_lag == "all":
    result_all = []
    length = int(Motif_file)
    for chrom_ID in np.unique(np.array(data_chrom)):
        #print(chrom_ID)
        aa = [index for index, value in enumerate(data_chrom) if (value == chrom_ID) ] 
        sub_position = [data_position[i] for i in aa] 
        sub_ref = [data_ref[i] for i in aa] 
        sub_A1 = [data_A1[i] for i in aa] 
        sub_A2 = [data_A2[i] for i in aa] 
        chrom_ID_lag = get_chr(species,chrom_ID,main_dic)
        
        
        if chrom_ID_lag ==  "Wrong":
            log.write("%s is not a standard chromosome. Please check the SBSA file all_chromosomes.csv or user-uploaded reference genome file for standard chromosome names presumed for the designated reference genome.\n\n" %(chrom_ID))                
        else:
            chrom_ID = chrom_ID_lag
            with gzip.open("%sdata_fa/%s/%s.fa.gz"%(main_dic,species,chrom_ID), "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        seqs = [record]
                        #print(len(seqs[0].seq))
            
            for nn in range(len(sub_position)):  
                p=int(sub_position[nn])
                ref_p = sub_ref[nn].replace(" ", "").upper() 
                ref_n = len(ref_p) - 1
                s_L=str(seqs[0].seq[(p-length-1):(p-1)]).lower()
                if ref_p != "-":
                    s_R=str(seqs[0].seq[(p+ref_n):(p+length+ref_n)]).lower()  
                if ref_p == "-":
                    s_R=str(seqs[0].seq[(p-1):(p+length-1)]).lower()     
                
                s = sub_A1[nn].upper()
                if s != ref_p:
                    seq = (s_L+s+s_R)
                    seq1 = (s_L+ref_p+s_R)                    
                    result_all.append([chrom_ID,str(p),seq, seq1,ref_p +">"+ str(s),"+"])
                s = sub_A2[nn].upper()       
                if s != ref_p and s != sub_A1[nn].upper():            
                    seq = (s_L+s+s_R)
                    seq1 = (s_L+ref_p+s_R) 
                    result_all.append([chrom_ID,str(p),seq, seq1,ref_p  +">"+ str(s),"+"])
                
     
    
    if len(result_all) > 0:
        log.write("Some somatic sequences were extracted and saved in %s%s"%(outnames,out_file_name) +".\n\n")
        df = pd.DataFrame(result_all)
        title = ["Chromosome","Location","Seq.alt","Seq.ref","Mutation","Strand"]
        df = pd.DataFrame({title[0]: df[0], title[1]: df[1],title[2]: df[2],title[3]: df[3],title[4]: df[4],title[5]: df[5]})
        title = ["Chromosome","Location","Mutation","Strand","Seq.ref","Seq.alt"]
        df[title].to_csv('%sresult/%s%s'%(main_dic,outnames,out_file_name), header=True, index=False)
        if MEME_lag == "Y":
            Convert_format(df,outnames,main_dic,Motif_lag)
    else:
        log.write("No somatic sequences were extracted." +"\n\n")



#log.write(str(datetime.datetime.now()))
log.close()


# os.system("tar -C %s  -zcvf  %s.tar.gz result/%s*"%(main_dic, outnames,outnames))
# c_dir = os.getcwd()
# os.system("rm %sresult/%s*"%(main_dic,outnames))
# os.system("mv %s/%s.tar.gz %sresult/%s.tar.gz"%(c_dir,outnames,main_dic,outnames))
