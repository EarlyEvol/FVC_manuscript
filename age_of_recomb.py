# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 14:50:07 2016

@author: earl
"""
import pandas as pd



#create header file
SAMPLE="headder"
outfile_one_by_background = open("OUT_DIR/recomb_background/%s.recombs.1.1.txt"  % SAMPLE, "w")
outline = ["Row","SAMPLE","Position","ERV-Mutations/ERV-Coverage","ERV Coverage","All Mutations/All Coverage","All Coverage"]
outline1 = '\t'.join([str(elem) for elem in outline]) + "\n"
outfile_one_by_background.write(outline1)

SAMPLES = ['22393',  '23039', '24050', '24072', '24088', '24443', 'BioClone']

#set threshholds frequencies for ERV variants to use for ERN read identification
ERV_SNP_thresh=0.05
SFFV_SNP_thresh=0.02
coverage = [0]*5000
ref=""
reference=open("REFS/SFFV_deNovo.fasta")
for line in reference:
    if line[0] == ">":
        None
    else:
        line= line.replace('\n','')
        ref+=line
        #print ("oranges")
print (ref[1977])

#dentify ERV variants and put them in a dict. Later, get intersection with SFFV variants while looping through samples
ERV_SNPS = {}
ERV_VCF=open("VCF_DIR/ERR216358_77.FVC.filter.vcf")
for line in ERV_VCF:
    if line[0] != "#":
        #print ("not header")
        l1 = line.replace('=','\t')
        l2 = l1.replace(';','\t')
        l3 = l2.split("\t")
        if float(l3[10]) >= ERV_SNP_thresh and l3[0] =="SFFV_deNovo":
            #print ("above threshold")
            ERV_SNPS[l3[1]] = l3[4]




for SAMPLE in SAMPLES:
    print (SAMPLE)
    
        
    #find intersection of ERV and SFFV snps
    ERV_SNPS_final = {}
    SFFV_VCF=open("VCF_DIR/%s.FVC.filter.vcf" % SAMPLE)
    for line in SFFV_VCF:
        if line[0] != "#":
            l1 = line.replace('=','\t')
            l2 = l1.replace(';','\t')
            l3 = l2.split("\t")
            if float(l3[10]) >= SFFV_SNP_thresh  and l3[0] =="SFFV_deNovo": 
                if l3[1] in ERV_SNPS:
                    if l3[4] == ERV_SNPS[l3[1]]:
                        ERV_SNPS_final[int(l3[1])-1] = ERV_SNPS[l3[1]] #fix an indexing error between samfile ingress and vcf. The vcf is correct, but it would take a while to fix and it doesnt change the analysis
                        
    
    #declare some variables to hold ERV/SFFV and novel variants
    is_erv = {} #ditionary with readname as key and [mut tpe, positon] as value
    novel_muts0 = []
    novel_muts1 = [[],[]] # keep 
    totals = {} #keep total coverages at reach copition


    
    
    #import mapped reads
    sam=open("ALIGNMENT_DIR/recomb_alignments/%s.realn.srt.rmd.sam" % SAMPLE)
    for line in sam:   
        l = line.replace("\n", "")
        l2 = l.split("\t")
        #only look at reads aligning to the concensus seq with a MQ > 100, and add aliases to make human readable
        if l2[2] == "SFFV_deNovo" and int(l2[4]) > 13:
            readname_bit=l2[0] + "_" + str(l2[1])
            seq=l2[9]
            pos=int(l2[3])-1
            CIGAR=str(l2[5])
            
            #correct reads to make them 1 to 1 alignment with ref
            if not "N" in CIGAR:
                CIG_str = ''
                holder = ''
                for l in range(0,len(CIGAR)):
                    if not CIGAR[l].isalpha():
                        holder += CIGAR[l]
                    else:
                        #if CIGAR[l] != "S":               
                        bases = int(holder)
                        CIG_str += bases*CIGAR[l]
                        holder = ''
                seq_correct = ''
                if True: #len(seq) != len(CIG_str):
                    j = 0 #keep track of number of Ds
                    for i in range(0,len(CIG_str)):
                        if CIG_str[i] == "M":
                            seq_correct += seq[i-j]
                        elif CIG_str[i] == "S":
                            None
                        elif CIG_str[i] == "I":
                            None
                        elif CIG_str[i] == "D":
                            seq_correct += "N"
                            j += 1
                          
                #create a list of [ref_pos,genome]'s for every concensus "N" covered by read.  genome = (0 if SFFV, and 1 if ERV) and 2 if novel (background)
                diff_pos=[]
                #keep track of above and novel mutants too.
                diff_pos_novel=[]
                
                "find if read has ERV_SNPs"
                i=1
                is_erv[readname_bit] = [[],[]]
                while i < len(seq_correct)-1:
                    if i+pos in ERV_SNPS_final:
                        if seq_correct[i] == "N":
                            None 
                            #print ("i+pos == 'N' ")
                            #print (readname_bit + "SFFV " + seq_correct[i-10:i+10] + " ref " + ref[i+pos-10:i+pos+10])
                        elif seq_correct[i] == ERV_SNPS_final[i+pos]:
                            
                            #if i+pos == 3209:
                                #print ("ERV  " + readname_bit + "SFFV " + seq_correct[i-5:i+5] + " ref " + ref[i+pos-10:i+pos+10])

                            is_erv[readname_bit][0].append(1)
                            is_erv[readname_bit][1].append(i+pos)
                        elif seq_correct[i] == ref[i+pos]:
                            #if i+pos == 3209:
                                #print ("SFFV   " + readname_bit + "SFFV " + seq_correct[i-5:i+5] + " ref " + ref[i+pos-10:i+pos+10])
                            is_erv[readname_bit][0].append(0)
                            is_erv[readname_bit][1].append(i+pos)
                        else:
                            #novel mutation, ignore for now
                            #print (readname_bit + " has novel at allele SNV site " + str(i+pos))
                            None
                    elif ref[i+pos] != seq_correct[i] and seq_correct[i] != "N": # test for novel mutant
                        is_erv[readname_bit][0].append(2)
                        is_erv[readname_bit][1].append(i+pos)
                    coverage[i+pos] +=1   
                    i += 1
            
    "quantify percent ERV_SNPs, then put in a category and count novel variants"
    for key in is_erv:
        zeros = is_erv[key][0].count(0)
        ones = is_erv[key][0].count(1)
        twos = is_erv[key][0].count(2)
        #if ones == (ones + zeros):
        
        j=1

        while j <= len(is_erv[key][0])-2:
            
            if is_erv[key][0][j] == 2:                
                if is_erv[key][0][j-1] == 1 and is_erv[key][0][j+1]:
                    novel_muts1.append([is_erv[key][0][j],is_erv[key][1][j]])
                if 0 not in is_erv[key][0][0:j] and 0 not in is_erv[key][0][j:len(is_erv[key][0])]:
                    if 1 in is_erv[key][0][0:j] and 1 in is_erv[key][0][j:len(is_erv[key][0])]:
                        novel_muts0.append([is_erv[key][0][j],is_erv[key][1][j]])
            
            j+=1
        
        #count totals at each position
        j=0
        while j <= len(is_erv[key][0])-1:
            if is_erv[key][1][j] not in totals:
                totals[is_erv[key][1][j]] = [0,0,0]
            
            if is_erv[key][0][j] == 2:
                totals[is_erv[key][1][j]][2] += 1
            elif is_erv[key][0][j] == 1:
                totals[is_erv[key][1][j]][1] += 1
            elif is_erv[key][0][j] == 0:
                totals[is_erv[key][1][j]][0] += 1
            j+=1

    "interpolate ERV frequency: replace zeros with linear interpolation of flanking values"            
    ERV_interpolate={}
    for I in range(0,4401):
        if I in totals and  totals[I][1] != 0:
            ERV_interpolate[I] = totals[I][1]
        else:
            ERV_interpolate[I] = None
    
    
                
    df = pd.DataFrame(list(ERV_interpolate.items()),columns = ['Position','value'])
    ERV_interpolate1 = df.interpolate(method ='linear', limit_direction ='forward') 
    df
    ERV_interpolate1
    ERV_int = ERV_interpolate1['value'].tolist()
    ERV_background = {}
    ERV_background1 = {}                
    for I in novel_muts0:
        if len(I) ==2: 
            if I[1] in ERV_background:
                ERV_background[I[1]] += 1
            else:
                ERV_background[I[1]] = 1
    
    
    for I in novel_muts1:
        if len(I) ==2: 
            if I[1] in ERV_background1:
                ERV_background1[I[1]] += 1
            else:
                ERV_background1[I[1]] = 1
    
    
    
    #should be in format:
    #["Row","SAMPLE","Position","ERV-Mutations/ERV-Coverage","ERV Coverage","All Mutations/All Coverage","All Coverage"]
    outfile_one_by_background = open("OUT_DIR/recomb_background/%s.recombs.1.1.txt"  % SAMPLE, "w")
    #outline = ["Position","ERV-Mutations/ERV-Coverage","ERV Coverage","All Mutations/All Coverage","All Coverage"]
    #outline1 = '\t'.join([str(elem) for elem in outline]) + "\n"
    #outfile_one_by_background.write(outline1)
    for I in range(1,4400):
        if I in ERV_background1:
            outline = [str(SAMPLE) + "_" +  str(I), SAMPLE, I,ERV_background1[I]/ERV_int[I],ERV_int[I],totals[I][2]/coverage[I],coverage[I]]
        elif I in totals:
            outline = [str(SAMPLE) + "_" +  str(I), SAMPLE, I,0,ERV_int[I],totals[I][2]/coverage[I],coverage[I]]
        else:
             outline = [str(SAMPLE) + "_" +  str(I), SAMPLE,I,0/ERV_int[I],ERV_int[I],0/coverage[I],0]
        outline1 = '\t'.join([str(elem) for elem in outline]) + "\n"
        outfile_one_by_background.write(outline1)
    outfile_one_by_background.close()
