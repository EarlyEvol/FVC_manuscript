# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 14:50:07 2016

@author: earl
"""


from scipy.stats import binom_test, binom

SAMPLES = ['22393', '22578',  '23039', '24050', '24072', '24088', '24443', 'BioClone', 'ERR216358_77']



ref=""
reference=open("/home/earlyevol/FVC_manuscript/REFS/SFFV_chr12_conc.fasta")
for line in reference:
    if line[0] == ">":
        None
    else:
        ref+=line
        print ("oranges")
#print ref[3163:3178]

SFFV=""
SFFV_deNovo_aln=open("/home/earlyevol/FVC_manuscript/REFS/SFFV_deNovo_recomb.fasta")
for line in SFFV_deNovo_aln:
    if line[0] == ">":
        None
    else:
        SFFV+=line
        print ("apples")

chr12=""        
chr12_recomb_aln=open("/home/earlyevol/FVC_manuscript/REFS/chr12_recomb.fasta")
for line in chr12_recomb_aln:
    if line[0] == ">":
        None
    else:
        chr12+=line    

for SAMPLE in SAMPLES:
    print (SAMPLE)
    sam=open("/home/earlyevol/FVC_manuscript/recomb_alignments/%s.realn.srt.rmd.dnsmpl.sam" % SAMPLE)
    
    outfile_one_by_recomb = open("/home/earlyevol/FVC_manuscript/recomb_alignments/recombs/%s.assmbld.recombs_one_one.out.4.10.17"  % SAMPLE, "w")
    outfile_two_by_two_recomb = open("/home/earlyevol/FVC_manuscript/recomb_alignments/recombs/%s.assmbld.recombs_two_two.out.4.10.17"  % SAMPLE, "w")
    outfile_whole_read = open("/home/earlyevol/FVC_manuscript/recomb_alignments/recombs/%s.assmbld.recombs_whole_read.out.4.10.17"  % SAMPLE, "w")
    outfile_background = open("/home/earlyevol/FVC_manuscript/recomb_alignments/recombs/%s.assmbld.recombs_background.out.4.10.17"  % SAMPLE, "w")
    
    
    
    
    position_totals = {}
    for i in range(2300,3780):
        if ref[i] == "N":
            position_totals[i] = [[[0,0,0]],[0,0],"",""]
            
    recomb_pos_two_by_two = []
    background_pos = []
    bkgrnd_two_by_two = []
    recomb_whole_read = {}
    recomb_one_by_pos = []
    novel_background = []
    
    #import mapped reads
    for line in sam:   
        l = line.replace("\n", "")
        l2 = l.split("\t")
        #only look at reads aligning to the concensus seq with a MQ > 100, and add aliases to make human readable
        if l2[2] == "SFFV_chr12_Ns" and int(l2[4]) > 13:
            readname_bit=l2[0] + "_" + str(l2[1])
            seq=l2[9]
            pos=int(l2[3])-1
            CIGAR=str(l2[5])
            
            #!!!!!!!this is inconsistent.  The if D in CIGAR....then I have elif D in CIGAR.  Currently this will just ignore deletions and leave the missalignments downstream
            #skip reads with splits or deletions
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
                if len(seq) != len(CIG_str):
                    None
                else:
                    for i in range(0,len(seq)):
                        if CIG_str[i] == "M":
                            seq_correct += seq[i]
                        elif CIG_str[i] == "S":
                            None
                        elif CIG_str[i] == "I":
                            None
                        elif CIG_str[i] == "D":
                            seq_correct += "N"
                
                
                            
                #create a list of [ref_pos,genome]'s for every concensus "N" covered by read.  genome = (0 if SFFV, and 1 if ERV)
                diff_pos=[]
                #keep track of above and novel mutants too.
                diff_pos_novel=[]
                
                i=1
                while i < len(seq_correct)-1:
                    if ref[i+pos] == "N":
                        if seq_correct[i] == "N":
                            None        
                        elif seq_correct[i] == ref[i+pos]:
                            print ("EXPLODE, ATGC should not == 'N'")
                            None
                        elif seq_correct[i] == SFFV[i+pos]:
                            diff_pos_novel.append([i+pos,0])
                            diff_pos.append([i+pos,0])
                            position_totals[i+pos][0][0][0] += 1
                        elif seq_correct[i] == chr12[i+pos]:
                            diff_pos_novel.append([i+pos,1])
                            diff_pos.append([i+pos,1])
                            position_totals[i+pos][0][0][1] += 1
            
                        else:
                            #novel mutation, ignore for now
                            diff_pos.append([i+pos,"novel_mut"])
                            #print seq_correct
                            #print CIG_str
                            position_totals[i+pos][0][0][2] += 1
                    
                    elif ref[i+pos] != seq_correct[i]: # test for novel mutant
                        diff_pos_novel.append([i+pos,2])
                    i += 1
            
                    
                #find various types of break points
                for w in range(0,len(diff_pos)-1):
                    
                    #Whole read breakpoint, this needs some work.  I think it skips some reads for no reason
                    recomb_whole_read[readname_bit] = []
                    if len(diff_pos) >=4:
                        if diff_pos[w][1] != diff_pos[w + 1][1] and not "novel_mut" in diff_pos:
                            if diff_pos[w][1] != "novel_mut" and diff_pos[w + 1][1] !=  "novel_mut":      
                                recomb_whole_read[readname_bit].append(str(diff_pos[w][0]) + ":" + str(diff_pos[w + 1][0]) + "\t" + str(diff_pos[w][1]) + "\t" + str(diff_pos[w + 1][1]))
                    
                    #one by one recomb
                    if diff_pos[w][1] != diff_pos[w + 1][1]:
                        outline = str(diff_pos[w][0]) + ":" + str(diff_pos[w + 1][0]) + "\t" + str(diff_pos[w][1]) + "\t" + str(diff_pos[w + 1][1]) + "\t" + readname_bit + "\n"
                        recomb_one_by_pos.append(outline)
                    #tally up one_by background, outfile_background.write(outline)
                    outline = str(diff_pos[w][0]) + ":" + str(diff_pos[w + 1][0]) + "\t" + str(diff_pos[w][1]) + "\t" + str(diff_pos[w + 1][1]) + "\t" + readname_bit + "\n"
                    background_pos.append(outline)
            
                    #two by two recombs
                    if w < len(diff_pos) - 3:
                        if diff_pos[w][1] == diff_pos[w + 1][1] and diff_pos[w + 2][1] == diff_pos[w + 3][1] and diff_pos[w + 1][1] != diff_pos[w + 2][1]:
                            outline = str(diff_pos[w + 1][0]) + ":" + str(diff_pos[w + 2][0]) + "\t" + str(diff_pos[w + 1][1]) + "\t" + str(diff_pos[w + 2][1]) + "\t" + readname_bit + "\n"
                            recomb_pos_two_by_two.append(outline) 
                        outline = (str(diff_pos[w + 1][0]) + ":" + str(diff_pos[w + 2][0]) + "\t" + str(diff_pos[w + 1][1]) + "\t" + str(diff_pos[w + 2][1]) + "\t" + readname_bit + "\n")
                        bkgrnd_two_by_two.append(outline)
                        if diff_pos[w + 1][0] < 3170:            
                            None
                            #print str(diff_pos[w + 1][0]) + ":" + str(diff_pos[w + 2][0])   
                    
                    
                for w in range(0,len(diff_pos_novel)-1):
                    if diff_pos[w][1] != diff_pos[w + 1][1] and not "novel_mut" in diff_pos:
                            if diff_pos[w][1] != "novel_mut" and diff_pos[w + 1][1] !=  "novel_mut":      
                                recomb_whole_read[readname_bit].append(str(diff_pos[w][0]) + ":" + str(diff_pos[w + 1][0]) + "\t" + str(diff_pos[w][1]) + "\t" + str(diff_pos[w + 1][1]))
                    novel_background.append
    ###############################################################################  
              
    #tally up allele frequencies for variable sites                
    for key in position_totals:
        SFFV_nuc = position_totals[key][0][0][0]
        chr12_nuc = position_totals[key][0][0][1] 
        total = SFFV_nuc + chr12_nuc + position_totals[key][0][0][2]
        if total != 0:
            position_totals[key][1] = [SFFV_nuc/float(total) , chr12_nuc/float(total)]
    
    
    keys = []
    for key in sorted(position_totals):
        keys.append(key)
        
    ###############################################################################
    #calculate probablilities of one_by_one (01 and 10) and two_two recombinations (0011) and (1100)
    probs_one_by_dict = {}
    probs_two_by_dict = {}    
    for i in range(3,len(keys) - 3):
        print (i)
        probs_one_by_dict[keys[i]] = [[],[],"",""]
        probs_one_by_dict[keys[i]][2] = position_totals[keys[i]][1][0]* position_totals[keys[i + 1]][1][1]
        probs_one_by_dict[keys[i]][3] = position_totals[keys[i]][1][1]* position_totals[keys[i + 1]][1][0]
        
        probs_two_by_dict[keys[i]] = [[],[],"",""]
        probs_two_by_dict[keys[i]][2] = position_totals[keys[i - 1]][1][0]* position_totals[keys[i]][1][0]* position_totals[keys[i + 1]][1][1]* position_totals[keys[i + 2]][1][1]
        probs_two_by_dict[keys[i]][3] = position_totals[keys[i - 1]][1][1]* position_totals[keys[i]][1][1]* position_totals[keys[i + 1]][1][0]* position_totals[keys[i + 2]][1][0]
    
    #########################################################
    
    #put recomb_pos_two_by_two into a ditrionary with readname_bit as key
    #find 0011 or 1100 recombs
    recomb_quant_two_two = {}      
    for line in recomb_pos_two_by_two:
        if not "novel_mut" in line:
        
            line = line.split("\t")
            if line[0] in recomb_quant_two_two:
                recomb_quant_two_two[line[0]][0].append(int(line[1]))
                recomb_quant_two_two[line[0]][1].append(int(line[2]))
            else:
                recomb_quant_two_two[line[0]] = [[int(line[1])], [int(line[2])]]
    
    #find any breakpoint (01 or 10)
    #this need more work
    recomb_one_by = {}
    for line in recomb_one_by_pos:
        if not "novel_mut" in line:
            line = line.split("\t")
            if line[0] in recomb_one_by:
                recomb_one_by[line[0]][0].append(int(line[1]))
                recomb_one_by[line[0]][1].append(int(line[2]))
            else:
                recomb_one_by[line[0]] = [[int(line[1])], [int(line[2])]]
    
    
    #find reads with only one "recombination" and then make a dict with postion as key and populate it with [[5' genomes],[3' genomes]]
    whole_read_one_recomb = []
    whole_read_no_recomb = []
    whole_read_one_genome = []
    for key in recomb_whole_read:
        if len(recomb_whole_read[key]) == 1:
            whole_read_one_recomb.append(str(key + "\t" + "\t".join(recomb_whole_read[key])))
        if len(recomb_whole_read[key]) == 0:
            whole_read_no_recomb.append(str(key + "\t" + "\t".join(recomb_whole_read[key])))
        if len(recomb_whole_read[key]) < 1:
            None
            
        
        
    whole_read_one_recomb_dict = {}
    for line in whole_read_one_recomb:
        line = line.split("\t")
        if line[1] in whole_read_one_recomb_dict:
                whole_read_one_recomb_dict[line[1]][0].append(int(line[2]))
                whole_read_one_recomb_dict[line[1]][1].append(int(line[3]))
        else:
            whole_read_one_recomb_dict[line[1]] = [[int(line[2])], [int(line[3])]]
           
        
          
    
    #get some background measure of coverage for one-by-one and two-by-two
    recomb_quant_bkgnd_one_by_one = {}      
    for line in background_pos:
        if not "novel_mut" in line:
            line = line.split("\t")
            if line[0] in recomb_quant_bkgnd_one_by_one:
                recomb_quant_bkgnd_one_by_one[line[0]][0].append(int(line[1]))
                recomb_quant_bkgnd_one_by_one[line[0]][1].append(int(line[2]))
            else:
                recomb_quant_bkgnd_one_by_one[line[0]] = [[int(line[1])], [int(line[2])]]
    
                
    recomb_quant_bkgnd_two_by_two = {}      
    for line in bkgrnd_two_by_two:
        if not "novel_mut" in line:
            line = line.split("\t")
            if line[0] in recomb_quant_bkgnd_two_by_two:
                recomb_quant_bkgnd_two_by_two[line[0]][0].append(int(line[1]))
                recomb_quant_bkgnd_two_by_two[line[0]][1].append(int(line[2]))
            else:
                recomb_quant_bkgnd_two_by_two[line[0]] = [[int(line[1])], [int(line[2])]]
                
   
    
    
    #########################################################            
    
             
    #output two by two recombs stuff (this is changing a bunch) 
    for key in sorted(recomb_quant_two_two):
        pos = key.split(":")
        prob0011 = probs_two_by_dict[int(pos[0])][2]
        prob1100 = probs_two_by_dict[int(pos[0])][3]
        outline = key + "\t" +  str(float(sum(recomb_quant_two_two[key][0]))) + "\t" + str(float(sum(recomb_quant_two_two[key][1]))) + "\t" + str(prob0011) + "\t" + str(prob1100) + "\n"
        outfile_two_by_two_recomb.write(outline)       
        
    #output stats for all recombs (01 and 10)
    for key in sorted(recomb_one_by):
        pos = key.split(":")
        if len(recomb_quant_bkgnd_one_by_one[key][0]) >= 300:
            prob01 = probs_one_by_dict[int(pos[0])][2]
            prob10 = probs_one_by_dict[int(pos[0])][3]
            outline = key + "\t" +  str(float(sum(recomb_one_by[key][1]))) + "\t" + str(float(sum(recomb_one_by[key][0]))) + "\t" + str(prob01) + "\t" + str(prob10) + "\t" + str(len(recomb_quant_bkgnd_one_by_one[key][0])) + "\n"
            outfile_one_by_recomb.write(outline)
            print (key + "\t" + str(binom_test( ((sum(recomb_one_by[key][0]))), (len(recomb_quant_bkgnd_one_by_one[key][0])), (prob10))))
    
    
    #output absolute number of recombinations (background hard to calc)
    print ("whole read test")
    for key in sorted(whole_read_one_recomb_dict):
        pos = key.split(":")
        prob0011 = probs_two_by_dict[int(pos[0])][2]
        prob1100 = probs_two_by_dict[int(pos[0])][3]
        outline = key + "\t" + str(sum(whole_read_one_recomb_dict[key][0])) + "\t" + str(sum(whole_read_one_recomb_dict[key][1])) + "\t" + str(len(recomb_quant_bkgnd_one_by_one[key][0]))  + "\t" + str(prob0011) + "\t" + str(prob1100) + "\n"
        outfile_whole_read.write(outline)
        bintest = binom_test(((sum(whole_read_one_recomb_dict[key][0]))), (len(recomb_quant_bkgnd_one_by_one[key][0])), (prob1100))
        print ("\t".join([str(sum(whole_read_one_recomb_dict[key][0])), str(len(recomb_quant_bkgnd_one_by_one[key][0])),str(prob1100)]))
        print (key + "\t" + str(bintest))
    
print ("DONE")

"""
for key in sorted(recomb_one_by):      
    print key + "\t" + str(binom_test([(sum(recomb_one_by[key][0])), (len(recomb_quant_bkgnd_one_by_one[key][0]))], (prob01),'greater'))
    
"""