#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 14:33:35 2020

@author: helloalina
"""
import os
import argparse
import sys
import operator
import filecmp
from contextlib import redirect_stdout
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from shutil import copy

        
def run_blast(database,cmd,query_file,outdir,basename,maxseqs,seq_lengths,summaryfile,title_map,outformat=5,e=0.1**30): 
    matched_list = [1]
    alignment_info = []
    version = 0
    num_representations = {} #seq id: ----##----
    while matched_list != None:
        blastn_cline = NcbiblastnCommandline(cmd,query=query_file, db=database, outfmt=outformat, out=os.path.join(outdir,basename), evalue=e, max_target_seqs=maxseqs, max_hsps=1)
        stdout, stderr = blastn_cline()
        matched_list, num_representations = parse_blast(os.path.join(outdir,basename),alignment_info,version+1,outdir,seq_lengths, num_representations)
        if matched_list != None:
            filename_by_part = str(os.path.basename(query_file)).rsplit(".",1)
            copy(query_file, filename_by_part[0]+"Version"+str(version)+ "." + filename_by_part[-1]) #save a copy of the original query file
            for matched in matched_list:
                mask(query_file,matched)
            combine_annotations(num_representations,seq_lengths,version+1,summaryfile,outdir,title_map)
            version +=1
    copy(query_file, filename_by_part[0]+ "_finished." + filename_by_part[-1]) #save the remaining unmatched sequence to
    os.rename(filename_by_part[0]+"Version0."+filename_by_part[-1],query_file) #convert the query file back to the original (saved as  version 0)
    
    
    
def mask(query_file,matched): #mask a matched sequence
    f = open(query_file,"r")
    temp = ''
    seq = {}
    for line in f:
        if line.startswith(">"): 
            name = line.split()[0] #save the sequence name/number
            seq[name]= ''
        else:
            temp = line.replace('\n','') #remove whitespaces in the sequence to help search for the matched sequence
            loc = temp.upper().find(matched)
            if loc != -1:
                seq[name] += temp[:loc]+ "N"*len(matched)+temp[(loc+len(matched)):] #replace the matched sequence to N's
            else:
                seq[name] += temp
                
    f.close()
    f = open(query_file,"w")
    for k in seq.keys():
        f.write(k)
        f.write("\n")
        f.write(seq[k])
        f.write("\n")
    f.close()


def parse_blast(input_file, alignment_info,round_num, outdir,seq_lengths,num_representations):
    result_handle = open(input_file)
    blast_records = NCBIXML.parse(result_handle)
    aligned_seq = []
    unblasted_seq = []
    for k in seq_lengths.keys():
        unblasted_seq.append(k)
    for record in blast_records:
        q = str(record.query.split()[0]) #only get the seq id
        if record.alignments: 
            best_alignment =(record.alignments)[0] 
            hsp =  (best_alignment.hsps)[0]
            if hsp.identities/(hsp.query_end-hsp.query_start) > 0.95: 
                annotation_file_name = str(record.query.split()[0]) + "_annotation.log"
                with open(os.path.join(outdir,annotation_file_name),'a') as f: #make a separate file for annotations of different queries
                    with redirect_stdout(f):
                        unblasted_seq.remove(q)
                        if ((not num_representations) or (q not in num_representations.keys())):
                            num_representations[q] = "-"*100 
                        start = int(100*int(hsp.query_start)/seq_lengths[str(q)])
                        end = int(100*int(hsp.query_end)/seq_lengths[str(q)])
                        new_rep = (num_representations[q])[0:start] + str(round_num)*(end-start)+ (num_representations[q])[end:]
                        num_representations[q] = new_rep
                        print(round_num, hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,hsp.identities)
                        print(best_alignment.title)
                        aligned_seq.append(hsp.query)
                        alignment_info.append([record.query,hsp.query_start,hsp.query_end,round_num])
    if not aligned_seq:
        return(None, num_representations)
    else:
        for seq in unblasted_seq:
            annotation_file_name = seq + "_annotation.log"
            with open(os.path.join(outdir,annotation_file_name),'a') as f:
                f.write("{}: no blast hits\n".format(round_num))
        return(aligned_seq,num_representations)



def store_seq_len(queryfile,outdir): #store seqs sorted by their lens and make separate files for annotations
    seq_len_map = {}
    for record in SeqIO.parse(queryfile, "fasta"):
        seq_len_map[record.id] = len(record)
        annotation_file_name = str(record.id) + "_annotation.log"
        f = open(os.path.join(outdir,annotation_file_name),'w')
        f.close()
    sorted_seq_len_map = dict(sorted(seq_len_map.items(), key=operator.itemgetter(1),reverse=True))
    return sorted_seq_len_map


def combine_annotations(number_reps,seq_len_map,round_num,summaryfile,outdir,title_map):
    f = open(summaryfile,"a")
    f.write("***Round {} ***\n".format(round_num,))
    for seqid in seq_len_map.keys():   
        f.write("\n>{} ".format(seqid))
        for s in title_map[seqid]:
            f.write("{} ".format(str(s)))
        annotation_file_name = str(seqid) + "_annotation.log"
        f.write("\n"+number_reps[seqid]+"\n")
        with open(os.path.join(outdir,annotation_file_name),'r') as g:
             anno = g.readlines()
        f.writelines(anno)
    f.write("\n\n\n")
    f.close()
            
    
def info_beyond_id(txtfile):
    with open(txtfile, 'r') as file:
        extra_info = {}
        for line in file:
            line = line.split()
            if not line:  # empty line
                continue
            extra_info[line[0][1:]] = line[1:]
    return extra_info




def main():
    #parse the commands from the script file
    parser = argparse.ArgumentParser() 
    parser.add_argument("--queryfile",help="The name (or path) of the file to search for as query sequences")
    parser.add_argument("--outdir",help="Write the output files in this path as opposed to the default of standard output")
    parser.add_argument("--max_seqs", type=int, help="Only report HSPs for the best <integer> different subject sequences")
    parser.add_argument("--blastdb",help="The name of the database to search against")
    parser.add_argument("--blastcmd",help="The command used to launch the blastall executable.")
    parser.add_argument("--blast_outputfile_basename", help = "The base name of the blast result file")
    parser.add_argument("--stage", help = "The stage name")
    parser.add_argument("--sample", help = "The sample name")
    parser.add_argument("--titles", help = "The txt file containing the titles of the sequences")
    args = parser.parse_args()
    try:
        title_extra_info = info_beyond_id(args.titles)
        summary_file_name = args.sample + "_" + args.stage+ "_" + "RESULTS.txt"
        f = open(summary_file_name, 'w') #create an empty file for the result
        f.close()
        seq_len_map = store_seq_len(args.queryfile,args.outdir)
        src = args.queryfile
        dst = os.path.join(args.outdir,os.path.basename(args.queryfile))
        if not filecmp.cmp(src,dst):
            copy(src,dst)
        run_blast(database=args.blastdb,cmd=args.blastcmd,query_file=os.path.join(args.outdir,os.path.basename(args.queryfile)), maxseqs=args.max_seqs,outdir=args.outdir,basename = args.blast_outputfile_basename, seq_lengths = seq_len_map,summaryfile=summary_file_name, title_map = title_extra_info) #run blast
        print("Blast finished")
        f = open(os.path.join(args.outdir,"EXIT_SUCCESS.txt"),'w')
        f.close()
    except: # if encounter any error
        exception_type, exception_object, exception_traceback = sys.exc_info()
        line_number = exception_traceback.tb_lineno
        print(line_number)
        print(sys.exc_info()) #is this part needed?
        g = open(os.path.join(args.outdir,"EXIT_FAILURE.txt"),'w')
        g.close()
        sys.exit(1)





if __name__ == "__main__": main()