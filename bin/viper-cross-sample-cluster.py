#!/usr/bin/env python3

import os, argparse, textwrap, shutil, sys
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.MetavarTypeHelpFormatter,
    description=f"""
    """)
    parser.add_argument('-i', '--input', dest='fasta', required=True, type=str, metavar='PATH', help="Fasta file.")
    parser.add_argument('-o', '--output', dest='output', required=True, type=str, metavar='STR', help="Output name for files.")
    parser.add_argument('-t', '--threads', dest='threads', type=int, metavar='INT', help="Number of threads to use.", default=1)
    #parser.add_argument('-l', '--length', dest='length', type=int, metavar='INT', help="Minimum length of the host or viral region.", default=1000)
    parser.add_argument('--min-identity', dest='pid', type=int, metavar='INT', help="Minimum average nucleotide identity (ANI) for sequences to be clustered.", default=95)
    parser.add_argument('--min-coverage', dest='cov', type=int, metavar='INT', help="Minimum coverage %% of the shortest sequence that should be covered before clustering.", default=85)
    #parser.add_argument('--keep-bed', dest='bed', action='store_true', help="Keep BED files with viral and host regions.", default=False)
    return vars(parser.parse_args())

#CheckV Copyright (c) 2020, The Regents of the University of California,
#through Lawrence Berkeley National Laboratory (subject to receipt of 
#any required approvals from the U.S. Dept. of Energy). All rights reserved.
#===>
def parse_blast(handle):
    for line in handle:
        r = line.split()
        yield {
            'qname':r[0],
            'tname':r[1],
            'pid':float(r[2]),
            'len':float(r[3]),
            'qcoords':sorted([int(r[6]), int(r[7])]),
            'tcoords':sorted([int(r[8]), int(r[9])]),
            'qlen':float(r[-2]),
            'tlen':float(r[-1]),
            'evalue':float(r[-4])
            }

def yield_alignment_blocks(handle):
    # init block with 1st record
    key, alns = None, None
    for aln in parse_blast(handle):
        if aln['qname'] == aln['tname']:
            continue
        key = (aln['qname'], aln['tname'])
        alns = [aln]
        break
    # loop over remaining records
    for aln in parse_blast(handle):
        # skip self hits
        if aln['qname'] == aln['tname']:
            continue
        # extend block
        elif (aln['qname'], aln['tname']) == key:
            alns.append(aln)
        # yield block and start new one
        else:
            yield alns
            key = (aln['qname'], aln['tname'])
            alns = [aln]
    yield alns

def prune_alns(alns, min_len=0, min_evalue=1e-3, min_pid=90):
    # remove short aligns
    alns = [aln for aln in alns if aln['len'] >= min_len and aln['evalue'] <= min_evalue and aln['pid'] >= min_pid]
    return alns

def compute_ani(alns):
    return round(sum(a['len'] * a['pid'] for a in alns)/sum(a['len'] for a in alns),2)

def compute_cov(alns):

    # merge qcoords
    coords = sorted([a['qcoords'] for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        # overlapping, update start coord
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        # non-overlapping, append to list
        else:
            nr_coords.append([start, stop])
    # compute qry_cov
    alen = sum([stop - start + 1 for start, stop in nr_coords])
    qcov = round(100.0*alen/alns[0]['qlen'],2)

    # merge tcoords
    coords = sorted([a['tcoords'] for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        # overlapping, update start coord
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        # non-overlapping, append to list
        else:
            nr_coords.append([start, stop])
    # compute qry_cov
    alen = sum([stop - start + 1 for start, stop in nr_coords])
    tcov = round(100.0*alen/ alns[0]['tlen'],2)

    return qcov, tcov
#<===

def anicalc(blast_input, pid=90, length=0):
    with open(blast_input, 'r') as input:
        fields = ['qname', 'tname', 'num_alns', 'pid', 'qcov', 'tcov']
        anicalc_dict={}
        index=0
        for alns in yield_alignment_blocks(input):
            alns = prune_alns(alns, min_pid=pid, min_len=length)
            if len(alns) == 0: 
                continue
            index+=1
            qname, tname = alns[0]['qname'], alns[0]['tname']
            ani = compute_ani(alns)
            qcov, tcov = compute_cov(alns)
            anicalc_dict[index] = {'qname': qname, 
                           'tname': tname, 
                           'num_alns': len(alns), 
                           'pid': ani, 
                           'qcov': qcov, 
                           'tcov': tcov}
        return anicalc_dict

def parse_seqs(path):
    handle = gzip.open(path) if path.split('.')[-1] == 'gz' else open(path)
    id = next(handle).split()[0][1:]
    seq = ''
    for line in handle:
        if line[0] == '>':
            yield id, seq
            id = line.split()[0][1:]
            seq = ''
        else:
            seq += line.rstrip()
    yield id, seq
    handle.close()

def aniclust(fasta, anicalc_dict, #outname, 
             minlength=0, min_qcov=0, 
             min_tcov=85, min_ani=95):
    # list seqs, sorted by length
    print("\nreading sequences...")
    seqs = {}
    
    for index, r in enumerate(parse_seqs(fasta)):
        id, seq  = r
        if len(seq) < minlength:
            continue
        else:
            seqs[id] = len(seq)
            
    seqs = [x[0] for x in sorted(seqs.items(), key=lambda x: x[1], reverse=True)]
    print("%s sequences retained from fna" % len(seqs))
    
    # store edges
    print("\nstoring edges...")
    num_edges = 0
    edges = dict([(x,[]) for x in seqs])
    
    for i in anicalc_dict:
        if anicalc_dict[i]['qname'] == anicalc_dict[i]['tname']:
            continue
        elif anicalc_dict[i]['qname'] not in edges or anicalc_dict[i]['tname'] not in edges:
            continue
        elif float(anicalc_dict[i]['qcov']) < min_qcov or float(anicalc_dict[i]['tcov']) < min_tcov or float(anicalc_dict[i]['pid']) < min_ani:
            continue
        edges[anicalc_dict[i]['qname']].append(anicalc_dict[i]['tname'])
        num_edges+=1
    

    print("%s edges retained from blastani" % num_edges)
    print("%s edges currently stored" % sum([len(_) for _ in edges.values()]))
    
    # cluster
    print("\nclustering...")
   
    clust_to_seqs = {}
    seq_to_clust = {}
    # loop over seqs in sorted order
    for seq_id in seqs:
    # seq already been assigned; cant be centroid
        if seq_id in seq_to_clust:
            continue
    # seq is centroid for new cluster
        else:
        # add self to cluster
            clust_to_seqs[seq_id] = [seq_id]
            seq_to_clust[seq_id] = seq_id
        # update with cluster members
            for mem_id in edges[seq_id]:
                if mem_id not in seq_to_clust:
                    clust_to_seqs[seq_id].append(mem_id)
                    seq_to_clust[mem_id] = seq_id
                    
    print("%s total clusters" % len(clust_to_seqs))
    
    # write
    print("\nwriting clusters...")
    keep=set()
    for seq_id, mem_ids in clust_to_seqs.items():
        keep.add(seq_id)
    #with open(outname, 'w') as out:
    #    for seq_id, mem_ids in clust_to_seqs.items():
    #        keep.add(seq_id)
    #        out.write(seq_id + '\t' + ','.join(mem_ids)+'\n')
    
    return keep

def clustering(fasta, output, threads, pid, cov):
    os.mkdir("blastdb")
    makedb = NcbimakeblastdbCommandline(dbtype="nucl",
                                   input_file=fasta,
                                   out="blastdb/"+output+"_db")
    blastn = NcbiblastnCommandline(query=fasta,
                              db="blastdb/"+output+"_db",
                             outfmt='6 std qlen slen',
                             max_target_seqs=10000,
                             perc_identity=90,
                             num_threads=threads,
                             out=output+".out")
    makedb()
    blastn()
    
    anicalc_dict=anicalc(output+".out")

    keep=aniclust(fasta, anicalc_dict, output, 
                  min_ani=pid, min_tcov=cov)
    
    with open(output+".fasta", "w") as f:
        for seq in SeqIO.parse(fasta, 'fasta'):
            if seq.id in keep:
                SeqIO.write(seq, f, 'fasta')

def main():
    args=parse_arguments()
    cluster_fasta=args['cluster_fasta']
    reinclude_fasta=args['reinclude_fasta']
    output=args['output']
    threads=args['threads']


    clustering(cluster_fasta, output, threads)

    NcbimakeblastdbCommandline(dbtype="nucl",
                                   input_file=reinclude_fasta,
                                   out="blastdb/"+output+"_reinclude_db")

    NcbiblastnCommandline(query=reinclude_fasta,
                              db="blastdb/"+output+"_reinclude_db",
                             outfmt='6 std qlen slen',
                             max_target_seqs=10000,
                             perc_identity=90,
                             num_threads=threads,
                             out=output+"_reinclude.out")
    
    
    anicalc_dict=anicalc(output+"_reinclude.out")

    qcov85=set()
    tcov85=set()
    singletons=set()
    
    for i in anicalc_dict:
        if anicalc_dict[i]['ani'] >= 95 and anicalc_dict[i]['qcov'] >= 85:
            if anicalc_dict[i]['qname'] not in qcov85:
                qcov85.add(anicalc_dict[i]['qname'])        
        elif anicalc_dict[i]['ani'] >= 95 and anicalc_dict[i]['tcov'] >= 85:
            if anicalc_dict[i]['qname'] not in qcov85.union(tcov85):
                tcov85.add(anicalc_dict[i]['qname'])
        else:
            if anicalc_dict[i]['qname'] not in qcov85.union(tcov85, singletons):
                singletons.add(anicalc_dict[i]['qname'])
    


