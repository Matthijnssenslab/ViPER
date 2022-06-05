#!/usr/bin/env python3

import os, argparse, textwrap, shutil, sys
import pandas as pd
import checkv, pysam, pybedtools
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.MetavarTypeHelpFormatter,
    description=f"""This script takes the unclustered contigs from ViPER's triple assembly as input, 
    runs CheckV (https://bitbucket.org/berkeleylab/checkv) to split proviruses in 'host' and 'viral' contigs, and
    subsequently clusters the contigs. This script also ensures that contigs with contamination, duplications or contigs that
    are longer than expected, as assessed by CheckV, do not end up as (wrong) cluster representatives.
    """)
    parser.add_argument('-i', '--input', dest='fasta', required=True, type=str, metavar='PATH', help="Fasta file.")
    parser.add_argument('-d', '--database', dest='db', required=True, type=str, metavar='PATH', help="CheckV database.")
    parser.add_argument('-o', '--output', dest='output', required=True, type=str, metavar='STR', help="Output name for files.")
    parser.add_argument('-t', '--threads', dest='threads', type=int, metavar='INT', help="Number of threads to use.", default=1)
    parser.add_argument('-l', '--length', dest='length', type=int, metavar='INT', help="Minimum length of the host or viral region.", default=1000)
    parser.add_argument('--min-identity', dest='pid', type=int, metavar='INT', help="Minimum average nucleotide identity (ANI) for sequences to be clustered.", default=95)
    parser.add_argument('--min-coverage', dest='cov', type=int, metavar='INT', help="Minimum coverage %% of the shortest sequence that should be covered before clustering.", default=85)
    parser.add_argument('--keep-bed', dest='bed', action='store_true', help="Keep BED files with viral and host regions.", default=False)
    return vars(parser.parse_args())

def make_bed(contamination, output, minlength, keepBed=False):
    df=pd.read_csv(contamination, sep="\t", dtype=object)

    col_list=['contig_id', 'region_types', 'region_lengths', 'region_coords_bp']
    data=[]
    provirus_count=0
    count=0
    for index, row in df.iterrows():
        if row['provirus'] == "Yes":
            provirus_count+=1
            if row['region_types'] == "viral,host":
                data.append(row[col_list])
            elif row['region_types'] == "host,viral":
                data.append(row[col_list])
            elif row['region_types'] == "host,viral,host":
                data.append(row[col_list])
            else:
                count=+1
                if count == 1:
                    print(f"WARNING: Other types of contamination:")
                print(f"{row['contig_id']}")
    if provirus_count == 0:
        print("\n")
        print(f"No proviruses found in the data.")
        return None, None
    elif count > 0:
        print("\n")
        print(f"{count} contigs had another type of contamination.")
            
    df1=pd.DataFrame(data)
    
    df2=pd.concat([df1['contig_id'], 
                   df1[col_list[1]].str.split(',', expand=True), 
                   df1[col_list[2]].str.split(',', expand=True), 
                   df1[col_list[3]].str.split(',', expand=True)], 
                  axis=1).reset_index(drop=True)
    
    data=[]
    for i in range(3):
        try:
            df2[i]
        except:
            if i == 0:
                print(f"No contamination in contigs.")
        else:
            tmpdf=pd.concat([df2['contig_id'], df2[i]], axis=1)
            tmpdf.columns=col_list
            data.append(tmpdf)
    df3=pd.concat(data).sort_index().reset_index(drop=True)
    
    df4=pd.concat([df3.drop('region_coords_bp', axis=1), df3['region_coords_bp'].str.split('-', expand=True)], axis=1)
    df4.rename(columns={0:'start', 1:'end'}, inplace=True)
    df4[['region_lengths','start', 'end']] = df4[['region_lengths','start', 'end']].apply(pd.to_numeric)
    df4[['start', 'end']] -= 1
    
    vdata=[]
    hdata=[]
    for index, row in df4.iterrows():
        if row['region_types'] == "viral" and row['region_lengths'] >= minlength:
            row['bed_name']=row['contig_id'].replace('_length','v_length')
            row['bed_name']=row['bed_name'].replace('_cov',"v"+str(row['region_lengths'])+"_cov")
            vdata.append(row[['contig_id', 'start', 'end', 'bed_name']])
        elif row['region_types'] == "host" and row['region_lengths'] >= minlength:
            row['bed_name']=row['contig_id'].replace('_length','h_length')
            row['bed_name']=row['bed_name'].replace('_cov',"h"+str(row['region_lengths'])+"_cov")
            hdata.append(row[['contig_id', 'start', 'end', 'bed_name']])
    viral=pd.DataFrame(vdata)
    host=pd.DataFrame(hdata)
    
    if keepBed:
        viral.to_csv(output+"_viral.bed", sep="\t", header=False, index=False)
        host.to_csv(output+"_host.bed", sep="\t", header=False, index=False)
    
    return viral, host

def bedtools(bed, fasta):
    a=pybedtools.BedTool(bed, from_string=True)
    a = a.sequence(fi=fasta, name=True)
    return open(a.seqfn).read()

def bedtools_fasta_dict(fasta_text):
    fasta_dict={} #make dictionary
    for i in fasta_text.split("\n"):
        if i.startswith(">"):
            key=i.split(":")[0]
            fasta_dict[key]=""
        elif i == '':
            continue
        else:
            fasta_dict[key]=i
    return fasta_dict

def quality_summary_selection(checkv_summary):
    df=pd.read_csv(checkv_summary, sep="\t", dtype={'warnings':'str'})
    df['warnings']=df['warnings'].fillna("") #find faster way to get rid of Nan
    
    include=set()
    exclude=set()
    warnings=['high kmer_freq may indicate large duplication','contig >1.5x longer than expected genome length']
    for index, row in df.iterrows():
        if row['provirus']=="Yes" or any(x in row['warnings'] for x in warnings):#['kmer', 'longer']):
            exclude.add(row['contig_id'])
            
        if any(x in row['warnings'] for x in warnings): #['kmer', 'longer']):
            include.add(row['contig_id'])
        
    return include, exclude

def write_fasta(dictionary, name):
    wrapper = textwrap.TextWrapper(width=60)
    with open(name, "w") as fh:
        for k, v in  dictionary.items():
            if k.startswith(">"):
                fh.write(k+"\n")
                text_list=wrapper.wrap(text=v)
                for i in text_list:
                    fh.write(i+"\n")
            else:
                print("ERROR: dictionary key does not start with >, your fasta file might be misformatted.")
                sys.exit()

def biopython_fasta(dictionary):
    biopython_dict={}
    for k, v in dictionary.items():
        k=k.lstrip(">")
        biopython_dict[k]=SeqRecord(
            Seq(v),
            id=k,
            name=k,
            description=''
            )
    return biopython_dict

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
    os.mkdir("tmp_clustering/blastdb")
    makedb = NcbimakeblastdbCommandline(dbtype="nucl",
                                   input_file=fasta,
                                   out="tmp_clustering/blastdb/"+output+"_db")
    blastn = NcbiblastnCommandline(query=fasta,
                              db="tmp_clustering/blastdb/"+output+"_db",
                             outfmt='6 std qlen slen',
                             max_target_seqs=10000,
                             perc_identity=90,
                             num_threads=threads,
                             out="tmp_clustering/"+output+".out")
    makedb()
    blastn()
    
    anicalc_dict=anicalc("tmp_clustering/"+output+".out")

    keep=aniclust(fasta, anicalc_dict, "tmp_clustering/"+output, 
                  min_ani=pid, min_tcov=cov)
    
    with open(output+".fasta", "w") as f:
        for seq in SeqIO.parse(fasta, 'fasta'):
            if seq.id in keep:
                SeqIO.write(seq, f, 'fasta')

def main():
    args=parse_arguments()
    fasta=args['fasta']
    
    if not os.path.exists("tmp_clustering"):
        os.mkdir("tmp_clustering")
    
    checkv_arguments={"input": fasta,
         "db": args['db'],
         "output": "tmp_clustering/checkv",
         "threads": args['threads'],
         "restart": True,
         "quiet": False,
         "remove_tmp": True}
    
    print(f"Running CheckV...")
    
    checkv.end_to_end.main(checkv_arguments)
    
    contamination='tmp_clustering/checkv/contamination.tsv'
    qsummary='tmp_clustering/checkv/quality_summary.tsv'
    output=args['output']
    minlength=args['length']
    bed=args['bed']
    
    print(f"\nMaking BED files..")
    
    viralbed, hostbed = make_bed(contamination, output, minlength, keepBed=bed)
    
    if viralbed is None and hostbed is None:
        shutil.rmtree('tmp_clustering')
        sys.exit()
    
    print(f"Splitting host sequence from viral contigs...")
    pysam.faidx(fasta)
    host=bedtools(hostbed.to_csv(header=None, index=False, sep='\t'), fasta)
    viral=bedtools(viralbed.to_csv(header=None, index=False, sep='\t'), fasta)

    hdict=bedtools_fasta_dict(host)
    vdict=bedtools_fasta_dict(viral)
    
    print(f"Excluding contigs with contamination, longer than expected and duplication issues...")
    include, exclude = quality_summary_selection(qsummary)
    
    viral_exclude=set()
    for index, row in viralbed.iterrows():
        if row['contig_id'] in include:
            list_index=include.index(row['contig_id'])
            include[list_index]=row['bed_name']
            viral_exclude.add(row['bed_name'])
            
    clean_v_dict = {k: v for k, v in vdict.items() if k.lstrip(">") not in viral_exclude}
    
    fasta_seqs={}
    with open(fasta, "r") as fh:
        lines = fh.readlines()
        for line in lines:
            if line.startswith(">"):
                key=line.strip()
                fasta_seqs[key]=""
            else:
                fasta_seqs[key]+=line.strip()
    clean_fasta_dict = {k: v for k, v in fasta_seqs.items() if k.lstrip(">") not in exclude}

    
    cluster_dict={**clean_v_dict, **hdict, **clean_fasta_dict}
    write_fasta(cluster_dict, "tmp_clustering/"+output+"_cluster.fasta")
    
    print(f"Writing fasta file with contigs to re-include after cross-sample clustering...")
    inclv_dict={**vdict, **fasta_seqs}
    reinclude_dict = {k: v for k, v in inclv_dict.items() if k.lstrip(">") in include}
    write_fasta(reinclude_dict, output+"_re-include.fasta")  
    
    print(f"\nClustering contigs...")
    clustering("tmp_clustering/"+output+"_cluster.fasta", output, 
               args['threads'], args['pid'], args['cov'])
    
    shutil.rmtree('tmp_clustering')
    os.remove(fasta+".fai")

if __name__ == '__main__':
    main()