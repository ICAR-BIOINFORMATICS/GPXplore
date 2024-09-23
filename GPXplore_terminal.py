#python GPX_2.py --genome_path "/home/samarth/Downloads/GPX_2/Exp_Cucurbit/Cmoschata_genome.fasta" --gene_path "/home/samarth/Downloads/GPX_2/Exp_Cucurbit/gene.fasta" --gff_path "/home/samarth/Downloads/GPX_2/Exp_Cucurbit/Cmoschata.gff3" --output_path "/home/samarth/Downloads/GPX_2/Exp_Cucurbit/output.csv" --upstream_len 50 --downstream_len 50

import argparse
import pandas as pd
from Bio import SeqIO

parser=argparse.ArgumentParser(description="Add path")
parser.add_argument("--genome_path")
parser.add_argument("--gene_path")
parser.add_argument("--gff_path")
parser.add_argument("--output_path")
parser.add_argument("--upstream_len", type=int)
parser.add_argument("--downstream_len", type=int)
args=parser.parse_args()

upstrm_len = 2000
dwnstrm_len = 2000
gene_list_path = args.gene_path
gff_path=args.gff_path
genome_path=args.genome_path

out_path = args.output_path
if args.upstream_len:
    upstrm_len = args.upstream_len
if args.downstream_len:
    dwnstrm_len = args.upstream_len

#obtain the list of genes from the gene fasta file

def extract_gene_list(gene_list_path):

    gene_list=[]
    with open(gene_list_path, "r") as f:
        for line in f:
            if line[0]=='>':
                gene_list.append(line[1:-1].split(' ')[0])

    return gene_list

print("\n\n*************************************************************************\n\n\
 	\t\t    Welcome to \n\n\
   \t\tGPXplore: Gene Promoter Extraction Tool\n\n\
*************************************************************************\n\n\
    Documentation Link: https://github.com/ICAR-BIOINFORMATICS/GPXplorer \n\n\
*************************************************************************\n\n\
    Institute: ICAR-NIPB & ICAR-IASRI, New Delhi, India\n\n\
    Developed by: \n\
    \tDr. Shbana Begam, ICAR-NIPB, New Delhi, India \n\
    \tDr. Samarth Godara, ICAR-IASRI, New Delhi, India\n\n\
*************************************************************************\n\n")

print("Extracting gene list...")

gene_list = extract_gene_list(gene_list_path)

print("\nNumber of gene IDs extracted:", len(gene_list))

print("\nExtracting information from .gff file...")

#obtain the gene information from the .gff file

rec_count=0
recs = []
with open(gff_path, "r") as f:
    for line in f:
        if line[0]!='#':
            att = line.split('\t')
            recs.append(att)

#create dataframe from the .gff file

chr_list =[]
type_list =[]
start_list =[]
end_list =[]
id_list =[]

for rec in recs:
    chr_list.append(rec[0])
    type_list.append(rec[2])
    start_list.append(rec[3])
    end_list.append(rec[4])
    id_list.append(rec[-1])

gff_data = pd.DataFrame()
gff_data['chr']=chr_list
gff_data['type']=type_list
gff_data['start']=start_list
gff_data['end']=end_list
gff_data['g_id']=id_list

print("\nNumber of entries in .gff file:", gff_data.shape[0])

print("\n\nWith raw gene id tags (samples): \n", gff_data[:5])

#clean the gene ids from the headers of the gene .fasta file

g_list=[]
for gene in gene_list:
    g_list.append(gene.replace(' ', '.').split('.')[0])

g_list=gene_list

print("Extracted gene list from (samples):\n", gene_list[:5])

print("\n\nExtracted gene list samples:\n", g_list[:5])

#clean the gene ids from the .gff file for matching

gff_data['id'] = gff_data['g_id']

#select the row from the .gff file whose ids are matched from the gene .fasta file

g_id_list = []
entry_count=0

def select_gid(rec):
    global g_id_list
    global entry_count
    global n_ids
    entry_count+=1
    
    if entry_count%10000==0:
    	print(round((entry_count/n_ids)*100,2),"%")
    
    rec_id = rec['id'].split(";")[0]
    for g in g_list:
        if g in rec_id:
            g_id_list.append(g)
            return True
    g_id_list.append("Unwanted")
    return False
    
print("\n\nFiltering gene ids from .gff file... Please wait, it may take a few minutes...")

n_ids = gff_data.shape[0]
print("\n\nTotal IDs in .gff file:",n_ids,"\n\nIDs processed:")

gff_data['select'] = gff_data.apply(select_gid, axis=1)
gff_data['id'] = g_id_list

print("\n\ngff data with selection mark:\n", gff_data.head())

print("\n\ngff selection mark distribution:\n", gff_data.select.value_counts())

#keep only the 'gene' and 'mRNA' type data from the .gff file

gff_data2 = gff_data[(gff_data['select']) & ((gff_data['type']=='gene') | (gff_data['type']=='mRNA'))]

print("\n\ngff data of the selected rows, i.e. gene and mRNA samples:\n", gff_data2.head())

#from the 'gene' and 'mRNA' record, select the 'gene' type ,if not present then select the 'mRNA' type.

id_grps = gff_data2.groupby(['id']) 

selected_list = []
for gid, grp in id_grps:
    selected_list.append(grp.values[0])

#create a new dataframe for the selected gene records from the .gff file    

chr_list =[]
type_list =[]
start_list =[]
end_list =[]
id_list =[]    

for item in selected_list:
    chr_list.append(item[0])
    type_list.append(item[1])
    start_list.append(item[2])
    end_list.append(item[3])
    id_list.append(item[5])

select_data = pd.DataFrame()
select_data['chr']=chr_list
select_data['type']=type_list
select_data['start']=start_list
select_data['end']=end_list
select_data['g_id']=id_list

print("\n\nSelected gene records from .gff file:\n", select_data.head())

print("\nNumber of genes to be investigated:",select_data.shape[0])

#extract the promotors from the genome file

upstrm_list = ['-' for i in range(select_data.shape[0])]
dwnstrm_list = ['-' for i in range(select_data.shape[0])]
#process chromosome wise
print("\nExtracting promoters' sequences...")

def extract_promo(id_str, seq):
    #global global_seq_idx
    print("Processing id:", id_str)
    print("Sequence length:", len(seq))
    
    s_data = select_data.values

    for i in range(len(upstrm_list)):
        #check if the chromosome is the correct one to pick the promotors from
        if s_data[i][0] in id_str:
            
            print("--Extrating Promotors for Gene id:", s_data[i][4])
            
            try:
                #pick the promotor
                up_start = int(s_data[i][2])-upstrm_len
                up_end = int(s_data[i][2])
                upstrm = seq[up_start:up_end]
                print("Starting Index:", up_start)
                print("Ending Index:", up_end)
                print("Upstream:",upstrm)
            except:
                upstrm = '-'
            upstrm_list[i] = upstrm

            try:
                dn_start = int(s_data[i][3])
                dn_end = int(s_data[i][3])+dwnstrm_len
                dwnstrm = seq[dn_start:dn_end]
                print("Starting Index:", dn_start)
                print("Ending Index:", dn_end)
                print("Downstream:",dwnstrm)
            except:
                dwnstrm ='-'
            dwnstrm_list[i] = dwnstrm
        
#genome file path
total_chars = len(open(genome_path, 'r').read())
read_chars = 0
seq= ''
id=''

print("\n\nIDs in the genome fasta file:\n")

for record in SeqIO.parse(genome_path, "fasta"):
    print("\n\nProcessing completed = ",round((read_chars/total_chars)*100, 2),"%")
    extract_promo(record.id, record.seq)
    read_chars += len(record.seq)

#store the upstream and downstream sequences in the file

select_data['upstream'] = upstrm_list
select_data['dwnstream'] = dwnstrm_list

#save the output file

print("\nSaving the extracted infromation in file:", out_path,"\n")

select_data.to_csv(out_path, index=False)
