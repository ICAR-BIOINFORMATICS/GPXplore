import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox
from tkinter import font
import pandas as pd
from Bio import SeqIO

global g_id_list
global entry_count
global n_ids

class FileProcessorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GPXplore : Gene Promoter Extraction Tool")

        self.file_paths = [None, None, None]
        self.file_labels = ['Whole Genome (.fasta)', 'Genes Sequences (.fasta)', 'General Features (.gff)']
        self.output_folder = None

        title_font = font.Font(size=25, weight='bold',underline=False)
        title_font2 = font.Font(size=15, weight='bold',underline=True)
        title_font3 = font.Font(size=12, weight='bold',underline=True)

        self.output_folder_label = tk.Label(root, text="----|    .GPXplore.    |----", font=title_font)
        self.output_folder_label.grid(row=0, column=0,columnspan=3, padx=10, pady=10)
        
        self.output_folder_label = tk.Label(root, text="Gene Promoter Extraction Tool", font=title_font2)
        self.output_folder_label.grid(row=1, column=0,columnspan=3, padx=10, pady=5)

        self.output_folder_label = tk.Label(root, text="Select Input Files", font=title_font3)
        self.output_folder_label.grid(row=2, column=0,columnspan=2, padx=10, pady=5)
        
        self.file_buttons = []
        button = tk.Button(root, text=self.file_labels[0], command=lambda i=0: self.select_file(0))
        button.grid(row=3, column=0, padx=10, pady=5)
        self.file_buttons.append(button)

        button = tk.Button(root, text=self.file_labels[1], command=lambda i=1: self.select_file(1))
        button.grid(row=3, column=1, padx=10, pady=5)
        self.file_buttons.append(button)

        button = tk.Button(root, text=self.file_labels[2], command=lambda i=2: self.select_file(2))
        button.grid(row=4, column=0, padx=10, pady=5)
        self.file_buttons.append(button)
        
        self.output_folder_button = tk.Button(root, text="Select Output Folder", command=self.select_output_folder)
        self.output_folder_button.grid(row=4, column=1, padx=10, pady=5)

        tk.Label(root, text="Output File Name:").grid(row=5, column=0, padx=10, pady=5)
        self.output_file_entry = tk.Entry(root)
        self.output_file_entry.grid(row=5, column=1, padx=10, pady=5)
        self.output_file_entry.insert(tk.END, "gpx_output.csv")
      
        tk.Label(root, text="Upstream Length:").grid(row=6, column=0, padx=10, pady=5)
        self.up_len_entry = tk.Entry(root)
        self.up_len_entry.grid(row=6, column=1, padx=10, pady=5)
        self.up_len_entry.insert(tk.END, "1000")
      
        tk.Label(root, text="Downstream Length:").grid(row=7, column=0, padx=10, pady=5)
        self.dwn_len_entry = tk.Entry(root)
        self.dwn_len_entry.grid(row=7, column=1, padx=10, pady=5)
        self.dwn_len_entry.insert(tk.END, "1000")
      
        self.process_button = tk.Button(root, text="Extract Promoters", command=self.process_files)
        self.process_button.grid(row=8, column=0, columnspan=2, padx=10, pady=5)
        
        self.log_box = scrolledtext.ScrolledText(root, width=60, height=20)
        self.log_box.grid(row=9, column=0, columnspan=2, padx=10, pady=5)
        self.log_box.config(state='disabled')

        tk.Label(root, text="Development Team:\n1. Dr. Shbana Begam, ICAR-NIPB, New Delhi, India\n2. Dr. Samarth Godara, ICAR-IASRI, New Delhi, India", justify='left').grid(row=10, column=0,columnspan=3, padx=10, pady=5, sticky='w')


    def gpx(self, genome_path, gene_list_path, gff_path, out_path, upstrm_len, dwnstrm_len):

        self.log("\n\n Whole Genome Path:"+genome_path)
        self.log("\n\n Genes File Path:"+gene_list_path)
        self.log("\n\n GFF File Path:"+gff_path)
        self.log("\n\n Output file Path:"+out_path)
        self.log("\n\n Upstream Length:"+str(upstrm_len))
        self.log("\n\n Downstream Length:"+str(dwnstrm_len))

        def extract_gene_list(gene_list_path):

            gene_list=[]
            with open(gene_list_path, "r") as f:
                for line in f:
                    if line[0]=='>':
                        gene_list.append(line[1:-1].split(' ')[0])

            return gene_list

        self.log("Extracting gene list...")

        gene_list = extract_gene_list(gene_list_path)

        self.log("\n Number of gene IDs extracted:"+str( len(gene_list)))

        self.log("\n Extracting information from .gff file...")

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

        self.log("\n Number of entries in .gff file:"+str(gff_data.shape[0]))

        #clean the gene ids from the headers of the gene .fasta file

        g_list=[]
        for gene in gene_list:
            g_list.append(gene.replace(' ', '.').split('.')[0])

        g_list=gene_list

        #clean the gene ids from the .gff file for matching

        gff_data['id'] = gff_data['g_id']

        #select the row from the .gff file whose ids are matched from the gene .fasta file

        global g_id_list
        global entry_count
        g_id_list = []
        entry_count=0

        def select_gid(rec):
            global g_id_list
            global entry_count
            global n_ids
            entry_count+=1
            
            if entry_count%10000==0:
            	self.log(str(round((entry_count/n_ids)*100,2))+"%", same_line=True)
            
            rec_id = rec['id'].split(";")[0]
            for g in g_list:
                if g in rec_id:
                    g_id_list.append(g)
                    return True
            g_id_list.append("Unwanted")
            return False
            
        self.log("\n\n Filtering gene ids from .gff file... \n Please wait, it may take a few minutes...")
        
        global n_ids
        n_ids = gff_data.shape[0]
        self.log("\n\n Total IDs in .gff file:"+str(n_ids)+"\n\n IDs processed:\n\n")

        gff_data['select'] = gff_data.apply(select_gid, axis=1)
        gff_data['id'] = g_id_list

        #keep only the 'gene' and 'mRNA' type data from the .gff file

        gff_data2 = gff_data[(gff_data['select']) & ((gff_data['type']=='gene') | (gff_data['type']=='mRNA'))]

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

        self.log("\n Number of genes to be investigated:"+str(select_data.shape[0]))

        #extract the promotors from the genome file

        upstrm_list = ['-' for i in range(select_data.shape[0])]
        dwnstrm_list = ['-' for i in range(select_data.shape[0])]

        #process chromosome wise
        self.log("\n Extracting promoters' sequences...")

        def extract_promo(id_str, seq):
            #global global_seq_idx
            self.log("Processing id:"+id_str)
            self.log("Sequence length:"+str(len(seq)))
            
            s_data = select_data.values

            for i in range(len(upstrm_list)):
                #check if the chromosome is the correct one to pick the promotors from
                if s_data[i][0] in id_str:
                    
                    self.log("--Extrating Promotors for Gene id:"+str(s_data[i][4]))
                    
                    try:
                        #pick the promotor
                        up_start = int(s_data[i][2])-upstrm_len
                        up_end = int(s_data[i][2])
                        upstrm = seq[up_start:up_end]
                    except:
                        upstrm = '-'
                    upstrm_list[i] = upstrm

                    try:
                        dn_start = int(s_data[i][3])
                        dn_end = int(s_data[i][3])+dwnstrm_len
                        dwnstrm = seq[dn_start:dn_end]
                    except:
                        dwnstrm ='-'
                    dwnstrm_list[i] = dwnstrm

        total_chars = len(open(genome_path, 'r').read())
        read_chars = 0

        for record in SeqIO.parse(genome_path, "fasta"):
            self.log("\n\n Processing completed = "+str(round((read_chars/total_chars)*100, 2))+"%")
            self.log("\n\n Processing ID now:"+record.id)
            extract_promo(record.id, record.seq)
            read_chars += len(record.seq)

        #store the upstream and downstream sequences in the file

        select_data['upstream'] = upstrm_list
        select_data['dwnstream'] = dwnstrm_list

        #save the output file

        self.log("\n Saving the extracted information in file:"+out_path+"\n")

        select_data.to_csv(out_path, index=False)
        
        self.log("\n Output File Saved...\n\n Extraction Completed...")

    def select_file(self, file_index):
        file_path = filedialog.askopenfilename(title="Select File")
        if file_path:
            self.file_paths[file_index] = file_path
            self.log(f"{self.file_labels[file_index]}: {file_path}\n")
    
    def select_output_folder(self):
        self.output_folder = filedialog.askdirectory(title="Select Output Folder")
        if self.output_folder:
            #self.output_folder_label.config(text=self.output_folder)
            self.log(f"Selected output folder: {self.output_folder}\n")
    
    def process_files(self):

        up_len = int(self.up_len_entry.get())
        dwn_len = int(self.dwn_len_entry.get())
        output_file_name = self.output_file_entry.get()
        if not all(self.file_paths):
            messagebox.showerror("Error", "Please select all three files.")
            return
        if not self.output_folder:
            messagebox.showerror("Error", "Please select the output folder.")
            return
        if not output_file_name:
            messagebox.showerror("Error", "Please enter the output file name.")
            return

        self.process_button.config(state=tk.DISABLED)
        
        output_file_path = f"{self.output_folder}/{output_file_name}"
        
        self.gpx(self.file_paths[0], self.file_paths[1], self.file_paths[2], output_file_path, up_len, dwn_len)
        self.process_button.config(state=tk.NORMAL)
    
    def log(self, message, same_line=False):
        self.log_box.config(state='normal')
        if same_line:
            last_line_index = self.log_box.index('end-2c linestart')
            self.log_box.delete(last_line_index, 'end-1c')
            self.log_box.insert(tk.END, ' '+message)
        else:
            self.log_box.insert(tk.END, ' '+message + '\n')
        self.log_box.see(tk.END)
        self.log_box.config(state='disabled')
        self.log_box.update()

if __name__ == "__main__":
    root = tk.Tk()
    app = FileProcessorApp(root)
    root.mainloop()
