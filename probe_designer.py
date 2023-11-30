import os
from pathlib import Path
import subprocess
import re
import csv
import ast
import pandas as pd


### This is a demonstration

def design_probes(genomes_dir,accession_list,output_dir,run_name,probe_props=None,max_nonspecific_match=18,selection_method='dp',is_allow_duplicates=True,min_pid=95, min_len_similarity_percent=90, min_percent_cov=90):
    """
    genomes_dir: where to find the sequences
    accession_list: what genomes to use (if more then one also ok)
    output_dir: where to place the output files
    run_name: the output dir name
    """
    run_dir = rf'{output_dir}/{run_name}'

    print('Collecting sequences and annotations')
    fna_db = collect_genomes_for_probe_design(genomes_dir,accession_list,output_dir,run_name)

    print('Run makeblastdb')
    run_makeblastdb(fna_db, 'nucl', 'merged_genes')

    print('Generate naive probes with sliding window')
    generate_naive_probes(run_dir,probe_props)
    
    print('blasting naive probes')
    blast_naive_probes_blast(run_dir,word_size=max_nonspecific_match)

    print('identify putative gene duplications')
    duplicated_genes_db = identify_duplicate_genes(run_dir,min_pid=95, min_len_similarity_percent=90, min_percent_cov=90)

    print('Determine probe specificity')
    identify_specific_probes(run_dir,duplicated_genes_db,max_nonspecific_match = max_nonspecific_match,is_allow_duplicates=True)

    return



def calculate_probes_per_gene(run_dir,min_distance=0, overlap=0, method='dp'):
    """" 
    method = dp (dynamic programming based), or heuristic (sort and select)
    """

    if not method == 'dp' and not method == 'heuristic':
         print('ERROR: method = dp (dynamic programming based), or heuristic (sort and select)')
         return

    selected_probes_dp = load_and_process(rf'{run_dir}/specific_probes_table.txt', min_distance=min_distance, overlap=overlap, method=method)
    save_selected_probes(rf'{run_dir}/specific_probes_{method}_min_{min_distance}_overlap_{overlap}.txt', selected_probes_dp)
    save_locus_probes_count(rf'{run_dir}/probes_per_gene_summary_{method}_min_{min_distance}_overlap_{overlap}.txt', rf'{run_dir}/merged_gene_info.txt', selected_probes_dp,run_dir)

    return

def run_makeblastdb(input_file, db_type, output_file):
    # Save the current working directory
    current_dir = os.getcwd()

    # Change the working directory to the location of the input file
    os.chdir(Path(input_file).parent)

    cmd = f'makeblastdb -in "{Path(input_file).name}" -dbtype {db_type} -out "{output_file}"'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    os.chdir(current_dir)
    return

def collect_genomes_for_probe_design(genomes_dir,assembly_list,output_dir,run_name):
    """
    genomes_dir: where to find the downloaded genomic data 
    assembly_list: list of assembly identifiers 
    """
    output_path = rf'{output_dir}/{run_name}'
    os.makedirs(rf'{output_dir}/{run_name}', exist_ok=True)
    downloaded_genomes_list = os.listdir(genomes_dir)
    output_fna_path = rf'{output_path}/merged_genes.fna'

    with open(output_fna_path,'w') as o:
        with open(rf'{output_path}/merged_gene_info.txt','w') as genes_o:
            for genome in downloaded_genomes_list:
                for assembly in assembly_list:
                    if genome.endswith(assembly):
                        with open(rf'{genomes_dir}/{genome}/merged_genes.fasta','r') as f:
                            for line in f:
                                o.write(line)
                        with open(rf'{genomes_dir}/{genome}/{assembly}.genes_info.txt','r') as f:
                            for line in f:
                                genes_o.write(line)
        
    return output_fna_path

def read_fasta_file(fasta_file):
	db = {}
	current_header = ''
	with open(fasta_file,'r') as fa:
		for line in fa:
			if re.search('^>',line):
				current_header = line.rstrip()
				db[current_header] = ''
			else:
				db[current_header] += line.rstrip()
	return db	

def calculate_probe_GC_content(probe_seq):
    AT, GC, Ns = 0, 0, 0
    for i in probe_seq:
        if i in 'AT': AT += 1
        elif i in 'GC': GC += 1
        elif i == 'N': Ns += 1 
    probe_GC_content = round(GC / (AT + GC) * 100, 1)  # Exclude 'N's from the calculation
    return probe_GC_content


def generate_naive_probes(run_dir,probe_props=None): #produce all possible probes with sliding window approach
    
    if probe_props is None:
        probe_props = {
        'probe_len' : 30,
        'gc_min' : 40,
        'gc_max' : 65,
        'max_base_rep' : 4
        }
        print('Using standard probe design paramters...')
        print(probe_props)

    fna_db = read_fasta_file(rf'{run_dir}/merged_genes.fna')
    naive_probes_path = rf'{run_dir}/naive_probes.fna'

    probe_len,gc_min,gc_max,max_base_rep = int(probe_props['probe_len']),int(probe_props['gc_min']),int(probe_props['gc_max']),int(probe_props['max_base_rep'])

    with open(naive_probes_path,'w') as probes_out_fh:
        for header,gene_seq in fna_db.items():
            header = re.sub('>','',header)
            locus_tag,contig,assembly = header.split(';')
            
            for i in range(0,len(gene_seq)-probe_len,1): #run over gene sequence, base by base and test each potential probe
                sub_seq = str(gene_seq[i:i+probe_len])
                sub_seq_gc = calculate_probe_GC_content(sub_seq)
                if sub_seq_gc >= gc_min and sub_seq_gc <= gc_max:
                    if not re.search("A{%s}|G{%s}|C{%s}|T{%s}"%(4*(max_base_rep,)),sub_seq):
                        probe_fr,probe_to = str(i+1),str(i+probe_len)
                        probe_header = f'{locus_tag};{contig};{assembly};{sub_seq_gc};{probe_fr}-{probe_to}'
                        probes_out_fh.write(f'>{probe_header}\n{sub_seq}\n')
    return 

def blast_naive_probes_blast(run_dir,word_size=18,evalue=100,n_threads=8,output_format=6):
    current_dir = os.getcwd()
    os.chdir(run_dir)

    query = './naive_probes.fna'
    db = './merged_genes'
    output_path = './naive_probes.blast_results'

    cmd = f'blastn -query {query} -db {db} -word_size {word_size} -evalue {evalue} -outfmt {output_format} -out {output_path}'

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    os.chdir(current_dir)
    return

def identify_duplicate_genes(run_dir,min_pid=95,min_len_similarity_percent=90,min_percent_cov=90): # all probes hitting these will be otherwise considered non-specif

	current_dir = os.getcwd()
	os.chdir(run_dir)
	fasta_db = read_fasta_file(rf'{run_dir}/merged_genes.fna')
	self_blast_path = r'gene_duplicate_identification.blast_res'
	query = r'merged_genes.fna'
	db = r'merged_genes'
	cmd = f'blastn -query {query} -db {db} -word_size 28 -evalue 100 -outfmt 6 -out {self_blast_path}'
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	
	duplicated_genes_db = {}
	with open(self_blast_path,'r') as blast_res_fh:
		for line in blast_res_fh:
			query_gene_header,hit_gene_header,pid,match_len,mismatches,gap_len,q_start,q_end,h_start,h_end,evalue,bitscore = line.rstrip().split('\t')
		
			if query_gene_header == hit_gene_header: continue # same gene hit

			q_locus,q_contig,q_assembly = query_gene_header.split(';') #org_id,locus_tag,gene,desc
			h_locus,h_contig,h_assembly = hit_gene_header.split(';')

			#don't allow cross species
			if not q_assembly == h_assembly:
				continue 
			if q_locus == h_locus:
				continue
			
			# Is the percent identity significant + is the gene size similar enough?
			is_duplicate = False
			if float(pid) >= min_pid:
				q_gene_len = len(fasta_db['>%s'%query_gene_header])
				h_gene_len = len(fasta_db['>%s'%hit_gene_header])
				gene_len_ratio = min(q_gene_len,h_gene_len) / max(q_gene_len,h_gene_len) * 100
				match_percent_cov = float(match_len) / min(q_gene_len,h_gene_len)  * 100 # percent cov of smallest gene	
				#genes of similar length?
				if gene_len_ratio >= min_len_similarity_percent and match_percent_cov >= min_percent_cov:					
					is_duplicate = True
			
			if is_duplicate:
				if not query_gene_header in duplicated_genes_db:
					duplicated_genes_db[query_gene_header] = [hit_gene_header]
				else:
					duplicated_genes_db[query_gene_header].append(hit_gene_header)
				
				with open('./duplicated_genes.txt','w') as f:
					for k,v in duplicated_genes_db.items():
						f.write(f'{k}\t{v}\n')
	os.chdir(current_dir)
	return	duplicated_genes_db
					


def identify_specific_probes(run_dir,duplicated_genes_db,max_nonspecific_match = 18,is_allow_duplicates=True):
    
    naive_probe_fasta_db = read_fasta_file(rf'{run_dir}/naive_probes.fna')
    probe_result_db = {}
    with open(rf'{run_dir}/naive_probes.blast_results','r') as f:
        for line in f: #BSU_00010;NC_000964.3;ASM904v1;40.0;6-35
            probe_header,hit_header,pid,match_len,mismatches,gap_len,q_start,q_end,h_start,h_end,evalue,bitscore = line.rstrip().split('\t')
            probe_locus,probe_contig,probe_assembly,probe_gc,probe_fr_to = probe_header.split(';')
            
            probe_fr,probe_to = probe_fr_to.split('-')

            hit_locus,hit_contig,hit_assembly = hit_header.split(';')
            is_duplicate_hit = False

            probe_header_split = probe_header.split(';') #'b1493;NC_000913.3;ASM584v2;43.3;22-51'.split(';') 
            probe_header_trim = ';'.join(probe_header_split[0:3]) 

            if hit_header in duplicated_genes_db and probe_header_trim in duplicated_genes_db[hit_header]:
                 is_duplicate_hit = True

            #override in case this option is set to False
            if is_allow_duplicates == False:
                 is_duplicate_hit = False

            # register the probe if not in the db
            if not probe_header in probe_result_db:
                probe_result_db[probe_header] = {
                    'probe_seq' : naive_probe_fasta_db[f'>{probe_header}'],
                    'number_of_non_specific_hits' : 0,
                    'all_non_specific_hits' : [],
                    'probe_assembly' : probe_assembly,
                    'probe_locus' : probe_locus,
                    'probe_gc' : probe_gc,
                    'probe_fr' : probe_fr,
                    'probe_to' : probe_to,
                    'probe_desc' : 'NA'
                }
            #is specific hit?
            if int(match_len) >= max_nonspecific_match:
                if not probe_locus == hit_locus and not is_duplicate_hit: #hit on another gene ok if its your double
                    probe_result_db[probe_header]['number_of_non_specific_hits'] += 1
                    probe_result_db[probe_header]['all_non_specific_hits'].append(hit_header)
            
    
    #save the specific probes 
    with open(fr'{run_dir}/specific_probes_table.txt','w') as o:
        for probe_header in probe_result_db:
            number_of_non_specific_hits = probe_result_db[probe_header]['number_of_non_specific_hits']
            # if probe is completely specific
            if number_of_non_specific_hits == 0:
                probe_assembly = probe_result_db[probe_header]['probe_assembly']
                probe_locus = probe_result_db[probe_header]['probe_locus']
                probe_fr = probe_result_db[probe_header]['probe_fr']
                probe_to = probe_result_db[probe_header]['probe_to']
                probe_gc = probe_result_db[probe_header]['probe_gc']
                probe_seq = probe_result_db[probe_header]['probe_seq']
                line = '\t'.join([probe_assembly,probe_locus,probe_fr,probe_to,probe_gc,probe_seq])
                o.write(line + '\n')
    return 


### select the probes 


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(seq))


def reconstruct_probes(probes, dp, p):
    i = len(dp) - 1
    selected_probes = []
    while i >= 0:
        if i == 0 or dp[i] > dp[i-1]:
            selected_probes.append(probes[i])
            i = p[i] - 1
        else:
            i -= 1
    return list(reversed(selected_probes))

def preprocess(probes, min_distance, overlap):
    p = [0]*len(probes)
    for i in range(len(probes)):
        for j in range(i-1, -1, -1):
            if probes[i][0] - probes[j][1] >= min_distance - overlap:
                p[i] = j + 1
                break
    return p

def max_probes(probes_list, min_distance=0, overlap=0):
    # Sort by end position
    probes_list.sort(key=lambda x: x[1])

    # Preprocess to get the list of previous non-conflicting probes
    p = preprocess(probes_list, min_distance, overlap)

    # Create an array to store the maximum number of probes for every index
    dp = [0]*len(probes_list)
    dp[0] = 1
    for i in range(1, len(probes_list)):
        # Maximum number of probes including the current probe
        include_probe = 1
        if p[i] != 0:
            include_probe += dp[p[i] - 1]

        # Maximum number of probes excluding the current probe
        exclude_probe = dp[i-1]

        # Take the maximum of these two possibilities
        dp[i] = max(include_probe, exclude_probe)

    return dp, p

def eft_probes(probes_list, min_distance=0, overlap=0):
    # Sort by end position
    probes_list.sort(key=lambda x: x[1])

    selected_probes = [probes_list[0]]
    for i in range(1, len(probes_list)):
        if probes_list[i][0] >= selected_probes[-1][1] + min_distance + 1:
            selected_probes.append(probes_list[i])

    return selected_probes

def load_and_process(filename, min_distance=0, overlap=0, method='dp'):
    data = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            genome, locus, start, end = line[:4]
            start = int(start)
            end = int(end)
            
            if locus in data:
                data[locus].append((start, end, line))
            else:
                data[locus] = [(start, end, line)]
    
    selected_probes = {}
    for locus, probes_list in data.items():
        if method == 'dp':
            dp, p = max_probes(probes_list, min_distance, overlap)
            selected_probes[locus] = reconstruct_probes(probes_list, dp, p)
        elif method == 'heuristic':
            selected_probes[locus] = eft_probes(probes_list, min_distance, overlap)
        else:
            raise ValueError(f"Invalid method: {method}")

    return selected_probes


def save_selected_probes(filename, selected_probes):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        # Write the header
        writer.writerow(['assembly', 'locus_tag', 'start', 'end', 'GC_content', 'sequence', 'reverse_complement'])
        for locus, probes in selected_probes.items():
            for probe in probes:
                # Unpack the metadata
                metadata = probe[2]
                sequence = metadata[-1] # sequence is the last item in metadata
                metadata.append(reverse_complement(sequence).lower()) # add reverse complement to the end
                writer.writerow(metadata)



def gene_duplicates(run_dir):
    df = pd.read_csv(rf'{run_dir}\duplicated_genes.txt',sep='\t')
    df.columns = ['gene','hits']

    db = {}
    for i in df.index:
        row = df.iloc[i]
        locus = row.gene.split(';')[0]
        db[locus] = ''
        for h in ast.literal_eval(row.hits):
            h_locus = h.split(';')[0]
            if db[locus] == '':
                db[locus] = h_locus
            else:
                db[locus] = db[locus] + ';' + h_locus
    return db


def save_locus_probes_count(summary_filename, annotation_filename, selected_probes,run_dir):
    
    duplicate_db = gene_duplicates(run_dir)
    
    annotation_data = {}
    with open(annotation_filename, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            line = line.strip().split('\t')
            locus = line[2]
            annotation_data[locus] = line
            
    with open(summary_filename, 'w') as file:
        # Write header
        header = ['contig', 'gene_type', 'locus_tag', 'gene_name', 'gene_fr', 'gene_to', 'gene_st', 'product', 'probe_count','duplicate_genes']
        file.write('\t'.join(header) + '\n')
        
        for locus, line in annotation_data.items():
            
            line.append(str(len(selected_probes.get(locus, []))))
            duplicates = ''
            if locus in duplicate_db:
                duplicates = duplicate_db[locus]
            line.append(duplicates)
            file.write('\t'.join(line) + '\n')

