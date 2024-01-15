import os
from pathlib import Path
import subprocess
import re
import csv
import ast
import pandas as pd


def findSpecificProbes(genome_database_dir: str, probes_dir: str, accessions: list, run_name: str, 
                 probe_props: {'probe_len': 30,'gc_min': 40,'gc_max': 65,'max_base_rep': 4},
                 max_nonspecific_match: int = 18, is_allow_gene_duplicates: bool = True,
                 min_duplicate_pid: int = 95, min_duplicate_len_similarity_per: int = 90, min_duplicate_percent_cov: int = 90) -> None:
    
    """
    genome_database_dir: where the downloaded genome directories are kept
    probes_dir: where to save the results
    accessions: the list of org+genomes to use in the design (can be multiple entries); accessions = {'ecoli_k12' : 'ASM584v2'}
    run_name: the name of the output result directory
    probe_props: dictionary with paramters. Default -> probe_props = {'probe_len': 30,'gc_min': 40, 'gc_max': 65, 'max_base_rep': 4}
    max_nonspecific_match: the maximum number of continuous perfect matches with an off target gene


    is_allow_gene_duplicates: how to handle gene duplicates.
     - True:  allow nonspecific binding between gene duplicates (same probe can bind multiple genes)
     - False: only probes that distinguish duplicates
    *** for example, to also get probes that binds all 16S rRNA copies set the paramter to True.

    Putative gene duplicates are identified by blasting gene DNA sequences (all vs. all). 
    Two or more genes are termed duplicates using 3 tunable paramers:
     - the percent identity (min_duplicate_pid = 95%)
     - the percent coverage (min_duplicate_percent_cov = 90%)
     - gene similarity in length (nt) (min_duplicate_len_similarity_per=90%)
    """

    #open a new dir?
    print('Designing probes...')
    run_dir = rf'{probes_dir}/{run_name}'
    os.makedirs(probes_dir, exist_ok=True)
    os.makedirs(run_dir, exist_ok=True)
    
    fna_db = collectGenomesForProbeDesign(genome_database_dir, accessions, run_dir)
    
    runMakeBlastDb(fna_db, 'nucl', 'merged_genes')

    print(' - generate naive probes')
    generateNaiveProbes(run_dir, probe_props)
    
    print(' - blasting naive probes')
    blastNaiveProbesBlast(run_dir, word_size=max_nonspecific_match)

    duplicated_genes_db = identifyDuplicateGenes(run_dir, min_pid=min_duplicate_pid, min_len_similarity_percent=min_duplicate_len_similarity_per, min_percent_cov=min_duplicate_percent_cov)
    
    print(' - determining probe specificities')
    identifySpecificProbes(run_dir, duplicated_genes_db, max_nonspecific_match=max_nonspecific_match,is_allow_gene_duplicates=is_allow_gene_duplicates)

    return None


def findMaxProbesPerGene(run_name: str, accessions: dict, genome_database_dir: str, probes_dir: str, min_probe_distance: int = 0, is_exon_only: bool = False) -> None:
    """Calculate maximum probes possible per gene, under the constraints.
    heuristic-greedy (sort and select)
    """
    #report the distance between probes
    if min_probe_distance < 0:
        print(f'Probe overlap enabled for: {min_probe_distance} bp')
    
    for org_name,accession in accessions.items():

        print(org_name,accession)
        selected_probes = loadAndProcess(
            genome_database_dir = genome_database_dir,
            specific_probes_file = rf'{probes_dir}/{run_name}/all_specific_probes.txt',
            accession = accession, # work on this accession only
            min_probe_distance = min_probe_distance,
            is_exon_only = is_exon_only
            )
        
        if is_exon_only:
            exons_str = '.exons'
        else:
            exons_str = ''
        selected_probe_output_path = rf'{probes_dir}/{run_name}/{org_name}{exons_str}.max_probes.txt'
        reads_per_gene_path = f'{org_name}{exons_str}.probes_per_gene.txt'
        
        saveSelectedProbes(
            output_filename = selected_probe_output_path,
            selected_probes = selected_probes,
            )
        
        saveLocusProbesCount(
            run_dir = rf'{probes_dir}/{run_name}',
            accession = accession,
            org_name = org_name,
            selected_probes = selected_probes,
            output_name = reads_per_gene_path
            )

    mergeProbesPerGene(run_dir_path = rf'{probes_dir}/{run_name}',acessions = accessions,is_exon_only= is_exon_only)

    return None



def mergeProbesPerGene(run_dir_path: str,acessions:dict,is_exon_only: bool) -> None:
    """in case of multispecies also produce a merged file for later analysis"""
    
    if len(acessions) < 2: # no need to merge
        return
    else:
        print('Merging outputs')

    if is_exon_only:
        exons_str = '.exons'
    else:
        exons_str = ''

    acc_counter = 1
    with open(f'{run_dir_path}/merged{exons_str}.max_probes.txt','w') as mp_o:
        with open(f'{run_dir_path}/merged{exons_str}.probes_per_gene.txt','w') as ppg_o:
            for org_name,acc in acessions.items():
                with open(rf'{run_dir_path}/{org_name}{exons_str}.max_probes.txt','r') as mp_fh:
                    header = mp_fh.readline()
                    if acc_counter == 1:
                        mp_o.write(header)
                    for line in mp_fh:
                        mp_o.write(line)
                with open(rf'{run_dir_path}/{org_name}{exons_str}.probes_per_gene.txt','r') as ppg_fh:
                    header = ppg_fh.readline()
                    if acc_counter == 1:
                        ppg_o.write(header)
                    for line in ppg_fh:
                        ppg_o.write(line)
                acc_counter = 0
    return

def runMakeBlastDb(input_file: str, db_type: str, output_file: str) -> None:
    """
    running the makeblastdb cmd in the directory.
    """
    current_dir = os.getcwd()

    os.chdir(Path(input_file).parent)
    
    cmd = f'makeblastdb -in "{Path(input_file).name}" -dbtype {db_type} -out "{output_file}"'
    subprocess.run(cmd, shell=True, capture_output=True, text=True)

    os.chdir(current_dir)
    return


def collectGenomesForProbeDesign(genome_database_dir: str, assembly_dict: dict, run_dir: str) -> str:
    """ Grab the genomic data and merge it if there is more than one organism
    """
    downloaded_genomes_list = os.listdir(genome_database_dir)
    output_fna_path = rf'{run_dir}/merged_genes.fna'
    header_count = 1
    with open(output_fna_path, 'w') as o:
        with open(rf'{run_dir}/merged_gene_info.txt', 'w') as genes_o:
            for genome in downloaded_genomes_list:
                for org,accession in assembly_dict.items():
                    if genome.endswith(accession):
                        with open(rf'{genome_database_dir}/{genome}/merged_genes.fasta', 'r') as f:
                            while(True):
                                h = f.readline()
                                s = f.readline()
                                if not h:
                                    break       
                                o.write(h)
                                o.write(s)
                                    
                        with open(rf'{genome_database_dir}/{genome}/{accession}.genes_info.txt', 'r') as f:
                            header = f.readline()
                            if header_count > 0:
                                genes_o.write('accession' + '\t' + header)
                                header_count = 0
                            for line in f:
                                genes_o.write(accession + '\t' + line)
    return output_fna_path


def readFastaFile(fasta_file: str) -> dict:
    """Opens and reads fasta file to dictionary line by line."""
    db = {}
    current_header = ''
    with open(fasta_file, 'r') as fa:
        for line in fa:
            if re.search('^>', line):
                current_header = line.rstrip()
                db[current_header] = ''
            else:
                db[current_header] += line.rstrip()
    return db


def calculateProbeGCContent(probe_seq: str) -> float:
    """Calculate the percentage of g or c bases."""
    at, gc, n = 0, 0, 0
    for i in probe_seq:
        if i in 'AT':
            at += 1
        elif i in 'GC':
            gc += 1
        elif i == 'N':
            n += 1
    
    if at + gc == 0:
        return 0
    elif n > 0: #remove N-containing sequences.
        return 0
    probe_gc_content = round(gc / (at + gc) * 100, 1)
    return probe_gc_content


def generateNaiveProbes(run_dir: str, probe_props: dict):
    """Produce all possible probes with sliding window approach."""

    fna_db = readFastaFile(rf'{run_dir}/merged_genes.fna')
    naive_probes_path = rf'{run_dir}/naive_probes.fna'

    probe_len, gc_min, gc_max, max_base_rep = int(probe_props['probe_len']), int(probe_props['gc_min']),\
        int(probe_props['gc_max']), int(probe_props['max_base_rep'])

    with open(naive_probes_path, 'w') as probes_out_fh:
        for header, gene_seq in fna_db.items():
            header = re.sub('>', '', header)
            locus_tag, contig, assembly = header.split(';')

            for i in range(len(gene_seq)-probe_len):
                # run over gene sequence, base by base and test each potential probe
                sub_seq = str(gene_seq[i:i+probe_len])
                sub_seq_gc = calculateProbeGCContent(sub_seq)
                if gc_max >= sub_seq_gc >= gc_min:
                    max_rep_tuple = (max_base_rep, max_base_rep, max_base_rep, max_base_rep)
                    if not re.search("A{%s}|G{%s}|C{%s}|T{%s}" % max_rep_tuple, sub_seq):
                        probe_fr, probe_to = str(i+1), str(i+probe_len)
                        probe_header = f'{locus_tag};{contig};{assembly};{sub_seq_gc};{probe_fr}-{probe_to}'
                        probes_out_fh.write(f'>{probe_header}\n{sub_seq}\n')
    del fna_db
    return None


def blastNaiveProbesBlast(run_dir: str, word_size: int = 18, evalue: int = 100, n_threads: int = 4,output_format: int = 6) -> None:
    """Run nucleotide blast on the sliding window naive probes."""
    current_dir = os.getcwd()
    os.chdir(run_dir)

    query = './naive_probes.fna'
    db = './merged_genes'
    output_path = './naive_probes.blast_results'

    cmd = f'blastn -query {query} -db {db} -word_size {word_size} -evalue {evalue} -outfmt {output_format}' \
          f' -out {output_path} -num_threads {n_threads}'

    subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

    os.chdir(current_dir)
    return None


def identifyDuplicateGenes(run_dir: str, min_pid: int = 95, min_len_similarity_percent: int = 90,min_percent_cov: int = 90) -> dict:
    """Find the potentially duplicated gene
        - similar length and nearly identical sequence across the majority of the genes
        - cross species duplicates are not allowed (always treated as different genes)
    """
    current_dir = os.getcwd()
    os.chdir(run_dir)
    fasta_db = readFastaFile(rf'{run_dir}/merged_genes.fna')
    self_blast_path = r'gene_duplicate_identification.blast_res'
    query = r'merged_genes.fna'
    db = r'merged_genes'
    cmd = f'blastn -query {query} -db {db} -word_size 28 -evalue 100 -outfmt 6 -out {self_blast_path}'
    subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

    duplicated_genes_db = {}
    with open(self_blast_path, 'r') as blast_res_fh:
        for line in blast_res_fh:
            query_gene_header, hit_gene_header, pid, match_len, mismatches, gap_len, q_start, q_end, h_start, h_end,evalue, bitscore = line.rstrip().split('\t')

            # same gene hit
            if query_gene_header == hit_gene_header:
                continue

            q_locus, q_contig, q_assembly = query_gene_header.split(';')
            # org_id,locus_tag,gene,desc
            h_locus, h_contig, h_assembly = hit_gene_header.split(';')

            # don't allow cross species
            if not q_assembly == h_assembly:
                continue
            if q_locus == h_locus:
                continue

            # Is the percent identity significant + is the gene size similar enough?
            is_duplicate = False
            if float(pid) >= min_pid:
                q_gene_len = len(fasta_db['>%s' % query_gene_header])
                h_gene_len = len(fasta_db['>%s' % hit_gene_header])
                gene_len_ratio = min(q_gene_len, h_gene_len) / max(q_gene_len, h_gene_len) * 100
                # percent cov of the smallest gene
                match_percent_cov = float(match_len) / min(q_gene_len, h_gene_len) * 100
                # genes of similar length?
                if gene_len_ratio >= min_len_similarity_percent and match_percent_cov >= min_percent_cov:
                    is_duplicate = True

            if is_duplicate:
                if query_gene_header not in duplicated_genes_db:
                    duplicated_genes_db[query_gene_header] = [hit_gene_header]
                else:
                    duplicated_genes_db[query_gene_header].append(hit_gene_header)

    with open('./duplicated_genes.txt', 'w') as f:
        for k, v in duplicated_genes_db.items():
            f.write(f'{k}\t{v}\n')
    os.chdir(current_dir)
    return duplicated_genes_db


def identifySpecificProbes(run_dir: str, duplicated_genes_db: dict, max_nonspecific_match: int = 18,is_allow_gene_duplicates: bool = True) -> None:
    """Keeps specific probes from all naive sliding window probes."""
    
    naive_probe_fasta_db = readFastaFile(rf'{run_dir}/naive_probes.fna')
    probe_result_db = {}
    with open(rf'{run_dir}/naive_probes.blast_results', 'r') as f:
        for line in f:  # BSU_00010;NC_000964.3;ASM904v1;40.0;6-35
            probe_header, hit_header, pid, match_len, mismatches, gap_len, q_start, q_end, h_start, h_end, evalue, bitscore = line.rstrip().split('\t')
            probe_locus, probe_contig, probe_assembly, probe_gc, probe_fr_to = probe_header.split(';')

            probe_fr, probe_to = probe_fr_to.split('-')

            hit_locus, hit_contig, hit_assembly = hit_header.split(';')
            is_duplicate_hit = False

            probe_header_split = probe_header.split(';') 
            probe_header_trim = ';'.join(probe_header_split[:3])

            if hit_header in duplicated_genes_db and probe_header_trim in duplicated_genes_db[hit_header]:
                is_duplicate_hit = True

            # override in case this option is set to False
            if not is_allow_gene_duplicates:
                is_duplicate_hit = False

            # register the probe if not in the db
            if probe_header not in probe_result_db:
                probe_result_db[probe_header] = {
                    'probe_seq': naive_probe_fasta_db[f'>{probe_header}'],
                    'number_of_non_specific_hits': 0,
                    'all_non_specific_hits': [],
                    'probe_assembly': probe_assembly,
                    'probe_locus': probe_locus,
                    'probe_gc': probe_gc,
                    'probe_fr': probe_fr,
                    'probe_to': probe_to,
                    'probe_desc': 'NA'
                }
            # is specific hit?
            if int(match_len) >= max_nonspecific_match:
                if not probe_locus == hit_locus and not is_duplicate_hit:  # hit on another gene ok if it's your double
                    probe_result_db[probe_header]['number_of_non_specific_hits'] += 1
                    probe_result_db[probe_header]['all_non_specific_hits'].append(hit_header)
    
    # save the specific probes
    with open(fr'{run_dir}/all_specific_probes.txt', 'w') as o:
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
                line = '\t'.join([probe_assembly, probe_locus, probe_fr, probe_to, probe_gc, probe_seq])
                o.write(line + '\n')
    return None


def reverseComplement(seq: str) -> str:
    """Create the reverse complement for the input DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement[base] for base in reversed(seq))


def grabMaxProbes(probes_list: list, min_probe_distance: int = 0) -> list:
    """Heuristic-greedy algorithm for finding large number of probes possible from a list
    under the overlap constraints
    """
    probes_list.sort(key=lambda x: x[1])

    selected_probes = [probes_list[0]]
    for i in range(1, len(probes_list)):
        if probes_list[i][0] >= selected_probes[-1][1] + min_probe_distance: # + 1
            selected_probes.append(probes_list[i])

    return selected_probes


def loadAndProcess(genome_database_dir: str, specific_probes_file: str, accession: str, min_probe_distance: int = 0,is_exon_only: bool = False) -> dict:
    """Loads and reads the specific probes table, and creates data dictionary for every locus.
    Then runs the algorithms that was chosen by the user to find large number of probes to produce.
    """
    specific_probes = {}
    with open(specific_probes_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            probe_accession, locus, start, end = line[:4]
            if not probe_accession == accession:
                continue
            start = int(start)
            end = int(end)
            
            if locus in specific_probes:
                specific_probes[locus].append((start, end, line))
            else:
                specific_probes[locus] = [(start, end, line)]
    
    if is_exon_only:
        genome_dir_name = [genome for genome in os.listdir(genome_database_dir) if genome.endswith(accession)][0]
        cds_regions = createCdsRegionsDict(rf'{genome_database_dir}/{genome_dir_name}/cds_exon_info.txt') 
        gene_info_df = pd.read_csv(rf'{genome_database_dir}/{genome_dir_name}/{accession}.genes_info.txt',sep='\t')        
        specific_probes = grabExonProbes(gene_info_df,specific_probes , cds_regions)

    selected_probes = {}
    for locus, probes_list in specific_probes.items():
        selected_probes[locus] = grabMaxProbes(probes_list, min_probe_distance)

    return selected_probes


def saveSelectedProbes(output_filename: str, selected_probes: dict) -> None:
    
    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        # Write the header
        writer.writerow(['assembly', 'locus_tag', 'start', 'end', 'GC_content', 'sequence', 'reverse_complement'])
        for locus, probes in selected_probes.items():
            for probe in probes:
                # Unpack the metadata
                metadata = probe[2]
                sequence = metadata[5]  # sequence is one before the last item in metadata
                #sequence = metadata[-2]
                metadata.append(reverseComplement(sequence).lower())  # add reverse complement to the end
                writer.writerow(metadata)
    return None


def getGeneDuplicates(run_dir: str) -> dict:
    """Handles gene duplicates."""
    try:
        df = pd.read_csv(rf'{run_dir}\duplicated_genes.txt', sep='\t')
        df.columns = ['gene','hits']
    except:
        df = pd.DataFrame(columns=['gene', 'hits'])
        
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


def saveLocusProbesCount(run_dir: str, accession: str, org_name: str,selected_probes: dict, output_name: str) -> None:
    duplicate_db = getGeneDuplicates(run_dir)
    annotation_data = {}
    with open(rf'{run_dir}/merged_gene_info.txt', 'r') as file:
        next(file) 
        for line in file:
            line = line.strip().split('\t')
            gene_accession = line[0]
            if gene_accession == accession:
                locus = line[3]
                annotation_data[locus] = line
            
    with open(rf'{run_dir}/{output_name}', 'w') as file: 
        header = ['contig', 'gene_type', 'locus_tag', 'gene_name', 'gene_fr', 'gene_to', 'gene_st', 'product','probe_count', 'duplicate_genes']
        file.write('\t'.join(header) + '\n')
        
        for locus, line in annotation_data.items():
            line.append(str(len(selected_probes.get(locus, []))))
            duplicates = ''
            if locus in duplicate_db:
                duplicates = duplicate_db[locus]
            line.append(duplicates)
            file.write('\t'.join(line) + '\n')
    return None


def grabExonProbes(gene_info_df, probe_dict, exon_data):
    filtered_probes = {}
    for locus_tag, probes in probe_dict.items():
        if locus_tag in exon_data:
            exons = exon_data[locus_tag]
            locus_start_pos = int(gene_info_df.loc[gene_info_df.locus_tag == locus_tag, 'gene_fr'].iloc[0])

            for probe in probes:
                start, end = probe[0], probe[1]
                norm_start = start + locus_start_pos - 1
                norm_end = end + locus_start_pos - 1

                # Check each exon for a match with the current probe
                for exon_start, exon_end in exons:
                    if exon_start <= norm_start <= exon_end and exon_start <= norm_end <= exon_end:
                        if locus_tag not in filtered_probes:
                            filtered_probes[locus_tag] = []
                        filtered_probes[locus_tag].append(probe)
                        break  # Break out of the inner loop if a match is found
        else:
            for probe in probes: #non spliced genes
                if locus_tag not in filtered_probes:
                            filtered_probes[locus_tag] = []
                filtered_probes[locus_tag].append(probe)

    return filtered_probes


def createCdsRegionsDict(cds_file_path):
    # Load the data into a DataFrame
    cds_df = pd.read_csv(cds_file_path, sep='\t')

    # Sort the DataFrame by locus_tag and start_position
    cds_df.sort_values(by=['locus_tag', 'start_position'], inplace=True)

    # Create a dictionary of sorted CDS regions
    cds_regions = {}
    for row in cds_df.itertuples(index=False):
        if row.locus_tag not in cds_regions:
            cds_regions[row.locus_tag] = []
        cds_regions[row.locus_tag].append((row.start_position, row.end_position))

    return cds_regions