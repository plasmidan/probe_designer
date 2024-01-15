import gzip 
from Bio import Entrez
import re
import string
from ftplib import FTP
import os
import gzip
import shutil
import pandas as pd
import os
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def downloadNcbiGenomes(accessions: dict, genome_database_dir: str, download_mode: str = 'basic', is_unzip: bool = True,
                   is_refseq_format: bool = True) -> None:
    """Downloading and unzipping (if needed) the file desired genome file.
    Creating from the gff file convenient annotation, and merging RNA and CDS.
    """

    os.makedirs(genome_database_dir, exist_ok=True)
    
    for org_name,accession in accessions.items():
        print(f'downloading: {org_name} -> {accession}')
        # get the NCBI url
        url = getGenomeUrl(accession, "your.email@example.com", is_refseq_format)

        # download the files and unzip
        genome_dir = downloadFtpDir(accession, url, genome_database_dir, download_mode, unzip=is_unzip)

        # generate a simple genes annotation file from gff
        parseGffFile(genome_dir, accession)

        # merge RNA and CDS fastas and rename using locus_tag;contig notation
        mergeAndRewriteFasta(genome_dir, accession)

    return None


def getGenomeUrl(assembly_name: str, email: str, is_refseq_format: bool = True) -> str:
    """Connect to Entrez system and return a download link.
    """
    Entrez.email = email  # Always tell NCBI who you are
    handle = Entrez.esearch(db="assembly", term=assembly_name)
    record = Entrez.read(handle)
    ids = record['IdList']

    for _id in ids:
        handle = Entrez.esummary(db="assembly", id=_id)
        summary = Entrez.read(handle)

        # Get FTP link for genome data
        if is_refseq_format:
            url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        else: 
            url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        if url == '':
            continue
        file_name = url.split('/')[-1]
        return f"{url}"  # Return the URL for the genome


def getOrganismName(ftp: FTP, file_name: str) -> str:
    """Read organism name from file, simplify and return it.
    """
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith("# Organism name:"):
                # Get the organism name and clean it
                organism_name = line.split(":")[1].strip()
                # Remove everything inside parentheses
                organism_name = re.sub(r'\(.*?\)', '', organism_name)
                # Remove leading and trailing spaces
                organism_name = organism_name.strip()
                # Remove all periods
                organism_name = organism_name.replace('.', '')
                # Convert only the first character to lowercase
                organism_name = organism_name[0].lower() + organism_name[1:]
                # Remove trailing underscores
                organism_name = organism_name.rstrip('_')
                return organism_name


def downloadFtpDir(accession: str, url: str, output_dir: str, download_mode: str = 'basic', unzip: bool = False) -> str:
    """Connect to server, then login as anonymous user.
    Download and unzip basic files if they are not exist.
    """
    # List of basic files to download
    basic_files = ['assembly_stats.txt', 'cds_from_genomic.fna.gz', 'genomic.gff.gz', 'rna_from_genomic.fna.gz']

    # Split the URL into the FTP server address and the path
    url_split = url.split('/')
    server = url_split[2]
    path = '/'.join(url_split[3:])
    del url_split
    with FTP(server) as ftp:
        ftp.login()  # Log in as anonymous user

        # Change to the directory on the FTP server
        ftp.cwd(path)

        # List all files in the directory
        files = ftp.nlst()

        # Download assembly report first to get organism name
        assembly_report_file = [file for file in files if file.endswith('assembly_stats.txt')][0]
        with open(assembly_report_file, 'wb') as f:
            ftp.retrbinary('RETR ' + assembly_report_file, f.write)

        # Get the organism name
        organism_name = getOrganismName(ftp, assembly_report_file)

        # Clean the organism name to make it a valid directory name
        valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
        organism_name = ''.join(c for c in organism_name if c in valid_chars)
        organism_name = organism_name.replace(' ', '_')  # Replace spaces with underscores

        output_name = organism_name + '_' + accession
        # Create a subdirectory in the output directory with the organism name
        output_dir = os.path.join(output_dir, output_name)
        os.makedirs(output_dir, exist_ok=True)


        for file in files:
            # If download mode is 'basic' and the file is not a basic file, skip it
            if download_mode == 'basic' and not any(file.endswith(ext) for ext in basic_files):
                continue

            # If file already exists or is the assembly report, skip
            local_file = os.path.join(output_dir, file)
            local_file_unzip = re.sub('.gz', '', local_file)
            if os.path.exists(local_file) or os.path.exists(local_file_unzip) or file == assembly_report_file:
                # print(f" - {file} exists, skipping...")
                continue

            # Download and unzip file
            for attempt in range(3):
                try:
                    # Download file
                    # print(f"Downloading {file}...")
                    with open(local_file, 'wb') as f:
                        ftp.retrbinary('RETR ' + file, f.write)

                    # Unzip file if necessary
                    if unzip and file.endswith('.gz'):
                        try:
                            with gzip.open(local_file, 'rb') as f_in:
                                with open(local_file[:-3], 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                        except Exception as e:
                            print(f"Exception occurred while decompressing {file}: {str(e)}")
                        os.remove(local_file)

                    # If the file was successfully downloaded and decompressed, break the loop
                    break
                except OSError:
                    print(f"Error occurred while downloading/unzipping {file}. Attempt {attempt + 1} of 3 failed.")
                    
                    # Delete the corrupted file if it exists
                    if os.path.exists(local_file):
                        os.remove(local_file)

                # If all three attempts failed, print an error message
                if attempt == 2:
                    print(f"Could not download/unzip {file} after 3 attempts.")
    return output_dir


def parseGffFile(genome_dir: str, accession: str) -> None:
    """Create more readable version for the genes info."""
    # Initialize empty list to hold the data
    data = []
    gff_file = glob.glob(os.path.join(genome_dir, f"*{'.gff'}"))
    with open(rf'{gff_file[0]}', 'r') as f:

        while True:
            line1 = f.readline().rstrip()
            if line1 == '': break

            line_1_sp = line1.rstrip().split('\t')
            try:
                if line_1_sp[2] == 'gene':
                    line2 = f.readline().rstrip()
                    line_2_sp = line2.rstrip().split('\t')
                    contig = line_2_sp[0]
                    gene_type = line_2_sp[2]
                    gene_fr,gene_to,gene_st = line_2_sp[3],line_2_sp[4],line_2_sp[6]
                    locus_tag = re.search(r'locus_tag=([^;]+)',line_2_sp[8]).group(1) #;product=(.+)
                    product = re.sub(' ','_',re.search(r'product=([^;]+)',line_2_sp[8]).group(1))
                    try:
                        gene_name = re.sub(' ','_',re.search(r'gene=([^;]+)',line_2_sp[8]).group(1))
                    except:
                        gene_name = ''
                    # append the data as a dictionary to the list
                    data.append({'contig': contig, 'gene_type': gene_type,'locus_tag': locus_tag,'gene_name': gene_name, 'gene_fr': gene_fr, 'gene_to': gene_to, 'gene_st': gene_st, 'product': product})
            except:
                pass

    # Convert list of dictionaries to pandas dataframe
    df = pd.DataFrame(data)

    # Write to csv
    df.to_csv(rf'{genome_dir}/{accession}.genes_info.txt', sep='\t', index=False)

    #prepare a file to guide exon analysis for euks
    cds_exon_info_path = os.path.join(genome_dir, 'cds_exon_info.txt')
    extractCdsInfo(df,rf'{gff_file[0]}', cds_exon_info_path)
    return None


def extractCdsInfo(genes_info_df,gff_file_path, output_file_path):
    with open(gff_file_path, 'r') as gff_file, open(output_file_path, 'w') as output_file:
        #gene_name = ''
        locus_tag = ''
        cds_count = 0
        output_file.write(f'locus_tag\tstart_position\tend_position\tcds_count\n')
        for line in gff_file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')

                if len(fields) >= 8:
                    feature_type = fields[2]
                    attributes = dict(item.split('=') for item in fields[8].split(';'))

                    if feature_type == 'CDS':
                        #current_gene_name = attributes.get('Name', '')
                        current_locus_tag = attributes.get('locus_tag', '')
                        if not genes_info_df['locus_tag'].isin([current_locus_tag]).any(): 
                            continue
                        if current_locus_tag != locus_tag:
                            locus_tag = current_locus_tag
                            cds_count = 1
                        else:
                            cds_count += 1

                        #cds_identifier = attributes.get('protein_id', '')
                        start_position, end_position = int(fields[3]), int(fields[4])

                        # Write to the output file
                        output_file.write(f'{locus_tag}\t{start_position}\t{end_position}\t{cds_count}\n')


def mergeAndRewriteFasta(directory: str, accession: str) -> None:
    """ Combine coding sequence and RNA to one fasta file."""
    genes_info_path = os.path.join(directory, accession + '.genes_info.txt')
    genes_info = pd.read_csv(genes_info_path, sep='\t')
    contig = genes_info.iloc[0]['contig']

    records = []
    
    cds_fasta = glob.glob(os.path.join(directory, f"*{'_cds_from_genomic.fna'}"))
    rna_fasta = glob.glob(os.path.join(directory, f"*{'_rna_from_genomic.fna'}"))

    for file_list in [cds_fasta, rna_fasta]:
        for file in file_list:
            with open(file, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    header_parts = record.description.split(" ")
                    # extract contig and locus_tag from header
                    try:
                        locus_tag = re.search(r'locus_tag=([^\]]+)', record.description).group(1)
                        
                        #make sure locus tag is in GFF - avoid pseudogenes
                        if not genes_info['locus_tag'].isin([locus_tag]).any(): 
                            continue

                    except (Exception,):
                        print('error check genomic seq file for corruption...')
                        exit()
                    # Rewrite the record id as 'contig;locus_tag'
                    record.id = f"{locus_tag};{contig};{accession}"
                    record.description = ''  # Empty the description to make the fasta header clean
                    records.append(record)
    
    output_fasta = rf'{directory}/merged_genes.fasta'
    with open(output_fasta, "w") as out_f:
        for record in records:
            # Create a new record with the sequence as a single string
            new_record = SeqRecord(Seq(str(record.seq)), id=record.id, description='')
            # Write the modified record to the file
            out_f.write(new_record.format("fasta-2line"))

    return None


#add parsers for Integrated Microbial Genomes)
