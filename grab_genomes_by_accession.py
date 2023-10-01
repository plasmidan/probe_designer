import gzip 
from Bio import Entrez
import string
from ftplib import FTP
import gzip
import shutil
import pandas as pd
import re
import os
import glob
from Bio import SeqIO
import hashlib

def download_genome(assembly_name,genome_database_dir,download_mode = 'basic',is_unzip=True):
    
    #get the NCBI url
    url = get_genome_url(assembly_name, "your.email@example.com")

    #download the files and unzip
    genome_dir = download_ftp_dir(assembly_name,url, genome_database_dir, download_mode, unzip=is_unzip)

    #generate a simple genes annotation file from gff
    parse_gff_file(genome_dir,assembly_name)

    #merge RNA and CDS fastas and rename using locus_tag;contig notation
    merge_and_rewrite_fasta(genome_dir,assembly_name)
    
    return


def get_genome_url(assembly_name, email):
    Entrez.email = email  # Always tell NCBI who you are
    handle = Entrez.esearch(db="assembly", term=assembly_name)
    record = Entrez.read(handle)
    ids = record['IdList']

    for _id in ids:
        handle = Entrez.esummary(db="assembly", id=_id)
        summary = Entrez.read(handle)

        # Get FTP link for genome data
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        file_name = url.split('/')[-1]
        link = f"{url}"
        return link  # Return the URL for the genome


def get_organism_name(ftp, file_name):
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



def calculate_checksum(file_path):
    sha256_hash = hashlib.sha256()
    with open(file_path,"rb") as f:
        for byte_block in iter(lambda: f.read(4096),b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def download_ftp_dir(accession, url, output_dir, download_mode='basic', unzip=False):
    # List of basic files to download
    basic_files = ['assembly_stats.txt', 'cds_from_genomic.fna.gz', 'genomic.gff.gz','rna_from_genomic.fna.gz','genomic.gtf.gz'] #,

    # Split the URL into the FTP server address and the path
    url_split = url.split('/')
    server = url_split[2]
    path = '/'.join(url_split[3:])

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
        organism_name = get_organism_name(ftp, assembly_report_file)

        # Clean the organism name to make it a valid directory name
        valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
        organism_name = ''.join(c for c in organism_name if c in valid_chars)
        organism_name = organism_name.replace(' ', '_')  # Replace spaces with underscores

        organism_name = organism_name + '_' + accession
        # Create a subdirectory in the output directory with the organism name
        output_dir = os.path.join(output_dir, organism_name)

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        for file in files:
            # If download mode is 'basic' and the file is not a basic file, skip it
            if download_mode == 'basic' and not any(file.endswith(ext) for ext in basic_files):
                continue

            # If file already exists or is the assembly report, skip
            local_file = os.path.join(output_dir, file)
            local_file_unzip = re.sub('.gz','',local_file)
            if os.path.exists(local_file) or os.path.exists(local_file_unzip) or file == assembly_report_file:
                print(f" - {file} exists, skipping...")
                continue

            # Download and unzip file
            for attempt in range(3):
                try:
                    # Download file
                    print(f"Downloading {file}...")
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
# Usage:



def parse_gff_file(genome_dir,accession):
    print('Parsing gff')
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
    
    return


def merge_and_rewrite_fasta(directory,accession):
    print('merge_and_rewrite_fasta gff')
    genes_info = pd.read_csv(rf'{directory}/{accession}.genes_info.txt',sep='\t')
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

                    locus_tag = re.search(r'locus_tag=([^\]]+)',record.description).group(1)

                    # Rewrite the record id as 'contig;locus_tag'
                    record.id = f"{locus_tag};{contig};{accession}"
                    record.description = ''  # Empty the description to make the fasta header clean
                    records.append(record)
    
    output_fasta = rf'{directory}/merged_genes.fasta'

    with open(output_fasta, "w") as out_f:
        SeqIO.write(records, out_f, "fasta")
    
    return
