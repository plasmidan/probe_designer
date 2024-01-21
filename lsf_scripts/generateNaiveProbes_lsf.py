import os,sys
from probe_designer import readFastaFile,generateNaiveProbes

self_name,run_dir,chunk_size,chunk_idx,probe_len,gc_min,gc_max,max_base_rep= sys.argv

chunk_idx = int(chunk_idx)
chunk_size = int(chunk_size)
probe_len = int(probe_len)
gc_min = int(gc_min)
gc_max = int(gc_max)
max_base_rep = int(max_base_rep)

probe_props = {'probe_len': probe_len,'gc_min': gc_min,'gc_max': gc_max,'max_base_rep': max_base_rep}

#get the right chunk
fna_db = readFastaFile(rf'{run_dir}/merged_genes.fna')  
fasta_headers = list(fna_db.keys())
header_chunk_list = [fasta_headers[i:i + chunk_size] for i in range(0, len(fasta_headers), chunk_size)]

genes_in_chunk = header_chunk_list[chunk_idx]

generateNaiveProbes(run_dir=run_dir, probe_props = probe_props ,genes_in_chunk = genes_in_chunk,chunk_idx=chunk_idx)
