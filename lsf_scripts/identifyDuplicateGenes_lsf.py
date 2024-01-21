import os,sys
from probe_designer import identifyDuplicateGenes

self_name,run_dir, min_pid, min_len_similarity_percent,min_percent_cov,chunk_idx= sys.argv

min_pid = int(min_pid)
min_len_similarity_percent = int(min_len_similarity_percent)
min_percent_cov = int(min_percent_cov)
chunk_idx = int(chunk_idx)

duplicated_genes_db = identifyDuplicateGenes(run_dir, min_pid=min_pid, min_len_similarity_percent=min_len_similarity_percent, min_percent_cov=min_percent_cov,chunk_idx=chunk_idx)

