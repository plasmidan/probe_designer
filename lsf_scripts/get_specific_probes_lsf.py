import os,sys
from probe_designer import identifySpecificProbes

self_name,run_dir,chunk_name,max_nonspecific_match,is_allow_gene_duplicates= sys.argv

is_allow_gene_duplicates = is_allow_gene_duplicates == 'True'
max_nonspecific_match = int(max_nonspecific_match)

identifySpecificProbes(run_dir=run_dir, max_nonspecific_match=max_nonspecific_match, is_allow_gene_duplicates=is_allow_gene_duplicates, is_parallel=True, chunk_name=chunk_name)
