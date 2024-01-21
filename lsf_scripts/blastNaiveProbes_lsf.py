import os,sys
import subprocess
from probe_designer import blastNaiveProbes

self_name,run_dir,word_size,chunk_name = sys.argv

word_size = int(word_size)

current_dir = os.getcwd()
os.chdir(run_dir)

query = rf'./lsf/{chunk_name}'
db = 'merged_genes'
output_path = rf'./lsf/{chunk_name}.blast_res'
cmd = f'blastn -query {query} -db {db} -word_size {word_size} -evalue 1000 -outfmt 6 -out {output_path} -num_threads 4'
subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
os.chdir(current_dir)