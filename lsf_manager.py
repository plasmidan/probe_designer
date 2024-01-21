import sys,os,re
import subprocess


class lsf_manager:

    def __init__(self,run_dir: str,lsf_scripts_dir = None, queue='new-all.q',sleep_duration_s = 10):
        
        #check that the script dir was supplied
        if lsf_scripts_dir is None:
            print(f'lsf_scripts_dir not defined = {lsf_scripts_dir}\n exiting...\n')
            exit() 
        
        #make the lsf dir
        lsf_dir = fr'{run_dir}/lsf'
        os.makedirs(lsf_dir,exist_ok=True)

        #save params
        self.lsf_dir = lsf_dir
        self.lsf_scripts_dir =lsf_scripts_dir
        self.queue = queue
        self.sleep_duration_s = sleep_duration_s
        return
        

    def wait_for_all_jobs_to_finish(self,message=''):
        #try before sleep
        out, err = subprocess.Popen('bjobs', stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()
        if out == b'No unfinished job found\n': return
        # print(f'waiting for all {message} jobs to finish...\n')
        while True:
            out, err = subprocess.Popen('bjobs', stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()
            if out == b'No unfinished job found\n':
                # print('done\n')
                break
            else:
                out, err = subprocess.Popen('sleep %s'%self.sleep_duration_s, stdout=subprocess.PIPE,shell=True).communicate()


    def identifyDuplicateGenes(self,run_dir, min_pid, min_len_similarity_percent,min_percent_cov,chunk_idx ,memory=4000,num_cores=1):
        #build command
        lsf_dir_path = fr'{run_dir}/lsf'
        job_name = f'chunk_{chunk_idx}_duplicates'
        output = rf'{lsf_dir_path}/{job_name}.o'
        error = rf'{lsf_dir_path}/{job_name}.e'
        code = rf'python {self.lsf_scripts_dir}/identifyDuplicateGenes_lsf.py {run_dir} {min_pid} {min_len_similarity_percent} {min_percent_cov} {chunk_idx}'
        #execute
        bsub_cmd = rf'bsub -J {job_name} -q {self.queue} -R "rusage[mem={memory}]"  -o {output} -e {error} "{code}"' #-n {num_cores}
        bsub_cmd_out,err = subprocess.Popen(bsub_cmd, stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True).communicate()

    def generateNaiveProbes(self,run_dir, chunk_size,chunk_idx,probe_len,gc_min,gc_max,max_base_rep,memory=4000,num_cores=1):
        #build command
        lsf_dir_path = fr'{run_dir}/lsf'
        job_name = f'chunk_{chunk_idx}_naive'
        output = rf'{lsf_dir_path}/{job_name}.o'
        error = rf'{lsf_dir_path}/{job_name}.e'
        code = rf'python {self.lsf_scripts_dir}/generateNaiveProbes_lsf.py {run_dir} {chunk_size} {chunk_idx} {probe_len} {gc_min} {gc_max} {max_base_rep}'
        #execute
        bsub_cmd = rf'bsub -J {job_name} -q {self.queue} -R "rusage[mem={memory}]"  -o {output} -e {error} "{code}"' #-n {num_cores}
        bsub_cmd_out,err = subprocess.Popen(bsub_cmd, stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True).communicate()

    def blastNaiveProbes(self,run_dir,word_size,chunk_name ,memory=4000,num_cores=1):
        #build command
        lsf_dir_path = fr'{run_dir}/lsf'
        job_name = f'{chunk_name}_blast'
        output = rf'{lsf_dir_path}/{job_name}.o'
        error = rf'{lsf_dir_path}/{job_name}.e'
        code = rf'python {self.lsf_scripts_dir}/blastNaiveProbes_lsf.py {run_dir} {word_size} {chunk_name}'
        #execute
        bsub_cmd = rf'bsub -J {job_name} -q {self.queue} -R "rusage[mem={memory}]"  -o {output} -e {error} "{code}"' #-n {num_cores}
        bsub_cmd_out,err = subprocess.Popen(bsub_cmd, stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True).communicate()

    def identifySpecificProbes(self,run_dir,chunk_name,max_nonspecific_match,is_allow_gene_duplicates,memory=4000,num_cores=1):
        #build command
        lsf_dir_path = fr'{run_dir}/lsf'
        job_name = f'{chunk_name}_specific'
        output = rf'{lsf_dir_path}/{job_name}.o'
        error = rf'{lsf_dir_path}/{job_name}.e'
        code = rf'python {self.lsf_scripts_dir}/get_specific_probes_lsf.py {run_dir} {chunk_name} {max_nonspecific_match} {is_allow_gene_duplicates}'
        #execute
        bsub_cmd = rf'bsub -J {job_name} -q {self.queue} -R "rusage[mem={memory}]"  -o {output} -e {error} "{code}"' #-n {num_cores}
        bsub_cmd_out,err = subprocess.Popen(bsub_cmd, stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True).communicate()



