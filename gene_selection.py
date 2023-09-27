import pandas as pd
import numpy as np

class sample_genes_for_par2fish:

    """
    Script for selecting genes - dividing the transcriptome into sets
    Allows to sample transcriptional units and to oversample them with a desired frequency
    """

    def __init__(self,probe_dir,probe_sum_table,set_size,num_of_sets,min_probe_per_gene,desired_probe_count,output_dir,genes_of_interest=[],genes_to_avoid=[],max_double_operons=0,biocyc_tu_table=None,allow_duplicate_genes=False):
        self.probe_dir = probe_dir
        self.probe_sum_table = probe_sum_table
        self.biocyc_tu_table = biocyc_tu_table
        self.genes_df = pd.read_csv(rf'{probe_dir}/{probe_sum_table}',sep='\t')

        #allow_duplicate_genes - allow selecting more then one of the duplicates? 
        # ok only if the design step treated them as different genes (probes would be specific to each copy)
        self.allow_duplicate_genes = allow_duplicate_genes

        #filter genes 
        cds_only = self.genes_df[self.genes_df.gene_type == 'CDS']
        selected_only = self.genes_df[np.isin(self.genes_df.locus_tag,genes_of_interest)]
        self.genes_df = pd.concat([selected_only,cds_only])

        self.set_size = set_size
        self.num_of_sets = num_of_sets
        self.min_probe_per_gene = min_probe_per_gene
        self.genes_of_interest = genes_of_interest
        self.genes_to_avoid = genes_to_avoid
        self.max_double_operons = max_double_operons
        self.desired_probe_count = desired_probe_count
        self.output_dir = output_dir
        
        if biocyc_tu_table is not None: #operon data, needs to be manually downloaded from biocyc
            self.add_biocyc_tu_info()

        #do the sampling
        self.sample_genes()

        self.save_gene_sets()


    def add_biocyc_tu_info(self):
        
        if self.biocyc_tu_table is None:
            print('no TU table supplied...')
            return

        #load table
        tu_df = pd.read_csv(rf'{self.probe_dir}/{self.biocyc_tu_table}',sep='\t')

        # Parsing genes and converting them into a frozenset for unique comparison
        tu_df['gene_list'] = tu_df['Genes of transcription unit'].str.split(' // ')

        # Parsing genes and converting them into a frozenset for unique comparison
        tu_df['gene_set'] = tu_df['gene_list'].apply(frozenset)

        # Drop duplicates based on the gene set
        tu_df = tu_df.drop_duplicates(subset='gene_set')

        # Convert DataFrame to dictionary
        tu_dict = dict(zip(tu_df['Object ID'], tu_df['gene_list']))

        # Calculate the average number of genes per TU
        
        gene_to_tu = {gene: tu for tu, genes in tu_dict.items() for gene in genes}

        # Create a dictionary to map genes to their TU ID
        # Mapping the gene to its TU ID using the gene_to_tu dictionary
        self.genes_df['TU'] = self.genes_df['gene_name'].map(gene_to_tu)

        # Mapping the TU ID to its size (number of genes) using the tu_dict
        self.genes_df['TU_size'] = self.genes_df['TU'].map(lambda tu: len(tu_dict[tu]) if tu in tu_dict else np.nan)

        # Update the 'TU' column for rows with NaN values
        nan_mask = self.genes_df['TU'].isna()
        num_nans = nan_mask.sum()
        unique_ids = ["tu_holder_" + str(i) for i in range(1, num_nans + 1)]
        self.genes_df.loc[nan_mask, 'TU'] = unique_ids
        self.genes_df.loc[nan_mask, 'TU_size'] = 1

        self.tu_dict = tu_dict
        self.gene_to_tu = gene_to_tu
    

    def sample_genes(self):
        set_size = self.set_size
        num_of_sets = self.num_of_sets
        min_probe_per_gene = self.min_probe_per_gene
        desired_probe_count = self.desired_probe_count
        genes_of_interest = self.genes_of_interest
        max_double_operons = self.max_double_operons
        """
        set_size: how many genes per group
        num_of_sets: how many sets to select
        min_probe_per_gene: min probes per selected gene
        desired_probe_count: if TU, then grab a gene with at least this many
        genes_of_interest: genes that must be selected in the first lib

        """
        filtered_genes = self.genes_df[self.genes_df['probe_count'] >= min_probe_per_gene].copy()
        filtered_genes = filtered_genes[~np.isin(filtered_genes.locus_tag,self.genes_to_avoid)]
        
        self.filtered_genes = filtered_genes
        df_shuffled = filtered_genes.sample(frac=1).reset_index(drop=True) #randomly suffle the gene rows
        for locus_tag in genes_of_interest: #set the genes of interest as first so they will be selected in the first library
            row_to_move = df_shuffled[df_shuffled['locus_tag'] == locus_tag]
            df_shuffled = df_shuffled.drop(row_to_move.index)
            df_shuffled = pd.concat([row_to_move, df_shuffled], ignore_index=True)
        df_shuffled= df_shuffled.drop_duplicates()
        selected_genes_df = pd.DataFrame(columns=df_shuffled.columns)
  
        #TU info missing
        if self.biocyc_tu_table is None:
            for set_i in range(0,num_of_sets):
                subset = df_shuffled.iloc[set_i*set_size:set_i*set_size+set_size].copy()
                subset['probe_set'] = set_i+1
                selected_genes_df = pd.concat([selected_genes_df,subset])
            self.selected_genes_df = selected_genes_df
            self.df_shuffled = df_shuffled
        #TU info added
        else:
            #params
            current_set_num = 1
            operons_seen = {}
            double_operon_count = 0            
            duplicate_genes_seen = []
            while num_of_sets > 0:
                #select another row in the shuffled list
                row = df_shuffled.iloc[0]
                tu_rows = df_shuffled[df_shuffled.TU == row.TU] #get the rows of TU members
                tu = tu_rows.TU.iloc[0] #tu name
                if not self.allow_duplicate_genes and row.locus_tag in duplicate_genes_seen:
                    df_shuffled = df_shuffled.iloc[1:]
                    continue
        
                #found a duplicate? record it so we can avoid labeling another member             
                if not pd.isna(row.duplicate_genes):
                    for duplicate in row.duplicate_genes.split(';'):
                        duplicate_genes_seen.append(duplicate)

  
                #allow sampling the same operon again twice at a predefined rates
                if not np.isin(row.locus_tag,genes_of_interest):
                    if not tu in operons_seen: # first time for this TU
                        operons_seen[tu] = 1
                    elif np.isin(list(tu_rows.locus_tag.values),genes_of_interest).any():
                        operons_seen[tu] += 1
                    else:
                        if max_double_operons == 0 or double_operon_count >= max_double_operons or operons_seen[tu] > 1:
                            df_shuffled = df_shuffled.iloc[1:]
                            continue
                        else:
                            operons_seen[tu] += 1
                            double_operon_count += 1
                    
                
                #select the gene in the TU
                # if a specific locus was specified get it - assuming no more then one (or how many are allowed for double dipping)
                if np.isin(row.locus_tag,genes_of_interest):
                    selected_row = df_shuffled.iloc[[0]].copy()
                    selected_row['probe_set'] = current_set_num
                    df_shuffled = df_shuffled.iloc[1:]
                    
                    selected_genes_df = pd.concat([selected_genes_df,selected_row])
  
                else: #not an important gene here, we can select the gene in the TU by num of probes
                    wishful_tu_rows = tu_rows[tu_rows.probe_count >= desired_probe_count] #are there genes with a lot of probes in the operon?
                    if len(wishful_tu_rows) > 0: 
                        selected_row = wishful_tu_rows.sample(1) #random
                    else: 
                        selected_row = tu_rows[tu_rows.probe_count == max(tu_rows.probe_count)] #take the gene with most probes
                    selected_row = selected_row.copy()
                    selected_row['probe_set'] = current_set_num

                    if selected_row.locus_tag.iloc[0] in selected_genes_df.locus_tag.values:
                        df_shuffled = df_shuffled.drop(selected_row.index[0])
                    else:
                        selected_genes_df = pd.concat([selected_genes_df,selected_row])

                if len(selected_genes_df[selected_genes_df.probe_set == current_set_num]) == set_size:
                    current_set_num += 1
                    num_of_sets-=1
            self.selected_genes_df = selected_genes_df
    
    def save_gene_sets(self):
        for set_i in range(1,self.num_of_sets+1):
            gene_set = self.selected_genes_df[self.selected_genes_df.probe_set == set_i]
            gene_set.to_csv(rf'{self.output_dir}/gene_set_{set_i}_{self.set_size}_genes.txt',sep='\t',index=False)