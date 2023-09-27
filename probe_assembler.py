import numpy as np
import pandas as pd
import re
import os

class ProbeAssembler:
    
    def __init__(self, primary_probe_path, path_to_gene_list,probe_per_gene_file,output_dir, is_amplification = True, fwd_linker='AA', primer_set = 1, num_of_flanking=2, max_probes_per_ch={'A647':25,'A550':25,'A488':29}, probes_to_remove = ['R6','R21','R36'],available_probes=range(1,160+1), reference_genes = [],negative_control_genes={},positive_control_genes={},genes_of_interest=[]):
        self.primary_probe_path = primary_probe_path
        self.path_to_gene_list = path_to_gene_list
        self.probe_per_gene_file = probe_per_gene_file
        self.output_dir = output_dir
        self.ro_df = get_ro_df(probes_to_remove=probes_to_remove)

        if is_amplification:
            self.t7_promoter_revcom = 'CCCTATAGTGAGTCGTATTA'
            if primer_set == 1:
                self.fwd_primer = 'TTTCGTCCGCGAGTGACCAG'
                self.rev_primer_revcom = 'CAACGTCCATGTCGGGATGC'
            elif primer_set == 4:
                self.fwd_primer = 'ATGCGCTGCAACTGAGACCG'
                self.rev_primer_revcom = 'TTGTGCCAGCCTTGGTCGAG'
            else:
                print(f'choose set 1 or 4... enetered set = {primer_set}... existing..') 
                exit()
            self.fwd_linker = fwd_linker
        else:
            self.t7_promoter_revcom = ''
            self.fwd_primer = ''
            self.rev_primer_revcom = ''
            self.fwd_linker = ''
        
        self.is_amplification = is_amplification
        self.num_of_flanking = num_of_flanking
        self.max_probes_per_ch = max_probes_per_ch
        self.probes_to_replace = probes_to_remove
        self.genes_of_interest = genes_of_interest
        self.reference_genes = reference_genes
        
        self.primer_set = primer_set
        self.available_probes = available_probes
        self.positive_control_genes = positive_control_genes
        self.negative_control_genes = negative_control_genes
        self.used_readouts = []
        self.set_num = 1

    def reset_primer_set(self,primer_set):
        self.primer_set = primer_set    
        if primer_set == 1:
            self.fwd_primer = 'TTTCGTCCGCGAGTGACCAG'
            self.rev_primer_revcom = 'CAACGTCCATGTCGGGATGC'
        elif primer_set == 4:
            self.fwd_primer = 'ATGCGCTGCAACTGAGACCG'
            self.rev_primer_revcom = 'TTGTGCCAGCCTTGGTCGAG'
        else:
            print(f'choose set 1 or 4... enetered set = {primer_set}... existing..') 
            exit()
 
    def read_files(self):
        # Ensure files exist
        for filepath in [self.primary_probe_path, self.path_to_gene_list, self.output_dir]:
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"{filepath} not found!")
        
        # Load files
        self.probe_per_gene = pd.read_csv(self.probe_per_gene_file,sep='\t')

        self.all_primary_df = pd.read_csv(self.primary_probe_path, sep='\t')
        if re.search('xlsx',self.path_to_gene_list):
            self.all_genes_df = pd.read_excel(self.path_to_gene_list)
        else:
            self.all_genes_df = pd.read_csv(self.path_to_gene_list,sep='\t')

        #get rows and probe counts
        ref_gene_df = self.probe_per_gene[self.probe_per_gene['locus_tag'].isin(self.reference_genes)]
        self.ref_gene_df = ref_gene_df.sort_values(by='probe_count', ascending=True)

        self.all_genes_df = self.probe_per_gene[self.probe_per_gene['locus_tag'].isin(self.all_genes_df['locus_tag'])]
        self.ref_gene_df = self.ref_gene_df.sort_values(by='probe_count', ascending=True)

        #number of genes and controls exceed the available set?
        total_num_of_targets = len(self.all_genes_df) + len(self.positive_control_genes)*3 + len(self.reference_genes)
        if len(self.negative_control_genes) > 0:
            total_num_of_targets += len(self.negative_control_genes['genes'])

        if total_num_of_targets > len(self.available_probes):
            print(f'total requested targets = {total_num_of_targets} exceeds available options ({self.available_probes})... existing')
            exit()
        else:
            self.total_num_of_targets = total_num_of_targets

    
    def process_primary_df(self):
        self.all_primary_df['readout'] = 'NA'
        self.all_primary_df['readout_seq'] = 'NA'
        self.all_primary_df['channel'] = 'NA'
        self.all_primary_df['full_probe_seq'] = 'NA'


    def assemble_positive_contols(self):
        # gene/operon marked with 3 sets of probes: one for each fluor (e.g., dnaA)
        positive_control_genes = self.positive_control_genes

        if len(positive_control_genes) == 0:
            return

        num_of_controls = len(positive_control_genes)
        probe_set = pd.DataFrame()

        for i in range(0,num_of_controls):
            
            gene_set_name = list(positive_control_genes.keys())[i]
            genes_locus_tags = self.positive_control_genes[gene_set_name]['genes']
            readouts_to_use_dict = self.positive_control_genes[gene_set_name]['readouts']
            readouts_to_use = [readouts_to_use_dict['A488'],readouts_to_use_dict['A550'],readouts_to_use_dict['A647']] # fluors reverse sorted by intensity (ensure more probes on A488)
            all_probes = self.all_primary_df.loc[np.isin(self.all_primary_df.locus_tag,genes_locus_tags)]
            all_probes.loc[:,'locus_tag'] = gene_set_name
            
            for ro_idx in readouts_to_use: #starting with R4, R5 and R6 (idx starts from 0)
                ro_idx = ro_idx - 1 + i
                ro_id,ro_seq,ro_fluor = self.ro_df.loc[ro_idx,'ro_id'],self.ro_df.loc[ro_idx,'readout_seq'],self.ro_df.loc[ro_idx,'fluor']
                probe_num = min(len(all_probes), self.max_probes_per_ch[ro_fluor])
                if probe_num < self.max_probes_per_ch[ro_fluor]:
                    print(f'*** positive set {gene_set_name} ch = {ro_fluor} has {probe_num} probes, but need {self.max_probes_per_ch[ro_fluor]}')

                selected_probes_df = select_probes(all_probes, probe_num)
                selected_probes_df['readout'] = int(re.sub('R','',ro_id))
                selected_probes_df['readout_seq'] = ro_seq
                selected_probes_df['channel'] = ro_fluor
                selected_probes_df = self.build_final_probe(selected_probes_df,ro_seq)
                all_probes = all_probes.drop(selected_probes_df.index) # remove the ones we alread selected
                probe_set = pd.concat([probe_set, selected_probes_df])
                self.used_readouts.append(ro_id)
        self.positive_control_probes = probe_set

        return

    def assemble_negative_contols(self):
        if len(self.negative_control_genes) == 0:
            return
        negative_control_genes = self.negative_control_genes['genes']
        readouts_to_use_dict = self.negative_control_genes['readouts']
        probe_set = pd.DataFrame()
        for i in range(0,3):
            ch = list(readouts_to_use_dict.keys())[i]
            ro_idx =  list(readouts_to_use_dict.values())[i] - 1
            ro_id,ro_seq,ro_fluor = self.ro_df.loc[ro_idx,'ro_id'],self.ro_df.loc[ro_idx,'readout_seq'],self.ro_df.loc[ro_idx,'fluor']

            gene_locus_tag = negative_control_genes[i]
            all_probes = self.all_primary_df.loc[self.all_primary_df.locus_tag == gene_locus_tag]
            probe_num = min(len(all_probes), self.max_probes_per_ch[ch])

            selected_probes_df = select_probes(all_probes, probe_num)
            selected_probes_df['readout'] = int(re.sub('R','',ro_id))
            selected_probes_df['readout_seq'] = ro_seq
            selected_probes_df['channel'] = ro_fluor
            selected_probes_df = self.build_final_probe(selected_probes_df,ro_seq)
            probe_set = pd.concat([probe_set, selected_probes_df])
            self.used_readouts.append(ro_id)
        
        self.negative_control_probes = probe_set
        return
        
    def prioritize_and_partition_genes(self):
        
        # Sort and prioritize
        ref_gene_part_size = len(self.ref_gene_df) // 3
        genes_of_interest_df = self.all_genes_df[self.all_genes_df['locus_tag'].isin(self.genes_of_interest)]
    
        self.all_genes_df = pd.concat([genes_of_interest_df, self.all_genes_df[~self.all_genes_df['locus_tag'].isin(self.genes_of_interest)]])
        
        rows = len(self.all_genes_df)
        partition_size = rows // 3
        
        self.channel_sets = {
            'A647': pd.concat([self.ref_gene_df.iloc[:ref_gene_part_size], self.all_genes_df.iloc[:partition_size]]),
            'A550': pd.concat([self.ref_gene_df.iloc[ref_gene_part_size:2*ref_gene_part_size], self.all_genes_df.iloc[partition_size:2*partition_size]]),
            'A488': pd.concat([self.ref_gene_df.iloc[2*ref_gene_part_size:], self.all_genes_df.iloc[2*partition_size:]])
        }
    
    def label_gene_class(self):
        self.probe_set['gene_class'] = ''
        self.probe_set.loc[np.isin(self.probe_set['locus_tag'],self.reference_genes),'gene_class'] = 'reference'

        self.probe_set.loc[np.isin(self.probe_set['locus_tag'],self.reference_genes),'gene_class'] = 'reference'

        for pos_set in self.positive_control_genes:
            self.probe_set.loc[np.isin(self.probe_set['locus_tag'],pos_set),'gene_class'] = 'positive_control'

        self.probe_set.loc[np.isin(self.probe_set['locus_tag'],self.negative_control_genes['genes']),'gene_class'] = 'negative_control'


    def assemble_probes(self):
        
        #assemble remaining probes
        probe_set = pd.DataFrame()
        try:
            probe_set = pd.concat([probe_set, self.negative_control_probes])
            probe_set = pd.concat([probe_set, self.positive_control_probes])
        except:
            pass
        for ch, df in self.channel_sets.items():
            #go over each gene
            for i in range(0,len(df)):
                locus_tag,num_of_possible_probes = df.iloc[i]['locus_tag'],df.iloc[i]['probe_count']
                ro_id,ro_fluor,ro_seq = allocate_readout(ch, self.used_readouts, self.ro_df)
                
                gene_probes_df = self.all_primary_df[self.all_primary_df['locus_tag'] == locus_tag]
                probe_num = min(num_of_possible_probes, self.max_probes_per_ch[ch])
                
                #select probes 
                selected_probes_df = select_probes(gene_probes_df, probe_num)
                selected_probes_df['readout'] = int(re.sub('R','',ro_id))
                selected_probes_df['readout_seq'] = ro_seq
                selected_probes_df['channel'] = ro_fluor

                selected_probes_df = self.build_final_probe(selected_probes_df,ro_seq)                
                probe_set = pd.concat([probe_set, selected_probes_df])
                self.used_readouts.append(ro_id)
        self.probe_set = probe_set

        return
    
    def build_final_probe(self,selected_probes_df,ro_seq):
        #assemble full probe sequence
        for idx in selected_probes_df.index:
            primary_seq = selected_probes_df.loc[idx,'reverse_complement']

            if self.num_of_flanking == 2:
                selected_probes_df.loc[idx,'full_probe_seq'] = f'{self.fwd_primer}{self.fwd_linker}{ro_seq}AA{primary_seq}AA{ro_seq}{self.t7_promoter_revcom}{self.rev_primer_revcom}'
            elif self.num_of_flanking == 4:
                selected_probes_df.loc[idx,'full_probe_seq'] = f'{self.fwd_primer}{self.fwd_linker}{ro_seq}AA{ro_seq}AA{primary_seq}AA{ro_seq}AA{ro_seq}{self.t7_promoter_revcom}{self.rev_primer_revcom}'
            else:
                print(f'only 2 or 4 flanking regions allowed... you asked for {self.num_of_flanking}...exiting')
                exit()
        return selected_probes_df

    def save_output(self):
        output_path = fr'{self.output_dir}/probes_set_{self.set_num}_primers_{self.primer_set}.txt'
        self.probe_set.to_csv(output_path,sep='\t',index=False)
        self.total_probes = len(self.probe_set)
    
    def design(self):
        self.used_readouts = [] #allow iterative applications
        self.read_files()
        self.assemble_positive_contols()
        self.assemble_negative_contols()
        self.prioritize_and_partition_genes()
        self.assemble_probes()
        self.label_gene_class()
        self.save_output()

def select_probes(df, max_probes_per_gene,is_overlap=True):
    # Convert data to a DataFrame

    # Sort by locus_tag and start position
    df = df.sort_values(by=['start'])

    # Initialize overlapping column with False values
    if is_overlap:
        df['overlapping'] = False

        # Check overlapping probes
        for i in range(1, len(df)):
            if df.iloc[i]['locus_tag'] == df.iloc[i-1]['locus_tag'] and df.iloc[i]['start'] <= df.iloc[i-1]['end']:
                df.at[i, 'overlapping'] = True

        selected_probes = []

        for tag, group in df.groupby('locus_tag'):
            non_overlapping = group[group['overlapping'] == False].sample(min(max_probes_per_gene, len(group[group['overlapping'] == False])))
            overlapping = group[group['overlapping'] == True]

            shortfall = max_probes_per_gene - len(non_overlapping)
            if shortfall > 0 and len(overlapping) > 0:
                selected_probes.append(pd.concat([non_overlapping, overlapping.sample(min(shortfall, len(overlapping)))]))
            else:
                selected_probes.append(non_overlapping)

        result = pd.concat(selected_probes)
        result = result.drop(columns='overlapping')
        result = result.sort_values(by=['start'])
    else:
        result = df.sample(min(max_probes_per_gene,len(df)))
    
    return result


def allocate_readout(channel, used_readouts, ro_df):
    # Filter the dataframe based on the channel
    channel_df = ro_df[ro_df['fluor'] == channel]
    
    # If A647, prioritize in order R1, R4, R7, etc.
    if channel == 'A647':
        channel_df = channel_df[channel_df['ro_id'].isin(['R' + str(i) for i in range(1, len(ro_df) + 1, 3)])]
    
    for index, row in channel_df.iterrows():
        if row['ro_id'] not in used_readouts:
            return row.ro_id,row.fluor,row.readout_seq
    
    return None

def get_ro_df(probes_to_remove=[]):
    data = [
    ('R1', 'A647', 'TAATCAAAGGCCGCA'),
    ('R2', 'A488', 'TGACGGTTGTAACGG'),
    ('R3', 'A550', 'GCGGTGATTGCACCA'),
    ('R4', 'A647', 'TGCTCCCGGATTACA'),
    ('R5', 'A488', 'AATGACGCGACATAG'),
    ('R6', 'A550', 'AGGGCCACGCACAAA'),
    ('R7', 'A647', 'CGCGGTGTTCAGGAT'),
    ('R8', 'A488', 'CATGAATTCTCGCCA'),
    ('R9', 'A550', 'TAGAGGTATGTCGGA'),
    ('R10', 'A647', 'GCAGGGTGAACCTAC'),
    ('R11', 'A488', 'AATTGCATCTTCGCG'),
    ('R12', 'A550', 'AGCTCGGCCTGTGAT'),
    ('R13', 'A647', 'GGGTGGAATTTCGTA'),
    ('R14', 'A488', 'GGCTTACTGAGTGAA'),
    ('R15', 'A550', 'GGGCCAGCGCAATAA'),
    ('R16', 'A647', 'CTGGTGATAACGCTA'),
    ('R17', 'A488', 'CCGCAAGTAATCTGT'),
    ('R18', 'A550', 'TATCGCAACCGACCT'),
    ('R19', 'A647', 'ATTCATCCGCAGGTC'),
    ('R20', 'A488', 'ACTCAAGGCCGGTAT'),
    ('R21', 'A550', 'GCCTTATACAGCGTA'),
    ('R22', 'A647', 'TCCGTTATGGGCGTA'),
    ('R23', 'A488', 'GTCCCGGCGTATCAA'),
    ('R24', 'A550', 'CACTCGTACAGCGTT'),
    ('R25', 'A647', 'GACTGATTGGATGCA'),
    ('R26', 'A488', 'GGCTAGCAGTCCGTC'),
    ('R27', 'A550', 'GCGCAATAAACCCTA'),
    ('R28', 'A647', 'GTCCTAGAGGTAGTT'),
    ('R29', 'A488', 'CCCTCGTACGCCTAT'),
    ('R30', 'A550', 'AATCCAGACGGGATG'),
    ('R31', 'A647', 'GCACGAGATTCACAT'),
    ('R32', 'A488', 'TCCGCCTGGATCCAA'),
    ('R33', 'A550', 'AAGCGGACGTCAGTT'),
    ('R34', 'A647', 'ACGGAAGCTCTAATC'),
    ('R35', 'A488', 'TAGAATTCGCGGCCA'),
    ('R36', 'A550', 'TGATACGCGAACTGT'),
    ('R37', 'A647', 'AATCAGGATACGGCG'),
    ('R38', 'A488', 'GGTCTCGGATATGCA'),
    ('R39', 'A550', 'TGCGTACGACAATGT'),
    ('R40', 'A647', 'ACAACGTGCGCCTGA'),
    ('R41', 'A488', 'TAAGATGTCGCCGTG'),
    ('R42', 'A550', 'TATCATACGACGTGG'),
    ('R43', 'A647', 'TGGATCATTCTACGG'),
    ('R44', 'A488', 'TCTAGTCTCACGCAA'),
    ('R45', 'A550', 'AGATTGGCGCGACAT'),
    ('R46', 'A647', 'AACGCGTAGGCTAGA'),
    ('R47', 'A488', 'TTGACGAAGCCTGAC'),
    ('R48', 'A550', 'ACTACGCGCGGATAT'),
    ('R49', 'A647', 'TTCCACTGCGTAGTG'),
    ('R50', 'A488', 'ATTGTGCGAGGTCCA'),
    ('R51', 'A550', 'TCAATTGGGCGCTCA'),
    ('R52', 'A647', 'CTCTAGGATCGTACG'),
    ('R53', 'A488', 'TGACTCGCACTGGGA'),
    ('R54', 'A550', 'TTATTACGGAGGCTG'),
    ('R55', 'A647', 'TAGCTCGACAGGGTT'),
    ('R56', 'A488', 'TGCGGACCTCTTATA'),
    ('R57', 'A550', 'GCTAGGCACCTGGAT'),
    ('R58', 'A647', 'GCGGCATTGACTAAG'),
    ('R59', 'A488', 'TCGGGCGGCTATAGA'),
    ('R60', 'A550', 'TGTTCGCCAGGATTG'),
    ('R61', 'A647', 'TGGTCGAGGTATTCG'),
    ('R62', 'A488', 'ACATTCGTTATCGCG'),
    ('R63', 'A550', 'ATATTCAGGCTAGCG'),
    ('R64', 'A647', 'ACGCTAGTAACGTCT'),
    ('R65', 'A488', 'CCTAGGCATATTCGA'),
    ('R66', 'A550', 'AGTACTCTTACGCGA'),
    ('R67', 'A647', 'TAACGAGTGCTAACG'),
    ('R68', 'A488', 'TATTCGTGAATCGCG'),
    ('R69', 'A550', 'ACGCGCGTCGTTATA'),
    ('R70', 'A647', 'TATTGGCTACGGAGC'),
    ('R71', 'A488', 'ATACGGACAATTCGG'),
    ('R72', 'A550', 'GTCGTAATAGCACGC'),
    ('R73', 'A647', 'TCTCGTGATACGGAA'),
    ('R74', 'A488', 'CTGCGAGACATTGCT'),
    ('R75', 'A550', 'CGCTAATTGCCGAGT'),
    ('R76', 'A647', 'TGCCTTGTCACGCAA'),
    ('R77', 'A488', 'AGCACCTTTACATCG'),
    ('R78', 'A550', 'CCGGTTAGCAATGCT'),
    ('R79', 'A647', 'CCACGATCGCTGATT'),
    ('R80', 'A488', 'TGGCCCATACTCGTT'),
    ('R81', 'A550', 'CTTCAGCGATAGGCA'),
    ('R82', 'A647', 'GTAAGCTATAAGCGG'),
    ('R83', 'A488', 'ATTGCGAGGCCACTA'),
    ('R84', 'A550', 'ACTTCAACGAGTGAC'),
    ('R85', 'A647', 'AGTCGTATACGCGTG'),
    ('R86', 'A488', 'CACTTAACTAGACGC'),
    ('R87', 'A550', 'CACAGACGCATATTC'),
    ('R88', 'A647', 'CGCGTCACGTTCAAT'),
    ('R89', 'A488', 'CTGATAAGCTTGGAC'),
    ('R90', 'A550', 'ATGCACGACTTCGGA'),
    ('R91', 'A647', 'GCTGTCGCCTATAAG'),
    ('R92', 'A488', 'TTAATCGGCCGAACG'),
    ('R93', 'A550', 'AATTCGACGGTCTAG'),
    ('R94', 'A647', 'GCGACTGGAACACTT'),
    ('R95', 'A488', 'TGATACGCTTGGCAT'),
    ('R96', 'A550', 'TGGCAGATCGCTAAT'),
    ('R97', 'A647', 'CGCGCAGACGTTAAT'),
    ('R98', 'A488', 'CATGATACAATCCGG'),
    ('R99', 'A550', 'GGCGTATCGAACATT'),
    ('R100', 'A647', 'ATTGCACAACCGATG'),
    ('R101', 'A488', 'CGTACTGCAGATCGA'),
    ('R102', 'A550', 'TACGACTCGGTCATT'),
    ('R103', 'A647', 'TCGTAGAGCCTAATC'),
    ('R104', 'A488', 'GGCGCACAACGAAAT'),
    ('R105', 'A550', 'TCGATGAATCGTCGA'),
    ('R106', 'A647', 'TGTGCGAACACTCGT'),
    ('R107', 'A488', 'AACGATCGCGGTGAT'),
    ('R108', 'A550', 'ATAGAATATCGGCCG'),
    ('R109', 'A647', 'TATGTGACTACGCAC'),
    ('R110', 'A488', 'CGTCGTATCCAGTAT'),
    ('R111', 'A550', 'CTATTATCGCCGAGA'),
    ('R112', 'A647', 'CGATTAGTCGTCACT'),
    ('R113', 'A488', 'ACTCCGAATGCTACG'),
    ('R114', 'A550', 'GGTTACACGCGACTA'),
    ('R115', 'A647', 'ATGTAACCAAGCGTC'),
    ('R116', 'A488', 'TCAGTTACCGGTGTA'),
    ('R117', 'A550', 'TCCAGCTTACGTTCG'),
    ('R118', 'A647', 'AGTACAACGGCATTC'),
    ('R119', 'A488', 'AACCTGCCAACGCTT'),
    ('R120', 'A550', 'ATCTAACCGGGACGT'),
    ('R121', 'A647', 'CTAGTGATTCGCACC'),
    ('R122', 'A488', 'CTGGTCATAACTCCC'),
    ('R123', 'A550', 'AATCAGGCGTGTTAC'),
    ('R124', 'A647', 'CGCTTGCGATTTCGA'),
    ('R125', 'A488', 'CTTTCGTACCTACGG'),
    ('R126', 'A550', 'ATTATCAACGAGACC'),
    ('R127', 'A647', 'CATGTACGACCGTAA'),
    ('R128', 'A488', 'AATCCTCGGTCGTAG'),
    ('R129', 'A550', 'GGGTATGATCGCAGT'),
    ('R130', 'A647', 'TGTGTTTCGCAACGA'),
    ('R131', 'A488', 'TCAGAGCTACGACTA'),
    ('R132', 'A550', 'AGTAGTCGCTAGGCT'),
    ('R133', 'A647', 'GACGAGTAACAACTG'),
    ('R134', 'A488', 'TTAAGGATCAGACCG'),
    ('R135', 'A550', 'TGACGAAGGGCTCTA'),
    ('R136', 'A647', 'CGTGATGAACCGCTT'),
    ('R137', 'A488', 'CTGTTAACGAGGCCT'),
    ('R138', 'A550', 'CGTACCGACATCGAA'),
    ('R139', 'A647', 'GAACGCCGCTATTAA'),
    ('R140', 'A488', 'GTGACTGGAGCGTTA'),
    ('R141', 'A550', 'ATAGCGTTCCACGGG'),
    ('R142', 'A647', 'AATACGCGAGCTCAG'),
    ('R143', 'A488', 'TTGGCGCCGCTAAAT'),
    ('R144', 'A550', 'CGCTTTTACGTACCA'),
    ('R145', 'A647', 'GTACAACCTCGGTTA'),
    ('R146', 'A488', 'AGACTTCGCCGTACA'),
    ('R147', 'A550', 'TGATTTGGAGGACGA'),
    ('R148', 'A647', 'CGTCGATAGGACTCA'),
    ('R149', 'A488', 'ATCTTGCGCGTTGGG'),
    ('R150', 'A550', 'TTCAGAGCCGGTAGA'),
    ('R151', 'A647', 'CACCGGTCTACGTAA'),
    ('R152', 'A488', 'AAAGTTCGCGACCTA'),
    ('R153', 'A550', 'TTGCTCGTCCTACGA'),
    ('R154', 'A647', 'AAGAACGGGTCACTC'),
    ('R155', 'A488', 'ATGCGTCGGTGCTAT'),
    ('R156', 'A550', 'GCTTGGGTCGAGCAA'),
    ('R157', 'A647', 'GTTACGGAACAGCGA'),
    ('R158', 'A488', 'GCGTGATCGGTGATT'),
    ('R159', 'A550', 'CGCCGAGCGATATAA'),
    ('R160', 'A647', 'TGCGTTGACGCAATA'),
    ('R161', 'A488', 'GTGTGCAGAATAGCG'),
    ('R162', 'A550', 'ACATTCGCAAGGAAC'),
    ('R163', 'A647', 'AGGATTCCAAGGCGA'),
    ('R164', 'A488', 'AACTGAAGCAACGGT'),
    ('R165', 'A550', 'AGTATCAGCTTGCAC'),
    ('R166', 'A647', 'CAGCCATAGATCTAG'),
    ('R167', 'A488', 'CCTCAAGTTATGCGC'),
    ('R168', 'A550', 'TCCATATACGAGGTG'),
    ('R169', 'A647', 'GGAGTAATACGACGC'),
    ('R170', 'A488', 'CCGTAAAGTCGACTG'),
    ('R171', 'A550', 'AACTATACCGTCCAA'),
    ('R172', 'A647', 'TATTGGCTGTTACGA'),
    ('R173', 'A488', 'TCTTTACCGACTCAT'),
    ('R174', 'A550', 'GCTGTTCGCGTGATA'),
    ('R175', 'A647', 'CCCGTAGATGACTTT'),
    ('R176', 'A488', 'ATTATCGGTTCCAGC'),
    ('R177', 'A550', 'TCGGCTTTTCGCATG'),
    ('R178', 'A647', 'GAAGTCAAATCCGTC'),
    ('R179', 'A488', 'CCCAATCGCTTTGGC'),
    ('R180', 'A550', 'ACATCCCGCGTATAT'),
    ('R181', 'A647', 'TTCAATCCGTCGACA'),
    ('R182', 'A488', 'GTACCCGTGTGGTAA'),
    ('R183', 'A550', 'CAAGGTGATTCGTGT'),
    ('R184', 'A647', 'TGAGTAAAGACCGCC'),
    ('R185', 'A488', 'AAACCGCTATATGCC'),
    ('R186', 'A550', 'TCTGCGTCGAGAGTT'),
    ('R187', 'A647', 'TACGAGTGTACCGAA'),
    ('R188', 'A488', 'CCCTTTACGACTACC'),
    ('R189', 'A550', 'GAGCGTTTGATTGAC'),
    ('R190', 'A647', 'GTCCCGCCTTATAAA'),
    ('R191', 'A488', 'CTAGCACTGGTCGTC'),
    ('R192', 'A550', 'ACATAACATCCCGGT'),
    ('R193', 'A647', 'TTGGGGACGTATCGA'),
    ('R194', 'A488', 'TAGACCGCCTATCAT'),
    ('R195', 'A550', 'TCTTCAACTTGCTCG'),
    ('R196', 'A647', 'TTGAACATTGCCGCC'),
    ('R197', 'A488', 'GCGAATACCATACGC'),
    ('R198', 'A550', 'CTATAGACGAATAGC'),
    ('R199', 'A647', 'CAAATCCGCGTGTTA'),
    ('R200', 'A488', 'GCTCTAAGTGTCCGA'),
    ('R201', 'A550', 'CGCCGTGACGTATTT'),
    ('R202', 'A647', 'AAGGTGTGTGAACGC'),
    ('R203', 'A488', 'TCTAGCGGACATGAC'),
    ('R204', 'A550', 'CCAGTTGCCGTTACA'),
    ('R205', 'A647', 'GTCGATTGATAAGGC'),
    ('R206', 'A488', 'TTATCATGCCGCCAC'),
    ('R207', 'A550', 'AAGACGTGCTCCAGA'),
    ('R208', 'A647', 'AAGTGAGGATTACGC'),
    ('R209', 'A488', 'TATTAACGTGCCTGC'),
    ('R210', 'A550', 'AGCGTCATCGAAATC'),
    ('R211', 'A647', 'AACTTCCCGACTCCG'),
    ('R212', 'A488', 'TTATGCCGGAGCAAA'),
    ('R213', 'A550', 'GATACTAGCAGCCGC'),
    ('R214', 'A647', 'ATCGCAATCGTACCG'),
    ('R215', 'A488', 'CTACAATCGAGCGAA'),
    ('R216', 'A550', 'TTGGCAGTCGCACAG'),
    ('R217', 'A647', 'CAATGGATTACCTCG'),
    ('R218', 'A488', 'TAGTCGCTCGAATAC'),
    ('R219', 'A550', 'GTCCGTATTGCGCTA'),
    ('R220', 'A647', 'CGCTGATCGGATAAT'),
    ('R221', 'A488', 'AAGAAGGCGAGACAA'),
    ('R222', 'A550', 'TGCAACCCCGCAACT'),
    ('R223', 'A647', 'AGATAACGATACGCA'),
    ('R224', 'A488', 'TCCTGATGTTACTGC'),
    ('R225', 'A550', 'CCGGACCGTGAAACA'),
    ('R226', 'A647', 'CCAGGTTCGAACTAC'),
    ('R227', 'A488', 'AGTTGGTCCGCGATA'),
    ('R228', 'A550', 'GGAAACTGATCGCGT'),
    ('R229', 'A647', 'GCCAAATAGCGTGCA'),
    ('R230', 'A488', 'CTATGTTCCTTCGGC'),
    ('R231', 'A550', 'GGGTTCATGACTAGA'),
    ('R232', 'A647', 'ATACGGCGGAAGATA'),
    ('R233', 'A488', 'TTGTCCCGCGTTGAA'),
    ('R234', 'A550', 'CACAAGGGCATCGTT'),
    ('R235', 'A647', 'AGGTAACCCATCGGA'),
    ('R236', 'A488', 'CCTGTATGACCGTAT'),
    ('R237', 'A550', 'GCGTCGCGTTAGACT'),
    ('R238', 'A647', 'AACGTCGCCTAAACT'),
    ('R239', 'A488', 'AGTGTTAGGCGTATG'),
    ('R240', 'A550', 'TCTGGTTTATCCGTG'),
]

    df = pd.DataFrame(data, columns=['ro_id', 'fluor', 'readout_seq'])
    # Filter out rows with 'R3' or 'R5'
    df_filtered = df[~df['ro_id'].isin(probes_to_remove)]
    # df = df_filtered.reset_index(drop=True)
            

    # return df
    return df_filtered