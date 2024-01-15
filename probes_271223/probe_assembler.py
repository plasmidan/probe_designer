import numpy as np
import pandas as pd
import re
import os


class ProbeAssembler:
    """Class for combining the readouts, promoters and linkers to the primary probe sequence."""
    def __init__(self, run_dir: str, primary_probes: str, selected_gene_list: str, probe_per_gene: str,
                 is_amplification: bool = True, num_of_flanking: int = 2, fwd_primer: str = 'TTTCGTCCGCGAGTGACCAG', rev_primer_revcom: str = 'CAACGTCCATGTCGGGATGC',t7_promoter_revcom: str = 'CCCTATAGTGAGTCGTATTA',
                 max_probes_per_ch: dict = {'fluor_640nm': 25, 'fluor_560nm': 25, 'fluor_488nm': 29},
                 probes_to_remove: list = ['R6', 'R21', 'R36'], reference_genes=[],
                 negative_control_genes = {}, positive_control_genes = {}, genes_of_interest = [] ) -> None:
        self.run_dir = run_dir
        self.selected_gene_list = selected_gene_list
        self.primary_probe_path = rf'{run_dir}/{primary_probes}'
        self.path_to_gene_list = rf'{run_dir}/{selected_gene_list}'
        self.probe_per_gene_file = rf'{run_dir}/{probe_per_gene}'
        self.run_dir = run_dir
        self.ro_df = getRoDf(probes_to_remove=probes_to_remove)

        self.is_amplification = is_amplification
        self.num_of_flanking = num_of_flanking
        self.fwd_primer = fwd_primer
        self.rev_primer_revcom = rev_primer_revcom
        self.t7_promoter_revcom = t7_promoter_revcom
        self.max_probes_per_ch = max_probes_per_ch
        self.probes_to_replace = probes_to_remove # 
        self.genes_of_interest = genes_of_interest
        self.reference_genes = reference_genes

        self.available_probes = len(self.ro_df) #
        self.positive_control_genes = positive_control_genes
        self.negative_control_genes = negative_control_genes
        self.used_readouts = []

    def readFiles(self) -> None:
        """Check existence and get probes from files."""
        # Ensure files exist
        for filepath in [self.primary_probe_path, self.path_to_gene_list, self.run_dir]:
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"{filepath} not found!")

        # Load files
        self.probe_per_gene = pd.read_csv(self.probe_per_gene_file, sep='\t')

        #load all primary probes
        self.all_primary_df = pd.read_csv(self.primary_probe_path, sep='\t')
        self.all_primary_df['readout'] = 'NA'
        self.all_primary_df['readout_seq'] = 'NA'
        self.all_primary_df['channel'] = 'NA'
        self.all_primary_df['full_probe_seq'] = 'NA'
        
        if re.search('xlsx', self.path_to_gene_list):
            self.all_genes_df = pd.read_excel(self.path_to_gene_list)
        else:
            self.all_genes_df = pd.read_csv(self.path_to_gene_list, sep='\t')

        # get rows and probe counts
        ref_gene_df = self.probe_per_gene[self.probe_per_gene['locus_tag'].isin(self.reference_genes)]
        self.ref_gene_df = ref_gene_df.sort_values(by='probe_count', ascending=True)

        self.all_genes_df = self.probe_per_gene[self.probe_per_gene['locus_tag'].isin(self.all_genes_df['locus_tag'])]
        self.ref_gene_df = self.ref_gene_df.sort_values(by='probe_count', ascending=True)

        # number of genes and controls exceed the available set?
        total_num_of_targets = len(self.all_genes_df) + len(self.positive_control_genes)*3 + len(self.reference_genes)
        if len(self.negative_control_genes) > 0:
            total_num_of_targets += len(self.negative_control_genes['genes'])

        if total_num_of_targets > self.available_probes:
            print(f'total requested targets = {total_num_of_targets} exceeds available options ({self.available_probes}'
                  f'). existing')
            exit()
        else:
            self.total_num_of_targets = total_num_of_targets
            return None

    def assemblePositiveControls(self) -> None:
        """Define positive control genes.
        Gene/operon marked with 3 sets of probes: one for each fluor (e.g., dnaA)
        """
        positive_control_genes = self.positive_control_genes

        if len(positive_control_genes) == 0:
            return None

        num_of_controls = len(positive_control_genes)
        probe_set = pd.DataFrame()

        for i in range(num_of_controls):

            gene_set_name = list(positive_control_genes.keys())[i]
            genes_locus_tags = self.positive_control_genes[gene_set_name]['genes']
            readouts_to_use_dict = self.positive_control_genes[gene_set_name]['readouts']
            readouts_to_use = [readouts_to_use_dict['fluor_488nm'], readouts_to_use_dict['fluor_560nm'], readouts_to_use_dict['fluor_640nm']]
            # fluors reverse sorted by intensity (ensure more probes on fluor_488nm)
            all_probes = self.all_primary_df.loc[np.isin(self.all_primary_df.locus_tag, genes_locus_tags)]
            all_probes.loc[:, 'locus_tag'] = gene_set_name


            for ro_idx in readouts_to_use:
                # starting with R4, R5 and R6 (idx starts from 0)
                ro_idx = ro_idx - 1 + i
                ro_id, ro_seq, ro_fluor = self.ro_df.loc[ro_idx, 'ro_id'], self.ro_df.loc[ro_idx, 'readout_seq'],\
                    self.ro_df.loc[ro_idx, 'fluor']
                probe_num = min(len(all_probes), self.max_probes_per_ch[ro_fluor])
                if probe_num < self.max_probes_per_ch[ro_fluor]:
                    print(f'*** positive set {gene_set_name} ch = {ro_fluor} has {probe_num} probes, but need'
                          f' {self.max_probes_per_ch[ro_fluor]}')

                selected_probes_df = selectProbes(all_probes, probe_num)
                selected_probes_df['readout'] = int(re.sub('R', '', ro_id))
                selected_probes_df['readout_seq'] = ro_seq
                selected_probes_df['channel'] = ro_fluor
                selected_probes_df = self.buildFinalProbe(selected_probes_df,ro_seq)
                all_probes = all_probes.drop(selected_probes_df.index)  # remove the ones we already selected
                probe_set = pd.concat([probe_set, selected_probes_df])
                self.used_readouts.append(ro_id)
        self.positive_control_probes = probe_set
        return None

    def assembleNegativeControls(self) -> None:
        """Define negative control genes."""
        if len(self.negative_control_genes) == 0:
            return None
        negative_control_genes = self.negative_control_genes['genes']
        readouts_to_use_dict = self.negative_control_genes['readouts']
        probe_set = pd.DataFrame()
        for i in range(3):
            ch = list(readouts_to_use_dict.keys())[i]
            ro_idx = list(readouts_to_use_dict.values())[i] - 1
            ro_id, ro_seq, ro_fluor = self.ro_df.loc[ro_idx, 'ro_id'], self.ro_df.loc[ro_idx, 'readout_seq'],\
                self.ro_df.loc[ro_idx, 'fluor']

            gene_locus_tag = negative_control_genes[i]
            all_probes = self.all_primary_df.loc[self.all_primary_df.locus_tag == gene_locus_tag]
            probe_num = min(len(all_probes), self.max_probes_per_ch[ch])

            selected_probes_df = selectProbes(all_probes, probe_num)
            selected_probes_df['readout'] = int(re.sub('R', '', ro_id))
            selected_probes_df['readout_seq'] = ro_seq
            selected_probes_df['channel'] = ro_fluor
            selected_probes_df = self.buildFinalProbe(selected_probes_df, ro_seq)
            probe_set = pd.concat([probe_set, selected_probes_df])
            self.used_readouts.append(ro_id)

        self.negative_control_probes = probe_set
        return None

    def prioritizeAndPartitionGenes(self) -> None:
        """Prioritize and set channels."""
        ref_gene_part_size = len(self.ref_gene_df) // 3
        genes_of_interest_df = self.all_genes_df[self.all_genes_df['locus_tag'].isin(self.genes_of_interest)]

        self.all_genes_df = pd.concat([genes_of_interest_df, self.all_genes_df[~self.all_genes_df['locus_tag'].isin(
            self.genes_of_interest)]])

        rows = len(self.all_genes_df)
        partition_size = rows // 3

        self.channel_sets = {
            'fluor_640nm': pd.concat([self.ref_gene_df.iloc[:ref_gene_part_size], self.all_genes_df.iloc[:partition_size]]),
            'fluor_560nm': pd.concat([self.ref_gene_df.iloc[ref_gene_part_size:2*ref_gene_part_size], self.all_genes_df.iloc[partition_size:2*partition_size]]),
            'fluor_488nm': pd.concat([self.ref_gene_df.iloc[2*ref_gene_part_size:], self.all_genes_df.iloc[2*partition_size:]])
        }
        return None

    def labelGeneClass(self) -> None:
        """Assign gene class to genes in probe_set."""
        self.probe_set['gene_class'] = ''
        self.probe_set.loc[np.isin(self.probe_set['locus_tag'], self.reference_genes), 'gene_class'] = 'reference'

        self.probe_set.loc[np.isin(self.probe_set['locus_tag'], self.reference_genes), 'gene_class'] = 'reference'

        for pos_set in self.positive_control_genes:
            self.probe_set.loc[np.isin(self.probe_set['locus_tag'], pos_set), 'gene_class'] = 'positive_control'

        self.probe_set.loc[np.isin(self.probe_set['locus_tag'], self.negative_control_genes['genes']),
                           'gene_class'] = 'negative_control'
        return None

    def assembleProbes(self) -> None:
        """Assemble remaining probes."""
        
        probe_set = pd.DataFrame()
        try:
            probe_set = pd.concat([probe_set, self.negative_control_probes])
        except (Exception,):
            pass
        try:
            probe_set = pd.concat([probe_set, self.positive_control_probes])
        except (Exception,):
            pass
        for ch, df in self.channel_sets.items():
            # go over each gene
            for i in range(len(df)):
                locus_tag, num_of_possible_probes = df.iloc[i]['locus_tag'], df.iloc[i]['probe_count']
                ro_id, ro_fluor, ro_seq = allocateReadout(ch, self.used_readouts, self.ro_df)

                gene_probes_df = self.all_primary_df[self.all_primary_df['locus_tag'] == locus_tag].copy()                
                probe_num = min(num_of_possible_probes, self.max_probes_per_ch[ch])

                # select probes
                selected_probes_df = selectProbes(gene_probes_df, probe_num)
                selected_probes_df['readout'] = int(re.sub('R', '', ro_id))
                selected_probes_df['readout_seq'] = ro_seq
                selected_probes_df['channel'] = ro_fluor
                # selected_probes_df['full_probe_seq'] = 'NA'
                selected_probes_df = self.buildFinalProbe(selected_probes_df, ro_seq)
                probe_set = pd.concat([probe_set, selected_probes_df])
                self.used_readouts.append(ro_id)
        self.probe_set = probe_set
        return None

    def buildFinalProbe(self, selected_probes_df: pd, ro_seq: str) -> pd:
        """Assemble full probe sequence."""
        for idx in selected_probes_df.index:
            primary_seq = selected_probes_df.loc[idx, 'reverse_complement']

            if self.num_of_flanking == 2:
                if self.is_amplification:
                    selected_probes_df.loc[idx, 'full_probe_seq'] = f'{self.fwd_primer}AA{ro_seq}AA{primary_seq}AA{ro_seq}{self.t7_promoter_revcom}{self.rev_primer_revcom}'
                else:
                    selected_probes_df.loc[idx, 'full_probe_seq'] = f'{ro_seq}AA{primary_seq}AA{ro_seq}'
            elif self.num_of_flanking == 4:
                if self.is_amplification:
                    selected_probes_df.loc[idx, 'full_probe_seq'] = f'{self.fwd_primer}AA{ro_seq}AA{ro_seq}AA{primary_seq}AA{ro_seq}AA{ro_seq}{self.t7_promoter_revcom}{self.rev_primer_revcom}'
                else:
                   selected_probes_df.loc[idx, 'full_probe_seq'] = f'{ro_seq}AA{ro_seq}AA{primary_seq}AA{ro_seq}AA{ro_seq}' 
            else:
                print(f'only 2 or 4 flanking regions allowed, you asked for {self.num_of_flanking}. exiting')
                exit()
        return selected_probes_df

    def saveOutput(self) -> None:
        """Save selected probes to txt file."""
        output_name = re.sub('.txt','',self.selected_gene_list)
        output_path = fr'{self.run_dir}/{output_name}.assembled_probes.txt'
        self.probe_set.to_csv(output_path, sep='\t', index=False)
        self.total_probes = len(self.probe_set)
        return None

    def design(self) -> None:
        """Main method."""
        self.used_readouts = []  # allow iterative applications
        self.readFiles()
        self.assemblePositiveControls()
        self.assembleNegativeControls()
        self.prioritizeAndPartitionGenes()
        self.assembleProbes()
        try:
            self.labelGeneClass()
        except(Exception,):
            pass
        self.saveOutput()
        return None


def selectProbes(df: pd, max_probes_per_gene: int, is_overlap: bool = True) -> pd:
    """Filter overlapping probes"""
    df = df.sort_values(by=['start'])

    # Initialize overlapping column with False values
    if is_overlap:
        df['overlapping'] = False

        # Check overlapping probes
        for i in range(1, len(df)):
            # check if they are overlapping' if yes, flag them as overlapping
            same_locus = (df.iloc[i]['locus_tag'] == df.iloc[i-1]['locus_tag'])
            overlap_pos = (df.iloc[i]['start'] <= df.iloc[i-1]['end'])
            if same_locus and overlap_pos:
                df.iat[i, df.columns.get_loc('overlapping')] = True

        selected_probes = []

        for tag, group in df.groupby('locus_tag'):
            non_overlapping = group[group['overlapping'] == False].sample(min(
                                                    max_probes_per_gene, len(group[group['overlapping'] == False])))
            overlapping = group[group['overlapping'] == True]

            shortfall = max_probes_per_gene - len(non_overlapping)
            if shortfall > 0 and len(overlapping) > 0:
                selected_probes.append(pd.concat([non_overlapping, overlapping.sample(
                                                    min(shortfall, len(overlapping)))]))
            else:
                selected_probes.append(non_overlapping)

        result = pd.concat(selected_probes)
        result = result.drop(columns='overlapping')
        result = result.sort_values(by=['start'])
    else:
        result = df.sample(min(max_probes_per_gene, len(df)))

    return result


def allocateReadout(channel: str, used_readouts: list, ro_df: pd) -> tuple[str, str, str] or None:
    """Filter the dataframe based on the channel.
    If fluor_640nm, prioritize in order R1, R4, R7, etc.
    Return readout specificities.
    """
    channel_df = ro_df[ro_df['fluor'] == channel]

    if channel == 'fluor_640nm':
        channel_df = channel_df[channel_df['ro_id'].isin(['R' + str(i) for i in range(1, len(ro_df) + 1, 3)])]

    for index, row in channel_df.iterrows():
        if row['ro_id'] not in used_readouts:
            return row.ro_id, row.fluor, row.readout_seq

    return None


def getRoDf(probes_to_remove: list = []) -> pd:
    """Removes problematic readouts."""
    data = [
        ('R1', 'fluor_640nm', 'TAATCAAAGGCCGCA'),
        ('R2', 'fluor_488nm', 'TGACGGTTGTAACGG'),
        ('R3', 'fluor_560nm', 'GCGGTGATTGCACCA'),
        ('R4', 'fluor_640nm', 'TGCTCCCGGATTACA'),
        ('R5', 'fluor_488nm', 'AATGACGCGACATAG'),
        ('R6', 'fluor_560nm', 'AGGGCCACGCACAAA'),
        ('R7', 'fluor_640nm', 'CGCGGTGTTCAGGAT'),
        ('R8', 'fluor_488nm', 'CATGAATTCTCGCCA'),
        ('R9', 'fluor_560nm', 'TAGAGGTATGTCGGA'),
        ('R10', 'fluor_640nm', 'GCAGGGTGAACCTAC'),
        ('R11', 'fluor_488nm', 'AATTGCATCTTCGCG'),
        ('R12', 'fluor_560nm', 'AGCTCGGCCTGTGAT'),
        ('R13', 'fluor_640nm', 'GGGTGGAATTTCGTA'),
        ('R14', 'fluor_488nm', 'GGCTTACTGAGTGAA'),
        ('R15', 'fluor_560nm', 'GGGCCAGCGCAATAA'),
        ('R16', 'fluor_640nm', 'CTGGTGATAACGCTA'),
        ('R17', 'fluor_488nm', 'CCGCAAGTAATCTGT'),
        ('R18', 'fluor_560nm', 'TATCGCAACCGACCT'),
        ('R19', 'fluor_640nm', 'ATTCATCCGCAGGTC'),
        ('R20', 'fluor_488nm', 'ACTCAAGGCCGGTAT'),
        ('R21', 'fluor_560nm', 'GCCTTATACAGCGTA'),
        ('R22', 'fluor_640nm', 'TCCGTTATGGGCGTA'),
        ('R23', 'fluor_488nm', 'GTCCCGGCGTATCAA'),
        ('R24', 'fluor_560nm', 'CACTCGTACAGCGTT'),
        ('R25', 'fluor_640nm', 'GACTGATTGGATGCA'),
        ('R26', 'fluor_488nm', 'GGCTAGCAGTCCGTC'),
        ('R27', 'fluor_560nm', 'GCGCAATAAACCCTA'),
        ('R28', 'fluor_640nm', 'GTCCTAGAGGTAGTT'),
        ('R29', 'fluor_488nm', 'CCCTCGTACGCCTAT'),
        ('R30', 'fluor_560nm', 'AATCCAGACGGGATG'),
        ('R31', 'fluor_640nm', 'GCACGAGATTCACAT'),
        ('R32', 'fluor_488nm', 'TCCGCCTGGATCCAA'),
        ('R33', 'fluor_560nm', 'AAGCGGACGTCAGTT'),
        ('R34', 'fluor_640nm', 'ACGGAAGCTCTAATC'),
        ('R35', 'fluor_488nm', 'TAGAATTCGCGGCCA'),
        ('R36', 'fluor_560nm', 'TGATACGCGAACTGT'),
        ('R37', 'fluor_640nm', 'AATCAGGATACGGCG'),
        ('R38', 'fluor_488nm', 'GGTCTCGGATATGCA'),
        ('R39', 'fluor_560nm', 'TGCGTACGACAATGT'),
        ('R40', 'fluor_640nm', 'ACAACGTGCGCCTGA'),
        ('R41', 'fluor_488nm', 'TAAGATGTCGCCGTG'),
        ('R42', 'fluor_560nm', 'TATCATACGACGTGG'),
        ('R43', 'fluor_640nm', 'TGGATCATTCTACGG'),
        ('R44', 'fluor_488nm', 'TCTAGTCTCACGCAA'),
        ('R45', 'fluor_560nm', 'AGATTGGCGCGACAT'),
        ('R46', 'fluor_640nm', 'AACGCGTAGGCTAGA'),
        ('R47', 'fluor_488nm', 'TTGACGAAGCCTGAC'),
        ('R48', 'fluor_560nm', 'ACTACGCGCGGATAT'),
        ('R49', 'fluor_640nm', 'TTCCACTGCGTAGTG'),
        ('R50', 'fluor_488nm', 'ATTGTGCGAGGTCCA'),
        ('R51', 'fluor_560nm', 'TCAATTGGGCGCTCA'),
        ('R52', 'fluor_640nm', 'CTCTAGGATCGTACG'),
        ('R53', 'fluor_488nm', 'TGACTCGCACTGGGA'),
        ('R54', 'fluor_560nm', 'TTATTACGGAGGCTG'),
        ('R55', 'fluor_640nm', 'TAGCTCGACAGGGTT'),
        ('R56', 'fluor_488nm', 'TGCGGACCTCTTATA'),
        ('R57', 'fluor_560nm', 'GCTAGGCACCTGGAT'),
        ('R58', 'fluor_640nm', 'GCGGCATTGACTAAG'),
        ('R59', 'fluor_488nm', 'TCGGGCGGCTATAGA'),
        ('R60', 'fluor_560nm', 'TGTTCGCCAGGATTG'),
        ('R61', 'fluor_640nm', 'TGGTCGAGGTATTCG'),
        ('R62', 'fluor_488nm', 'ACATTCGTTATCGCG'),
        ('R63', 'fluor_560nm', 'ATATTCAGGCTAGCG'),
        ('R64', 'fluor_640nm', 'ACGCTAGTAACGTCT'),
        ('R65', 'fluor_488nm', 'CCTAGGCATATTCGA'),
        ('R66', 'fluor_560nm', 'AGTACTCTTACGCGA'),
        ('R67', 'fluor_640nm', 'TAACGAGTGCTAACG'),
        ('R68', 'fluor_488nm', 'TATTCGTGAATCGCG'),
        ('R69', 'fluor_560nm', 'ACGCGCGTCGTTATA'),
        ('R70', 'fluor_640nm', 'TATTGGCTACGGAGC'),
        ('R71', 'fluor_488nm', 'ATACGGACAATTCGG'),
        ('R72', 'fluor_560nm', 'GTCGTAATAGCACGC'),
        ('R73', 'fluor_640nm', 'TCTCGTGATACGGAA'),
        ('R74', 'fluor_488nm', 'CTGCGAGACATTGCT'),
        ('R75', 'fluor_560nm', 'CGCTAATTGCCGAGT'),
        ('R76', 'fluor_640nm', 'TGCCTTGTCACGCAA'),
        ('R77', 'fluor_488nm', 'AGCACCTTTACATCG'),
        ('R78', 'fluor_560nm', 'CCGGTTAGCAATGCT'),
        ('R79', 'fluor_640nm', 'CCACGATCGCTGATT'),
        ('R80', 'fluor_488nm', 'TGGCCCATACTCGTT'),
        ('R81', 'fluor_560nm', 'CTTCAGCGATAGGCA'),
        ('R82', 'fluor_640nm', 'GTAAGCTATAAGCGG'),
        ('R83', 'fluor_488nm', 'ATTGCGAGGCCACTA'),
        ('R84', 'fluor_560nm', 'ACTTCAACGAGTGAC'),
        ('R85', 'fluor_640nm', 'AGTCGTATACGCGTG'),
        ('R86', 'fluor_488nm', 'CACTTAACTAGACGC'),
        ('R87', 'fluor_560nm', 'CACAGACGCATATTC'),
        ('R88', 'fluor_640nm', 'CGCGTCACGTTCAAT'),
        ('R89', 'fluor_488nm', 'CTGATAAGCTTGGAC'),
        ('R90', 'fluor_560nm', 'ATGCACGACTTCGGA'),
        ('R91', 'fluor_640nm', 'GCTGTCGCCTATAAG'),
        ('R92', 'fluor_488nm', 'TTAATCGGCCGAACG'),
        ('R93', 'fluor_560nm', 'AATTCGACGGTCTAG'),
        ('R94', 'fluor_640nm', 'GCGACTGGAACACTT'),
        ('R95', 'fluor_488nm', 'TGATACGCTTGGCAT'),
        ('R96', 'fluor_560nm', 'TGGCAGATCGCTAAT'),
        ('R97', 'fluor_640nm', 'CGCGCAGACGTTAAT'),
        ('R98', 'fluor_488nm', 'CATGATACAATCCGG'),
        ('R99', 'fluor_560nm', 'GGCGTATCGAACATT'),
        ('R100', 'fluor_640nm', 'ATTGCACAACCGATG'),
        ('R101', 'fluor_488nm', 'CGTACTGCAGATCGA'),
        ('R102', 'fluor_560nm', 'TACGACTCGGTCATT'),
        ('R103', 'fluor_640nm', 'TCGTAGAGCCTAATC'),
        ('R104', 'fluor_488nm', 'GGCGCACAACGAAAT'),
        ('R105', 'fluor_560nm', 'TCGATGAATCGTCGA'),
        ('R106', 'fluor_640nm', 'TGTGCGAACACTCGT'),
        ('R107', 'fluor_488nm', 'AACGATCGCGGTGAT'),
        ('R108', 'fluor_560nm', 'ATAGAATATCGGCCG'),
        ('R109', 'fluor_640nm', 'TATGTGACTACGCAC'),
        ('R110', 'fluor_488nm', 'CGTCGTATCCAGTAT'),
        ('R111', 'fluor_560nm', 'CTATTATCGCCGAGA'),
        ('R112', 'fluor_640nm', 'CGATTAGTCGTCACT'),
        ('R113', 'fluor_488nm', 'ACTCCGAATGCTACG'),
        ('R114', 'fluor_560nm', 'GGTTACACGCGACTA'),
        ('R115', 'fluor_640nm', 'ATGTAACCAAGCGTC'),
        ('R116', 'fluor_488nm', 'TCAGTTACCGGTGTA'),
        ('R117', 'fluor_560nm', 'TCCAGCTTACGTTCG'),
        ('R118', 'fluor_640nm', 'AGTACAACGGCATTC'),
        ('R119', 'fluor_488nm', 'AACCTGCCAACGCTT'),
        ('R120', 'fluor_560nm', 'ATCTAACCGGGACGT'),
        ('R121', 'fluor_640nm', 'CTAGTGATTCGCACC'),
        ('R122', 'fluor_488nm', 'CTGGTCATAACTCCC'),
        ('R123', 'fluor_560nm', 'AATCAGGCGTGTTAC'),
        ('R124', 'fluor_640nm', 'CGCTTGCGATTTCGA'),
        ('R125', 'fluor_488nm', 'CTTTCGTACCTACGG'),
        ('R126', 'fluor_560nm', 'ATTATCAACGAGACC'),
        ('R127', 'fluor_640nm', 'CATGTACGACCGTAA'),
        ('R128', 'fluor_488nm', 'AATCCTCGGTCGTAG'),
        ('R129', 'fluor_560nm', 'GGGTATGATCGCAGT'),
        ('R130', 'fluor_640nm', 'TGTGTTTCGCAACGA'),
        ('R131', 'fluor_488nm', 'TCAGAGCTACGACTA'),
        ('R132', 'fluor_560nm', 'AGTAGTCGCTAGGCT'),
        ('R133', 'fluor_640nm', 'GACGAGTAACAACTG'),
        ('R134', 'fluor_488nm', 'TTAAGGATCAGACCG'),
        ('R135', 'fluor_560nm', 'TGACGAAGGGCTCTA'),
        ('R136', 'fluor_640nm', 'CGTGATGAACCGCTT'),
        ('R137', 'fluor_488nm', 'CTGTTAACGAGGCCT'),
        ('R138', 'fluor_560nm', 'CGTACCGACATCGAA'),
        ('R139', 'fluor_640nm', 'GAACGCCGCTATTAA'),
        ('R140', 'fluor_488nm', 'GTGACTGGAGCGTTA'),
        ('R141', 'fluor_560nm', 'ATAGCGTTCCACGGG'),
        ('R142', 'fluor_640nm', 'AATACGCGAGCTCAG'),
        ('R143', 'fluor_488nm', 'TTGGCGCCGCTAAAT'),
        ('R144', 'fluor_560nm', 'CGCTTTTACGTACCA'),
        ('R145', 'fluor_640nm', 'GTACAACCTCGGTTA'),
        ('R146', 'fluor_488nm', 'AGACTTCGCCGTACA'),
        ('R147', 'fluor_560nm', 'TGATTTGGAGGACGA'),
        ('R148', 'fluor_640nm', 'CGTCGATAGGACTCA'),
        ('R149', 'fluor_488nm', 'ATCTTGCGCGTTGGG'),
        ('R150', 'fluor_560nm', 'TTCAGAGCCGGTAGA'),
        ('R151', 'fluor_640nm', 'CACCGGTCTACGTAA'),
        ('R152', 'fluor_488nm', 'AAAGTTCGCGACCTA'),
        ('R153', 'fluor_560nm', 'TTGCTCGTCCTACGA'),
        ('R154', 'fluor_640nm', 'AAGAACGGGTCACTC'),
        ('R155', 'fluor_488nm', 'ATGCGTCGGTGCTAT'),
        ('R156', 'fluor_560nm', 'GCTTGGGTCGAGCAA'),
        ('R157', 'fluor_640nm', 'GTTACGGAACAGCGA'),
        ('R158', 'fluor_488nm', 'GCGTGATCGGTGATT'),
        ('R159', 'fluor_560nm', 'CGCCGAGCGATATAA'),
        ('R160', 'fluor_640nm', 'TGCGTTGACGCAATA'),
        ('R161', 'fluor_488nm', 'GTGTGCAGAATAGCG'),
        ('R162', 'fluor_560nm', 'ACATTCGCAAGGAAC'),
        ('R163', 'fluor_640nm', 'AGGATTCCAAGGCGA'),
        ('R164', 'fluor_488nm', 'AACTGAAGCAACGGT'),
        ('R165', 'fluor_560nm', 'AGTATCAGCTTGCAC'),
        ('R166', 'fluor_640nm', 'CAGCCATAGATCTAG'),
        ('R167', 'fluor_488nm', 'CCTCAAGTTATGCGC'),
        ('R168', 'fluor_560nm', 'TCCATATACGAGGTG'),
        ('R169', 'fluor_640nm', 'GGAGTAATACGACGC'),
        ('R170', 'fluor_488nm', 'CCGTAAAGTCGACTG'),
        ('R171', 'fluor_560nm', 'AACTATACCGTCCAA'),
        ('R172', 'fluor_640nm', 'TATTGGCTGTTACGA'),
        ('R173', 'fluor_488nm', 'TCTTTACCGACTCAT'),
        ('R174', 'fluor_560nm', 'GCTGTTCGCGTGATA'),
        ('R175', 'fluor_640nm', 'CCCGTAGATGACTTT'),
        ('R176', 'fluor_488nm', 'ATTATCGGTTCCAGC'),
        ('R177', 'fluor_560nm', 'TCGGCTTTTCGCATG'),
        ('R178', 'fluor_640nm', 'GAAGTCAAATCCGTC'),
        ('R179', 'fluor_488nm', 'CCCAATCGCTTTGGC'),
        ('R180', 'fluor_560nm', 'ACATCCCGCGTATAT'),
        ('R181', 'fluor_640nm', 'TTCAATCCGTCGACA'),
        ('R182', 'fluor_488nm', 'GTACCCGTGTGGTAA'),
        ('R183', 'fluor_560nm', 'CAAGGTGATTCGTGT'),
        ('R184', 'fluor_640nm', 'TGAGTAAAGACCGCC'),
        ('R185', 'fluor_488nm', 'AAACCGCTATATGCC'),
        ('R186', 'fluor_560nm', 'TCTGCGTCGAGAGTT'),
        ('R187', 'fluor_640nm', 'TACGAGTGTACCGAA'),
        ('R188', 'fluor_488nm', 'CCCTTTACGACTACC'),
        ('R189', 'fluor_560nm', 'GAGCGTTTGATTGAC'),
        ('R190', 'fluor_640nm', 'GTCCCGCCTTATAAA'),
        ('R191', 'fluor_488nm', 'CTAGCACTGGTCGTC'),
        ('R192', 'fluor_560nm', 'ACATAACATCCCGGT'),
        ('R193', 'fluor_640nm', 'TTGGGGACGTATCGA'),
        ('R194', 'fluor_488nm', 'TAGACCGCCTATCAT'),
        ('R195', 'fluor_560nm', 'TCTTCAACTTGCTCG'),
        ('R196', 'fluor_640nm', 'TTGAACATTGCCGCC'),
        ('R197', 'fluor_488nm', 'GCGAATACCATACGC'),
        ('R198', 'fluor_560nm', 'CTATAGACGAATAGC'),
        ('R199', 'fluor_640nm', 'CAAATCCGCGTGTTA'),
        ('R200', 'fluor_488nm', 'GCTCTAAGTGTCCGA'),
        ('R201', 'fluor_560nm', 'CGCCGTGACGTATTT'),
        ('R202', 'fluor_640nm', 'AAGGTGTGTGAACGC'),
        ('R203', 'fluor_488nm', 'TCTAGCGGACATGAC'),
        ('R204', 'fluor_560nm', 'CCAGTTGCCGTTACA'),
        ('R205', 'fluor_640nm', 'GTCGATTGATAAGGC'),
        ('R206', 'fluor_488nm', 'TTATCATGCCGCCAC'),
        ('R207', 'fluor_560nm', 'AAGACGTGCTCCAGA'),
        ('R208', 'fluor_640nm', 'AAGTGAGGATTACGC'),
        ('R209', 'fluor_488nm', 'TATTAACGTGCCTGC'),
        ('R210', 'fluor_560nm', 'AGCGTCATCGAAATC'),
        ('R211', 'fluor_640nm', 'AACTTCCCGACTCCG'),
        ('R212', 'fluor_488nm', 'TTATGCCGGAGCAAA'),
        ('R213', 'fluor_560nm', 'GATACTAGCAGCCGC'),
        ('R214', 'fluor_640nm', 'ATCGCAATCGTACCG'),
        ('R215', 'fluor_488nm', 'CTACAATCGAGCGAA'),
        ('R216', 'fluor_560nm', 'TTGGCAGTCGCACAG'),
        ('R217', 'fluor_640nm', 'CAATGGATTACCTCG'),
        ('R218', 'fluor_488nm', 'TAGTCGCTCGAATAC'),
        ('R219', 'fluor_560nm', 'GTCCGTATTGCGCTA'),
    ]

    df = pd.DataFrame(data, columns=['ro_id', 'fluor', 'readout_seq'])
    # Filter out rows with 'R3' or 'R5'
    df_filtered = df[~df['ro_id'].isin(probes_to_remove)]
    # df = df_filtered.reset_index(drop=True)

    # return df
    return df_filtered
