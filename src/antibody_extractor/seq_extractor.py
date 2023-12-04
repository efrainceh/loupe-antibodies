import os
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class SeqExtractor:

    """
    
        Class for extracting the protein sequences from fasta DNA files. It extracts clonotype ids from the input csv files, then gets the DNA sequences associated
        with those ids from the fasta file. Finally, it translates the DNA to their protein sequence.
        
        Inputs: 
        
            - folder where all the clonotype csv files are located. The program will read all csv files in the folder. The example below is from loupe V(D)J
            but the only requirement is for a "clonotype_id" column
            - fasta file with the DNA sequences (placed in the same folder as clonotype csv)

        Outputs:

            - A list protein SeqRecords
            - A list of all clonotypes in the csv inputs

        loupe csv file example:

            clonotype_id	consensus_ids	frequency	proportion
            clonotype139	clonotype139_consensus_293;clonotype139_consensus_294	1	0.00018511662347278800
            clonotype211	clonotype211_consensus_432	1	0.00018511662347278800
            clonotype214	clonotype214_consensus_435	1	0.00018511662347278800
            clonotype243	clonotype243_consensus_464	1	0.00018511662347278800
            clonotype312	clonotype312_consensus_533	1	0.00018511662347278800
            clonotype939	clonotype939_consensus_1708;clonotype939_consensus_1709	1	0.00018511662347278800

        fasta file example (sequences were shortened for readability): 

            >clonotype1_consensus_1
            TGGGGATCTCCTCACTAGAGCCCCCATCAGAGCATGGCTGTCCTGGTGCTGTTCCTCTGCCTGGTTGCATTTCCAAGCTGTGTCCTGTCCCAGGTGCAGCTGAAGGAGTCAGGACCTGGCCTGGTGGCGCCCTCACAGAG
            >clonotype1_consensus_2
            TTTTTGGGGGACCAATATTGAAAATAATAGACTTGGTTTGTGAATTATGGCCTGGACTTCACTTATACTCTCTCTCCTGGCTCTCTGCTCAGGAGCCAGTTCCCAGGCTGTTGTGACTCAGGAATCTGCACTCACCACAT
            >clonotype2_consensus_1
            AGTAACCACTGTGTCTGATATGGGGAACCGACGATCAGTGTCCTCTCCAAAGTCCCTGAACACACTGACTCTAACCATGGAATGGAGTTGGATATTTCTCTTTCTCCTGTCAGGAACTGCAGGTGTCCACTCTGAGGTCCA
            >clonotype2_consensus_2
            GAAATGCATCACACCAGCATGGGCATCAAAATGGAGTCACAGATTCAGGTCTTTGTATTCGTGTTTCTCTGGTTGTCTGGTGTTGACGGAGACATTGTGATGACCCAGTCTCACAAATTCATGTCCACATCAGTAGGAGAC
            >clonotype3_consensus_1
            TTGTTGGGGGACCAATATTGAAAATAATAGACTTGGTTTGTGAATTATGGCCTGGACTTCACTTATACTCTCTCTCCTGGCTCTCTGCTCAGGAGCCAGTTCCCAGGCTGTTGTGACTCAGGAATCTGCACTCACCACATC
    
    """

    def __init__(self, folder, fasta):

        self.folder = folder
        self.fasta = fasta
        self.clonotypes = []

    def get_sequences(self):

        # Get a list of clonotypes from all the csv files in folder
        self.clonotypes = self._get_clonotypes(self.folder)

        # Get the DNA sequences associated with each clonotype.
        dna_sequences = self._extract_dna_sequences_from_fasta(self.clonotypes, os.path.join(self.folder, self.fasta))

        # Translate the DNA sequences to protein.
        protein_sequences = self._translate_sequences(dna_sequences)

        return protein_sequences
    
    def get_clonotypes(self):

        if self.clonotypes:

            return self.clonotypes
        
        else:

            return self._get_clonotypes(self.folder)

    def _get_clonotypes(self, folder):

        # Get all csv files in folder
        files = [os.path.join(folder, file) for file in os.listdir(folder) if file.endswith(".csv")]

        return [clonotype for file in files for clonotype in self._clonotypes_from_csv(file)]
    
    def _clonotypes_from_csv(self, file):

        # Return a list of all clonotypes in the csv file
        CLONOTYPE_COL = "clonotype_id"
        df = pd.read_csv(file)

        return df[CLONOTYPE_COL].tolist()
    
    def _extract_dna_sequences_from_fasta(self, clonotypes, fasta):

        # Return a list of sequece records that are found in clonotypes
        records = SeqIO.parse(fasta, "fasta")

        return [record for record in records if self._in_list(record.id, clonotypes)]

    def _in_list(self, id, clonotypes):

        # Clonotype id is found before the first '_' in record.id
        short_id = id.split('_', 1)[0]

        # Needs to be an exact match, otherwise clonotype_1 is found in clonotype_11, clonotype_12...
        return any([(clonotype == short_id) for clonotype in clonotypes])

    def _translate_sequences(self, sequences):

        return [self._translate(sequence) for sequence in sequences]

    def _translate(self, dna_record):

        dna = dna_record.seq

        # Find start and end position of coding region, then translate
        start, end = self._get_cds(dna, 0)
        protein = dna[start:end].translate()

        # If protein has a stop codon, then loop until it finds a sequence without one
        i = 0
        while protein.find('*') != -1 and i < 3:

            start, end = self._get_cds(dna, start + 1)
            protein = dna[start:end].translate()
            i += 1

        return SeqRecord(
                Seq(protein),
                id=dna_record.id,
                name=dna_record.name,
                description=dna_record.description,
                )

    def _get_cds(self, dna, start):

         # Because the DNA sequence doesn't generally start with ATG, just calling translate(from BioPython) returns the wrong protein sequence. 
        start = dna.find("ATG", start)
                
        # The DNA sequence needs to be a multiple of 3 to avoid a warning. Find last position from the right that it's a multiple of 3
        end = len(dna) - len(dna[start:]) % 3

        return start, end
    