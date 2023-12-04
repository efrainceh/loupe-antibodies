from abnumber import Chain
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class AntibodyPairer:

    """
    
        Creates a list of antibody pairs, where each pair is represented by a dictionary. The protein sequences are represented by Seq
        objects so that their CDRs can be easily found.

        Inputs:

           - A list of protein sequence SeqRecords, such as 
                Seq(protein sequence),
                    id=protein id,
                    name=protein name,
                    description=protein description,
                    )

        Outputs:

             - A list of antibody pairs represented as dictionary of type 
                { 
                    id : clonotype,
                    chains: 1 o 2,
                    heavy : protein Seq object of heavy chain
                    light : protein Seq object of light chain
                }
            - A fasta file with either all the clonotype protein sequences or just those with heavy-light chain pairings
            
    """

    ID = "id"
    CHAINS = "chains"
    HEAVY = "heavy_chain"
    LIGHT = "light_chain"

    def __init__(self, proteins):
        
        self.proteins = proteins
        self.chain_pairs = []

    def get_chain_pairs(self):

        # Create a dictionary of type { clonotype : [protein Seq objects]}
        protein_dict = self._separate_by_clonotypes(self.proteins)

        for clonotype, sequences in protein_dict.items():

            # Get all heavy-light chain pairings possible for that clonotype
            per_clonotype_pairs = self._get_pairs(sequences, clonotype)
            self.chain_pairs.extend(per_clonotype_pairs)

        # Add "_number" to clonotypes ids that are repeated
        self._rename_chain_pairs_ids()

        return self.chain_pairs
    
    def to_fasta(self, output_file, only_pairs=False):

        # Outputs the antibody pairs in a fasta format. If only_pairs = True, then only clonotypes with heavy-light chain pairings
        # are saved. Otherwise it prints all clonotypes.

        records = []

        for pair in self.chain_pairs:

            if only_pairs and pair[self.CHAINS] == 2:
                    
                heavy_chain = self._create_record(pair[self.HEAVY], '{}_VH'.format(pair[self.ID]))
                light_chain = self._create_record(pair[self.LIGHT], '{}_VL'.format(pair[self.ID]))
                records.extend([heavy_chain, light_chain])
                
            elif not only_pairs:

                if pair[self.HEAVY]:

                    heavy_chain = self._create_record(pair[self.HEAVY], '{}_VH'.format(pair[self.ID]))
                    records.append(heavy_chain)

                if pair[self.LIGHT]:

                    light_chain = self._create_record(pair[self.LIGHT], '{}_VL'.format(pair[self.ID]))
                    records.append(light_chain)

        SeqIO.write(records, output_file, "fasta")

    def _separate_by_clonotypes(self, proteins):

        # Return a dictionary of clonotypes associated with their protein sequences as Seq object
        dict = {}

        for protein in proteins:

            # Protein clonotype id  example = clonotype2718_consensus_1. We are only interested in the clonotype  
            id = protein.id.split('_', 1)[0]

            if id in dict.keys():

                dict[id].append(protein.seq)
                
            else:

                dict[id] = [protein.seq]

        return dict

    def _get_pairs(self, sequences, clonotype):

        # Get a dictionary with the heavy and light chains in the clonotype
        heavy_and_light_chains = self._find_heavy_and_light_chains(sequences)

        # Create all the heavy-light pairings possible
        heavy_chains = heavy_and_light_chains[self.HEAVY]
        light_chains = heavy_and_light_chains[self.LIGHT]

        if len(heavy_chains) > 0 and len(light_chains) > 0:
            
            pairs = [self._create_pair_dict(clonotype, 2, heavy_chain, light_chain) for light_chain in light_chains for heavy_chain in heavy_chains]

        elif heavy_chains:

            pairs = [self._create_pair_dict(clonotype, 1, heavy_chain, "") for heavy_chain in heavy_chains]

        elif light_chains:

            pairs = [self._create_pair_dict(clonotype, 1, "", light_chain) for light_chain in light_chains]

        return pairs

    def _find_heavy_and_light_chains(self, sequences):

        # In some cases sequences have more than two protein sequences in it, representing more than one heavy-light chain pairing.
        # That's why we use lists and append
        temp_dict = { self.HEAVY : [], self.LIGHT : [] }

        for seq in sequences:

            chain  = Chain(seq, scheme='imgt')
            
            if chain.is_heavy_chain():

                temp_dict[self.HEAVY].append(chain)

            if chain.is_light_chain():

                temp_dict[self.LIGHT].append(chain)

        return temp_dict
    
    def _create_pair_dict(self, clonotype, chains, heavy_chain, light_chain):

        return {
            self.ID : clonotype,
            self.CHAINS : chains,
            self.HEAVY : heavy_chain,
            self.LIGHT : light_chain
        }
    
    def _create_record(self, sequence, id):

        return SeqRecord(
                Seq(str(sequence)),
                id=id,
                name=id,
                description=id,
            )
    
    def _rename_chain_pairs_ids(self):

        ix = 0
        while ix < len(self.chain_pairs) - 1:

            current_id = self.chain_pairs[ix][self.ID]
            next_id = self.chain_pairs[ix + 1][self.ID]
            ids_to_change = []

            while current_id == next_id:

                ix += 1
                ids_to_change.append(ix)
                
                if ix == len(self.chain_pairs) - 1:
                    break
                
                current_id = self.chain_pairs[ix][self.ID]
                next_id = self.chain_pairs[ix + 1][self.ID]

            # Rename IDs, except the first antibody pair with the specific ID
            self._rename_ids(ids_to_change)

            ix += 1

    def _rename_ids(self, indexes):

        # Add "_number" to clonotype ids, starting at number = 2
        count = 2
        for ix in indexes:

            id = self.chain_pairs[ix][self.ID]
            self.chain_pairs[ix][self.ID] = '{}_{}'.format(id, count)
            count += 1