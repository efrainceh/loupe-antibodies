import pandas as pd

class CDRFinder:

    """

        Finds the variable regions and CDRs of each antibody pair and returns them in a dataframe 

        Inputs:

            - A list of antibody pairs represented as dictionary of type 
                { 
                    id : clonotype,
                    chains: 1 o 2,
                    heavy : protein Seq object of heavy chain
                    light : protein Seq object of light chain
                }

        Outputs:

            - A dataframe with COLUMNS: 

                clonotype       HV          H_CDR1      H_CDR2      H_CDR3          LV          L_CDR1      L_CDR2      L_CDR3                                                                                                                                                                   
                clonotype1      EVQLV...    GVTFCRNA    IGSGARYT    AIHKGKDWFDY     DIQMTQS...  QYIGTW      DTT         QPLYSSPPT
                clonotype2                                                          DIQMNQS...  QNIYVW      KWS         QLGQLYPLT
                clonotype3                                                          DIQMTQS...  WNIYSN      WAT         QHFLLTPFT

    """

    COLUMNS = ["clonotype", "HV", "H_CDR1", "H_CDR2", "H_CDR3", "LV", "L_CDR1", "L_CDR2", "L_CDR3"]
    ID = "id"
    CHAINS = "chains"
    HEAVY = "heavy_chain"
    LIGHT = "light_chain"

    def __init__(self, antibodies):
        
        self.antibodies = antibodies

        # Create an empty dataframe with COLUMNS
        self.df = pd.DataFrame(columns = self.COLUMNS)
        

    def find_cdrs(self):

        for antibody in self.antibodies:

            self._add_to_df(antibody)

        self.df.set_index("clonotype", inplace=True)

        return self.df

    def _add_to_df(self, antibody):

        # Find CDRs and return the data as a dictionary
        row_dict = self._create_cdr_dict(antibody)

        # Add the dictionary to the dataframe
        new_row = pd.DataFrame(row_dict)
        self.df = pd.concat([self.df, new_row], ignore_index = True)

    def _create_cdr_dict(self, antibody):

        h_seq, h_cdr1, h_cdr2, h_cdr3 = self._get_seq_data(antibody[self.HEAVY])
        l_seq, l_cdr1, l_cdr2, l_cdr3 = self._get_seq_data(antibody[self.LIGHT])

        return { 
            self.COLUMNS[0] : [antibody[self.ID]], 
            self.COLUMNS[1] : [h_seq], self.COLUMNS[2] : [h_cdr1], self.COLUMNS[3] : [h_cdr2], self.COLUMNS[4] : [h_cdr3],
            self.COLUMNS[5] : [l_seq], self.COLUMNS[6] : [l_cdr1], self.COLUMNS[7] : [l_cdr2], self.COLUMNS[8] : [l_cdr3]
            }
    
    def _get_seq_data(self, seq):

        # Start with empty values
        chain_seq = ""
        cdr1 = ""
        cdr2 = ""
        cdr3 = ""

        # Ensure we are not trying to get data from an empty seq or a string.
        if not isinstance(seq, str) and bool(seq):

            chain_seq = seq.seq
            cdr1 = seq.cdr1_seq
            cdr2 = seq.cdr2_seq
            cdr3 = seq.cdr3_seq

        return chain_seq, cdr1, cdr2, cdr3