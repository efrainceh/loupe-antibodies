import os
import subprocess
import pandas as pd 

class BioPhiWrapper:

    """
    
        A wrapper for BioPhi (https://github.com/Merck/BioPhi/blob/main/README.md).

        Inputs:

           - A fasta file of antibody pairs. Each antibody ID needs to be a pair where the heavy chain is labelled with the suffix "_HV"
            and the light chain with "LV"

        Outputs:

            - The BioPhi xlsx output file
            - A dataframe with COLUMNS, with the values converted to percentages for readability

                                        OASis Identity  Heavy OASis Identity  Light OASis Identity
                    Antibody                                                                 
                    clonotype1220            61.0                  63.0                  58.0
                    clonotype1241            43.0                  38.0                  49.0
                    clonotype1257            51.0                  51.0                  52.0
                    clonotype1263            48.0                  47.0                  49.0
                    clonotype1360            49.0                  49.0                  48.0

    """

    COLUMNS = ["Antibody", "OASis Identity", "Heavy OASis Identity", "Light OASis Identity"]
    flags = {
        "program" : "biophi", 
        "oasis" : "oasis",
        "fasta" : "",
        "db_flag" : "--oasis-db", 
        "db_path" : "", 
        "output_flag" : "--output",
        "output_file" : ""
        }

    def __init__(self, fasta_file):
        
        self.fasta_file = fasta_file

    def humannes_report(self, database):

        # If the run is succesful, it returns the name of the xlsx output file
        biophi_file = self._run_biophi(self.fasta_file, database)

        return self._extract_humanness_percentages(biophi_file)

    def _run_biophi(self, file, database):

        # Update flag dictionary
        self.flags["fasta"] = file
        self.flags["db_path"] = database
        output_file = '{}.xlsx'.format(os.path.splitext(file)[0])
        self.flags["output_file"] = output_file
        flags_lst = [value for key, value in self.flags.items()]

        # Run BioPhi
        if subprocess.run(flags_lst).returncode == 0:

            return output_file
        
        else:

            print("ERROR: Something went wrong with BioPhi")
            exit()
    
    def _extract_humanness_percentages(self, file):

        # Read the output file created by BioPhi
        xlsx = pd.ExcelFile(file)
        df = xlsx.parse("Overview")

        # Select columns of interests, multiply by 100 for csv readability
        df = df[self.COLUMNS]
        df.set_index("Antibody", inplace=True)
        df = df.mul(100).round(0)

        return df