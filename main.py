import argparse
import os

from src.antibody_extractor import AntibodyPairer, CDRFinder, SeqExtractor
from src.biophi_wrapper import BioPhiWrapper

def main(args):

    # Creating folder and file structure
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    runs_folder = os.path.join(ROOT_DIR, "runs")
    if not os.path.isdir(runs_folder):
        os.mkdir(runs_folder)
    input_folder = os.path.join(runs_folder, args.folder)
    results_folder = os.path.join(input_folder, "results")
    if not os.path.isdir(results_folder):
        os.mkdir(results_folder)
    output_fasta = os.path.join(results_folder, '{}.fasta'.format(args.folder))
    output_csv = os.path.join(results_folder, '{}.csv'.format(args.folder))

    # Extract protein sequences as a SeqRecord list from clonotype csv and fasta DNA files
    print("Extracting clonotype sequences from fasta file...")
    sq = SeqExtractor(input_folder, args.fasta)
    proteins = sq.get_sequences()

    # Get a list of antibody pairs
    print("Pairing antibodies...")
    ap = AntibodyPairer(proteins)
    pairs = ap.get_chain_pairs()
    print("Writing output fasta file...")
    ap.to_fasta(output_fasta, True)

    # Find chains CDRs
    print("Finding CDRs...")
    cf = CDRFinder(pairs)
    cdr_df = cf.find_cdrs()


    if args.biophi:

        if args.oasis_db == "":

            print("ERROR: A path for the oasis database is required")
            exit()

        # Calculate humannes with BioPhi
        print("Calculating humanness with BioPhi...")
        bw = BioPhiWrapper(output_fasta)
        humanness_df = bw.humannes_report(args.oasis_db)
        cdr_df = cdr_df.join(humanness_df)

    # Save to output file
    print("Creating csv output file...")
    cdr_df.to_csv(output_csv)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Script creates a csv file of antibody pairs (variable region and CDRs) from a list of clonotypes and a fasta DNA file. It can optionally also run BioPhi"
    )
    parser.add_argument("--folder", required=True, type=str,
                        help='folder where the inputs are located: csv clonotypes and fasta file')
    parser.add_argument("--fasta", required=True, type=str,
                        help='name of fasta file with clonotype DNA sequences, should be located in folder')
    parser.add_argument("--biophi", required=False, default=False, type=bool,
                        help='whether BioPhi should be run on the extracted antibody sequences')
    parser.add_argument("--oasis_db", required=False, default="", type=str,
                        help='path to the oasis database')
    args = parser.parse_args()
    
    main(args)
