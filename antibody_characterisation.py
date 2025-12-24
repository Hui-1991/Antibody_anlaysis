######################################### Import Package/Module###################################################
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np 
import argparse
import os
import sys

######################################### Main Code ##############################################################

def main(args):
    # Initialize an empty list to store the data from each antibody sequence       
    results = []

    # Iterate over each sequence record in the input FASTA file
    for record in SeqIO.parse(args.fastaFile,'fasta'): 
        # Extract antibody name
        ab_name = record.id.split("_")[0]
        
        # Retrieve the FASTA description line
        desc = record.description
        # In case description does not exist in 'FASTA file 
        target_info = np.nan
        chain = np.nan
        # Extract each information from the description 
        try: 
            target_info = desc.split("Target:")[1].split()[0]
            chain_info = desc.split("Chain:")[1].split()[0]
            isotype = desc.split("Isotype:")[1].split()[0]
        except IndexError:
            # If any variable is missing, keep the default NaN vaules
            pass   
        
        # Ensure the amino acid sequence is uppercase 
        seq = record.seq.upper()
        # Perform protein pysicochemical analysis 
        analysed_seq = ProteinAnalysis(seq)            

        # Store all properties in a dictionary 
        data = {
                "Therapeutics": ab_name,
                "Target": target_info,
                "Chain": chain,
                "Isotype": isotype,
                "Length": len(seq),
                "Molecular Weight": round(analysed_seq.molecular_weight(), 2),
                "Isoelectric_Point": round(analysed_seq.isoelectric_point(), 2),
                "Instability Index": round(analysed_seq.instability_index(), 2),
                "GRAVY Score": round(analysed_seq.gravy(), 2),
                "Aromaticity": round(analysed_seq.aromaticity(), 2),
                f"Net Charge at {args.pH_value}": round(analysed_seq.charge_at_pH(args.pH_value), 2)
            }
    
        # Append the data into result list 
        results.append(data)
        
    
    
    # Convert the list of dictionaries into DataFrame
    df = pd.DataFrame(results)
    # Write the DataFrame into a CSV file 
    df.to_csv(args.csvFile, index=False)

    


######################################### Script execution  ######################################################

# Run the main function when the script is executed directly. 
if __name__ == "__main__":    

    # Create an argument parser for command-line execution 
    parser = argparse.ArgumentParser(
        description='Characterise antibody sequences',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    ) 

    # Define the input FASTA file argument
    parser.add_argument('--fastaFile', required=True, help='Path to the input FASTA file')         
    # Define the output CSV file argument
    parser.add_argument('--csvFile', required=True, type=str, help='Name for output CSV file')
    # Define an optional pH value for net charge calculation (default = 7.4)
    parser.add_argument('--pH_value', default=7.4, 
                        type=float, 
                        help='pH value for net charge calculation')

    # Parse the command-line arugemnts 
    args = parser.parse_args()

    # Check if the input FASTA file exsit and execute the main function 
    if not os.path.isfile(args.fastaFile):
        sys.exit(f"ERROR: FASTA file '{args.fastaFile}' not found.")
    else: 
        main(args)