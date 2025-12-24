# Antibody Sequence Characterisation Tool

Hello! 

I've started putting together these scripts for my personal projects, 
mainly to keep track of what I've learned and make sure I don't forget it all!

This is a simple Python tool to help analyse antibody sequences from a FASTA file.
It calculates various physicochemical properties and stores them into a CSV file.

## Features
* **Molecular Weight** calculation
* **Isoelectric Point** determination
* **Net Charge** at a specific pH
* **GRAVY Score** and **Aromaticity**

## How to run it
First, make sure you've got the necessary bits installed:
```bash
pip install biopython pandas numpy
```

Then, run the script from your terminal like so:
```bash
python your_script_name.py --fastaFile input.fasta --csvFile results.csv --pH_value 7.4
```
(Note: --pH_value is optional and defaults to 7.4 if not specified)


## Example FASTA file 
These are sample data, which includes amino acid sequence of therapeutics used for breast cancer. 

only Fv regions: [Ab_breastcancer_seq_Fv.fasta](./Ab_breastcancer_seq_Fv.fasta)<br>
Full sequences : [Ab_breastcancer_seq_full.fasta](./Ab_breastcancer_seq_full.fasta)

The sequences were retrieved from Thera-SAbDab and KEGG DRUG

