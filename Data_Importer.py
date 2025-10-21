# Data_importer.py
import pandas as pd
import argparse

def load_samples(samples_file):
    """Reads the CSV and returns grouped VCF file paths as a dict"""
    samples = pd.read_csv(samples_file)

    groups = {}
    for group, df in samples.groupby("group"):
        groups[group] = df["vcf_file"].tolist()
    return groups

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True, help="CSV file with vcf_file and group columns")
    args = parser.parse_args()

    groups = load_samples(args.samples)
    print(groups)
