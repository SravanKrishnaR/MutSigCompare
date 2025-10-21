import pysam
from Bio.Seq import Seq
import pandas as pd
import itertools
import argparse
import os
from Data_Importer import load_samples
from PCA_plotter import plot_pca
from COSMIC_NNLS import cosmic_nnlss


def SNP_parsing(vcf_files,group_name,ref_file):
    ref = pysam.FastaFile(ref_file) 

    for vcf_file in vcf_files:   # loop through multiple input files
        vcf = pysam.VariantFile(vcf_file)   # open as VariantFile    
        count = {}

        #iterate through the vcf file
        for rec in vcf.fetch(): 
            chrom = rec.chrom
            pos = rec.pos
            
            #check if the length of ref and it's alt is single base 
            if len(rec.ref) == 1 and all(len(alt) == 1 for alt in rec.alts):
                #iterate through the alts(if multi-allelic)
                for alt in rec.alts: 

                    #Normalization with pyramidines(C,T)
                    if rec.ref in ["G", "A"]:
                        #fetch the trinucleotide and reverse complement it
                        trinuc = ref.fetch(chrom, pos-2, pos+1).upper()
                        rc_trinuc = str(Seq(trinuc).reverse_complement())
                    
                        #reverse complement the reference and the alt if purines
                        rc_ref=str(Seq(rec.ref).reverse_complement())
                        rc_alt=str(Seq(alt).reverse_complement())

                        mutation = f"{rc_trinuc[0]}[{rc_ref}>{rc_alt}]{rc_trinuc[2]}"
                        count[mutation] = count.get(mutation, 0) + 1

                    elif rec.ref in ["C", "T"]:
                        trinuc = ref.fetch(chrom, pos-2, pos+1).upper()
                     
                        mutation = f"{trinuc[0]}[{rec.ref}>{alt}]{trinuc[2]}"
                        count[mutation] = count.get(mutation, 0) + 1
                        
        #convert the counts dictionary to pands dataframe
        mutation_counts = pd.DataFrame(count.items(), columns=["Type", "count"])

        #all possible 96 mutations
        substitutions = ["C>A","C>G","C>T","T>A","T>C","T>G"]
        bases = ["A","C","G","T"]

        rows = []

        for substitution in substitutions:
            for lb in bases:
                for rb in bases:
                    tris = f"{lb}[{substitution}]{rb}"
                    rows.append([tris, 0])    # store as [mutation, count]

        all_mutations = pd.DataFrame(rows, columns=["Type", "count"])

        #Merge the possible 96 mutations dataframe all_mutations with the counts dataframe mutation_counts
        new_data = pd.merge(all_mutations, mutation_counts, on="Type", how="left")
        new_data["count"] = new_data["count_y"].fillna(new_data["count_x"])
        new_data = new_data[["Type", "count"]]

        #convert to csv file
        out_file = os.path.splitext(os.path.basename(vcf_file))[0] + "_mutation_counts.txt"
        new_data.to_csv(out_file, index=False)
        print(f"Processed {vcf_file} → {out_file}")

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="To input the data")
    parser.add_argument("--samples", required=True,type=str,help="Input the samples CSV file")
    parser.add_argument("--ref", required=True,type=str,help="Input the reference sequence")
    args=parser.parse_args()

    groups = load_samples(args.samples)
    group_matrices = {}

    for group_name, vcf_list in groups.items():
        SNP_parsing(vcf_list, group_name, args.ref)

        dfs = []
        for vcf_file in vcf_list:
            counts_file = os.path.splitext(vcf_file)[0] + "_mutation_counts.txt"
            if not os.path.exists(counts_file):
                raise FileNotFoundError(f"Expected counts file: {counts_file}")
            df = pd.read_csv(counts_file)
            df = df.set_index("Type")
            sample_col = f"{group_name}_{os.path.splitext(os.path.basename(vcf_file))[0]}"
            dfs.append(df["count"].rename(sample_col))

        group_matrix = pd.concat(dfs, axis=1).fillna(0).sort_index()
        group_matrices[group_name] = group_matrix

    full_matrix = pd.concat(group_matrices.values(), axis=1).fillna(0)
    matrix_rel = full_matrix.div(full_matrix.sum(axis=0), axis=1)

    sample_groups = [g for g, df in group_matrices.items() for _ in range(df.shape[1])]

    # For PCA use normalized data:
    plot_pca(matrix_rel, sample_groups, output_file="mutational_PCA.jpeg")
    
    cosmic_data="DATA/COSMIC_v3.4_SBS_GRCh37.txt"

    cosmic_nnlss(cosmic_data, matrix_rel)

