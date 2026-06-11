#!/urs.bin/env python3
import pandas as pd
import biom
from biom.table import Table
import argparse

def tsv_to_biom(input_tsv, output_biom):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(input_tsv, sep='\t', index_col=0)
    # Convert the DataFrame to a biom.Table object and make sure both observation and sample IDs are strings
    table = Table(df.values, 
              observation_ids=[str(i) for i in df.index], 
              sample_ids=[str(i) for i in df.columns])
    # Write the biom.Table object to a BIOM file
    with biom.util.biom_open(output_biom, 'w') as f:
        table.to_hdf5(f, "Converted from TSV")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_tsv", required=True, help="Input TSV file")
    parser.add_argument("--output_biom", required=True, help="Output BIOM file")
    args = parser.parse_args()
    tsv_to_biom(args.input_tsv, args.output_biom)