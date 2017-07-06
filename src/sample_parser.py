import pandas, re, os, numpy, time, sys

# Parse the ground truth file into a DataFrame
def parse_gt(excel_file):
    # Create data frame from Excel file
    df = pandas.read_excel(
            excel_file, 'Sheet2', index_col=None, na_values=['NA']
    )
    # These two columns are the only ones we care about
    ids = df['SAMPLE ID']
    variants = df['"REPORTABLE" VARIANTS CHECKED in PPMP 0.3.1']

    sample_ids = []
    genes = []
    dna_changes = []
    protein_changes = []

    for index, sample in enumerate(variants):
        sample_variants = sample.split('\n')
        sample_id = ids[index]
        # Split each line of the last column into gene, dna, and protein info
        for variant in sample_variants:
            if (not re.search('^NOTE', variant)
                    and not variant == ''
                    and not re.search('^Missed', variant)):
                genes.append(variant.split(' ')[0])
                sample_ids.append(sample_id)
                # Get DNA change
                try:
                    dna_change = re.search('c\..*[A-Z]:', variant).group(0)\
                            .split(':')[0]
                except AttributeError:
                    dna_change = ''
                dna_changes.append(dna_change)
                # Get protein change
                try:
                    protein_change = re.search('p\..*$', variant).group(0)
                except AttributeError:
                    protein_change = ''
                protein_changes.append(protein_change)

    # Convert all sample IDs to strings (some are ints)
    sample_ids = [str(_id) for _id in sample_ids]
    output = pandas.DataFrame({ 'SAMPLE_ID'      : sample_ids,
                                'GENE'           : genes,
                                'DNA_CHANGE'     : dna_changes,
                                'PROTEIN_CHANGE' : protein_changes })
    """
    Example output:
            DNA_CHANGE    GENE  PROTEIN_CHANGE   SAMPLE_ID
    0         c.364C>T  CDKN2A       p.Arg122*  1530804642
    1         c.406G>A   CRLF2     p.Val136Met  1530804642
    2        c.1682G>A     NF1       p.Trp561*  1530804642
    ...
    6         c.740A>G    ABL1     p.Lys247Arg  1601912547
    7  c.384_386dupGCA  ARID1B     p.Gln129dup  1601912547
    ...
    """
    return output

def parse_bed(bed):
    parsed = pandas.read_csv(
            bed, names=['CHROMOSOME', 'START', 'END', 'GENE'], sep='\t'
    )
    genes = [
            re.search('_[A-Z0-9]*_', gene_string).group(0).split('_')[1]
            for gene_string in parsed['GENE']
    ]
    parsed['GENE'] = genes
    return parsed

# TODO does the 'end' position include the last covered position?
# Returns a DataFrame of individual positions, genes, and chromosomes covered
def split_panel(panel):
    positions = []
    genes = []
    chromosomes = []
    for i in range(0, panel.shape[0]):
        pos_list = list(range(panel['START'][i], panel['END'][i]))
        positions += pos_list
        genes += len(pos_list) * [panel['GENE'][i]]
        chromosomes += len(pos_list) * [panel['CHROMOSOME'][i]]
    return pandas.DataFrame({
            'POSITION' : positions,
            'GENE' : genes,
            'CHROMOSOME' : chromosomes
    }).reset_index(drop=True)

# Takes a list of sample IDs and ground truth DataFrame
# These will be slightly different, so this method matches IDs from the first
# list to names from the second.
# Returns a list of tuples, with the first string in the tuple being the sample
# folder name, and the second being the name in the ground truth DataFrame
def find_matches(sample_ids, gt_parsed):
    # Initialize dictionary, with
    # first item = name of sample folder,
    # second item = ID in ground truth
    match_list = []
    # Search for similar ID in ground truth
    candidate_ids = set(gt_parsed['SAMPLE_ID']) # List of IDs in ground truth
    for sample_id in sample_ids:
        sample = sample_id.lower().replace(' ','').split('_')[-1]
        for candidate_id in candidate_ids:
            candidate = candidate_id.lower().replace(' ','')
            candidate_list = candidate_id.lower().split(' ')
            # Check if candidate and sample IDs are the same,
            # or if all the parts of the candidate ID are present in the sample
            if candidate in sample or all(x in sample for x in candidate_list):
                match_list.append((sample_id, candidate_id))
    return match_list

# Splits ground truth parsed data frame into separate data frames for each ID
# Only generates those for which there are target_transcripts files
def split_gt(gt_parsed, match_list):
    # List of data frames
    gt_df_list = []
    # List of unique IDs in the ground truth
    unique_ids = [sample[1] for sample in match_list]

    for unique_id in unique_ids:
        # Find all indices with this ID
        indices = gt_parsed[gt_parsed['SAMPLE_ID'] == unique_id].index.tolist()
        # Create a data frame
        sample_gt = pandas.DataFrame(gt_parsed[indices[0]:indices[-1]+1])
        gt_df_list.append(sample_gt.reset_index(drop=True))

    return gt_df_list

# Returns the names of the variant callers
def get_caller_names(sample):
    headers = list(sample)
    return [h for h in headers if re.search('^GT_', h)]

# Returns a list of DataFrames of samples for which matches exist in the ground
# truth
def get_samples(samples_dir, match_list):
    # List of sample DataFrames
    df_list = []
    for sample in match_list:
        # Get IDs
        sample_id = sample[0]
        gt_id = sample[1]
        # Get sample file path
        sample_path = os.path.join(samples_dir, sample_id)
        tab_path = os.path.join(
                sample_path, next(os.walk(sample_path))[1][0],
                'vDEMO1', 'r1', 'hg19', 'panel', 'ann.target_transcripts.tab'
        )
        # Import sample CSV
        sample_df = pandas.read_csv(tab_path, sep='\t')
        # Get variant caller names
        caller_names = get_caller_names(sample_df)
        # List of columns to keep
        col_list = [
                'CHROMOSOME', 'POSITION', 'TOTAL_CALLERS',
                'SNPEFF_ANNOTATED_GENE', 'SNPEFF_ANNOTATED_DNA_CHANGE',
                'SNPEFF_ANNOTATED_PROTEIN_CHANGE', 'EXAC_FREQ_ESTIMATE',
                '1KG_FREQ_ESTIMATE', 'TARGET_ALLELE_COUNT',
        ] + caller_names
        # New sample DataFrame with relevant columns only
        altered_df = sample_df[col_list].rename(columns={
                'SNPEFF_ANNOTATED_GENE': 'GENE',
                'SNPEFF_ANNOTATED_DNA_CHANGE': 'DNA_CHANGE',
                'SNPEFF_ANNOTATED_PROTEIN_CHANGE': 'PROTEIN_CHANGE'
        })
        # Add SAMPLE_ID column (for when this DataFrame will be merged with others)
        altered_df['SAMPLE_ID'] = sample_df.shape[0] * [gt_id]
        # Add sample DataFrame to list
        df_list.append(altered_df.reset_index(drop=True))
    return df_list

def combine(gt_file, bed_file, samples_dir):
    # Parse ground truth Excel file
    gt_parsed = parse_gt(gt_file)
    # Parse .bed file
    bed = parse_bed(bed_file)
    # Find sample ID matches
    match_list = find_matches(next(os.walk(samples_dir))[1], gt_parsed)
    # Combine ground truth DataFrames
    gt = pandas.concat(split_gt(gt_parsed, match_list)).reset_index(drop=True)
    # List of sample DataFrames
    samples_list = get_samples(samples_dir, match_list)
    # List of variant caller names
    callers = get_caller_names(samples_list[0])

    # Combine all sample DataFrames into one big DataFrame
    df = pandas.concat(samples_list).reset_index(drop=True)
    # List of variants covered by the small panel
    covered_list = df.shape[0] * [False]
    # List of variants reported in the ground truth
    reportables_list = df.shape[0] * [False]

    # Set up progress indicator
    sys.stdout.write('Generating combined DataFrame:  0%')
    sys.stdout.flush()

    # Cycle through every variant
    for s in range(0, df.shape[0]):
        # Progress indicator
        progress = int(100 * (s / df.shape[0]))
        if progress < 10:
            sys.stdout.write('\b\b' + str(progress) + '%')
        else:
            sys.stdout.write('\b\b\b' + str(progress) + '%')
        sys.stdout.flush()

        # Find variants covered by the small panel
        # If the gene is covered by the small panel
        if df['GENE'][s] in set(bed['GENE']):
            covered_genes = [
                    x for x, gene in enumerate(bed['GENE'])
                    if gene == df['GENE'][s]
            ]
            # List of possible locations in small panel coverage file
            for n in covered_genes:
                # If variant position is in a location covered by the small panel
                if (df['CHROMOSOME'][s] == bed['CHROMOSOME'][n]
                        and int(df['POSITION'][s]) >= bed['START'][n]
                        and int(df['POSITION'][s]) <= bed['END'][n]):
                # if (df['CHROMOSOME'][s] == bed['CHROMOSOME'][n]
                #         and int(df['POSITION'][s]) == bed['POSITION'][n]):
                    covered_list[s] = True

        # Find variants reported in the ground truth
        # This if statement weeds out unnecessary checks
        if (df['GENE'][s] in gt['GENE'].tolist()
                and df['DNA_CHANGE'][s] in gt['DNA_CHANGE'].tolist()
                and df['PROTEIN_CHANGE'][s] in gt['PROTEIN_CHANGE'].tolist()):
            # Cycle through variants in ground truth
            for g in range(0, gt.shape[0]):
                # If the variant info matches between the DataFrame and ground truth
                if (df['SAMPLE_ID'][s] == gt['SAMPLE_ID'][g]
                        and df['GENE'][s] == gt['GENE'][g]
                        and df['DNA_CHANGE'][s] == gt['DNA_CHANGE'][g]
                        and df['PROTEIN_CHANGE'][s] == gt['PROTEIN_CHANGE'][g]):
                    reportables_list[s] = True

    sys.stdout.write('\b\b\b100%')
    sys.stdout.flush

    # Add covered and reportable lists to DataFrame
    df['COVERED'] = covered_list
    df['REPORTABLE'] = reportables_list
    # Create combined tab file
    df.to_csv(
            os.path.join(samples_dir, 'combined.tab'), sep='\t',
            encoding='utf-8', index=False
    )
    print('\nOutput to file ' + os.path.join(samples_dir, 'combined.tab'))

if __name__ == "__main__":
    import argparse

    # Parser / command line instructions
    # Description of script
    arg_parser = argparse.ArgumentParser(
            description='Compares sample variants to\
            ground truth file and produces false negatives, false positives,\
            and true positives.'
    )
    # First argument: ground truth Excel file
    arg_parser.add_argument('file1', help='ground truth .xlsx', action='store')
    # Second argument: small panel gene list, bed file
    arg_parser.add_argument('file2', help='.bed', action='store')
    # Third argument: directory containing sample directories
    arg_parser.add_argument(
            'directory', help='directory containing all sample directories',
            action='store'
    )

    # Parse arguments
    args = arg_parser.parse_args()
    gt_file = args.file1
    bed_file = args.file2
    samples_dir = args.directory

    combine(gt_file, bed_file, samples_dir)
