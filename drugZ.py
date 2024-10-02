# DRUGZ (modified version): Script for the identification of drug-gene interactions in
# paired sample genomic perturbation screens
# original version: https://github.com/hart-lab/drugz
# Last modified 19.11.2023
# ------------------------------------
import sys
import numpy as np
import pandas as pd
import scipy.stats as stats
import argparse
import logging as log

log.basicConfig(level=log.INFO)
log_ = log.getLogger(__name__)

pd.options.mode.chained_assignment = None
# default='warn' - printing a warning message
# None - ignoring the warning
# "raise" - raising an exception

# ------------------------------------
# constants
norm_value = 1e7
min_reads_thresh = 1


# ------------------------------------
# Parameters

# infile = input read counts matrix (first column is expected to contain the gene name, file is expected to be
# tab-delimited)
# drugz_output_file = output filename in which you will write the drugz results
# fc_outfile = output filename in which you will write the log2 fold-changes computed by drugz
# target_samples = the names of target samples
# reference_samples = the names of reference/control samples
# remove_genes = Genes that should be removed from the table
# unpaired = boolean, if True - compares mean(treated samples) to mean(control samples), otherwise calculates
#   fold-change individually and ultimately sums up the values for the drugz computation
# pseudocount = counts added to the observed read counts for normalization, default = 5
# half_window_size = size of the first bin and half the size of the initial sample (window) to estimate std,
# default = 500 (for whole genome screens)
# replicates = boolean, defines whether the dataset contains replicates or not
# ------------------------------------

def check_input(reads, reference_sample, target_sample, replicates, unpaired):
    """
    This function verifies the input to avoid errors that occur during the analysis.

    :param reads: The dataset with the reads
    :param reference_sample: A list with the name of the reference_sample
    :param target_sample: A list with the name of the target_sample
    :param replicates: A boolean indicating the presence of replicates in the dataframe reads
    :param unpaired: A boolean indicating whether the analysis is run in an unpaired (True) or paired (False) mode

    :return: new_reference_sample, new_target_sample: The names of the replicate columns for
    reference_sample and target_sample
    """

    # Check if the element in "target_sample" and "reference_sample" match
    # column names of the dataset when the notation after the last "_" is removed
    if replicates:

        # Store all replicate columns corresponding to the target and reference sample in a new list
        new_reference_sample = [column for column in reads.columns if reference_sample == column.rsplit('_', 1)[0]]
        new_target_sample = [column for column in reads.columns if target_sample == column.rsplit('_', 1)[0]]
    else:

        # Store the columns that correspond to the name in target_sample and reference_sample
        new_reference_sample = [column for column in reads.columns if column in reference_sample]
        new_target_sample = [column for column in reads.columns if column in target_sample]

    # Raise a ValueError if either reference_sample or target_sample could not be found
    if not new_reference_sample:
        raise ValueError(
            f"The element '{reference_sample}' from the reference_sample list cannot be found in the given dataframe.")
    if not new_target_sample:
        raise ValueError(
            f"The element '{target_sample}' from the target_sample list cannot be found in the given dataframe.")

    # Check if the lengths of the target_sample and reference_sample are the same in case unpaired is False
    if len(new_reference_sample) != len(new_target_sample) and unpaired is False:
        raise ValueError("The number of target samples and reference samples must be the same to "
                         "run the paired DrugZ analysis.")

    return new_reference_sample, new_target_sample


def load_reads(filepath, index_column, genes_to_remove):
    """
    Function to load the input file and remove genes if indicated by genes_to_remove.

    :param filepath: The path to the file to be loaded
    :param index_column: The column to use as an index (column with sgRNA names)
    :param genes_to_remove: A string of comma separated gene names that should be removed

    :return: reads: A dataframe containing the read counts per gene
    """

    # Load the input file with guide IDs as row IDs
    reads = pd.read_csv(filepath, index_col=index_column, delimiter='\t')

    # Remove genes specified in genes_to_remove
    if genes_to_remove:
        gene_column = reads.columns.values[0]
        reads = reads.loc[~reads[gene_column].isin(genes_to_remove), :]

    return reads


def normalize_readcounts(reads, target_sample, reference_sample):
    """
    Function to normalise input read counts using the global variable norm_value.

    :param reads: Dataframe containing reads counts (guide level)
    :param target_sample: List of columns names for the samples in the target group
    :param reference_sample: List of column names for the samples in the reference/control group

    :return: normalised_counts: A dataframe containing normalised read counts
    """

    # Select the treatment and reference columns
    reads_to_normalize = reads[reference_sample + target_sample]

    # Normalize raw read counts using norm_value (1e7)
    normalized_counts = (norm_value * reads_to_normalize) / reads_to_normalize.sum().values
    return normalized_counts


def calculate_fold_change(reads, normalized_counts, target_sample, reference_sample, pseudocount, replicate):
    """
    Function to create a dataframe with index as guide ids and calculate
    the log2 ratio/fold-change between treated and reference samples

    :param reads: Dataframe containing read counts
    :param normalized_counts: Dataframe containing normalized read counts
    :param target_sample: List of control sample names
    :param reference_sample: List of treated sample names
    :param pseudocount: Constant value added to all reads (default 5) - prevents log(0) problems
    :param replicate: integer, specifying the current replicate number

    :return: A dataframe with the calculated fold-change for each replicate and
    initialized columns for the guide estimated variance and z-score
    """

    fold_change = pd.DataFrame(index=reads.index.values)
    fold_change['GENE'] = reads[reads.columns.values[0]]

    # Generate fold-change, estimated_variance, and fold-change zscore column ids for each replicate
    fc_replicate_id = 'fc_{replicate}'.format(replicate=replicate)
    fc_zscore_id = 'zscore_' + fc_replicate_id
    empirical_bayes_id = 'eb_std_{replicate}'.format(replicate=replicate)

    # Get the control and treatment sample ids for each replicate
    target_replicate = target_sample[replicate]
    reference_replicate = reference_sample[replicate]

    # Add the control sample column to the fold change dataframe and sort by this column
    fold_change[reference_replicate] = reads[reference_replicate]
    fold_change.sort_values(reference_replicate, ascending=False, inplace=True)

    # Extract the number of rows (number of guides) of the reads dataframe
    no_of_guides = reads.shape[0]

    # Fill in the estimated_variance and fold-change_zscore columns with 0s for now
    fold_change[empirical_bayes_id] = np.zeros(no_of_guides)
    fold_change[fc_zscore_id] = np.zeros(no_of_guides)

    # Calculate the log2 ratio of treatment normalised read counts to control - fold-change
    fold_change[fc_replicate_id] = np.log2(
        (normalized_counts[target_replicate] + pseudocount) / (normalized_counts[reference_replicate] + pseudocount))

    return fold_change


def empirical_bayes(fold_change, half_window_size, no_of_guides, fc_replicate_id, empirical_bayes_id, fc_zscore_id):
    """
    Calculate the variation present in fold-change between treatment and control read counts for bins of the data,
    smoothing the variation during this process thus ensuring it increases or remains the same as the estimate for
    the previous bin. The estimate of variation is first calculated for 2x half_window_size guides (sorted by reads in
    corresponding control sample). The default for half_window_size is 500. So, for the first 1000 values we calculate
    the estimate of variation and the set the 1st bin (0 to half_the_window_size - i.e. 0 to 500 by default) equal to
    this estimate. For each following bin between the half-window-size and n - (where n is the number of rows (guides))
    we then calculate the estimate of variance for this bin:
        if the estimate is greater than for the previous bin, we keep this calculated estimated variance
        otherwise (estimate is less than for the previous bin) we use the estimate variance of the previous bin.
    This smooths the variance, i.e. estimated variance only ever increases or remains flat between each bin.
    The final bin (n-half_window_size : n) is then set to the variance of the previous bin (e.g. the last estimate
    to be calculated).

    :param fold_change: A dataframe containing fold-change (log2 ratio of treatment read counts to control read counts)
    :param half_window_size: An integer value equal to the size of the first bin and half the size of the initial sample
    (window) to estimate StDev. Default is 500.
    :param empirical_bayes_id: Column header under which to store the variance estimates.
    :param no_of_guides: An integer value equal to the number of rows (guides) of the data frame.
    :param fc_replicate_id: The ID that is assigned to each replicate
    :param fc_zscore_id: Column header under which the fold-change values are stored for the current comparison
     (replicate)

    :return: fold_change: Updated instance of the input dataframe
    """

    # Calculate the standard deviation of fold-change based on a 2 * define window size range
    std_dev = fold_change.iloc[0: half_window_size * 2][fc_replicate_id].std()
    fold_change[empirical_bayes_id][0: half_window_size] = std_dev

    # Iterate in a range(half_window_size, n-half_window_size, 25) where n is the number of guides
    for i in range(half_window_size, no_of_guides - half_window_size + 25, 25):
        # in every bin calculate stdev
        std_dev = fold_change.iloc[i - half_window_size:i + half_window_size][fc_replicate_id].std()

        # If the current variation is greater than the one for previous bin then set variation equal to this
        if std_dev >= fold_change[empirical_bayes_id].iloc[i - 1]:
            fold_change[empirical_bayes_id][i:i + 25] = std_dev  # set new std in whole step size (25)
        # Otherwise, set it equal to the variation of the previous bin
        # This allows variation estimate for each bin to only increase or stay the same as the previous
        else:
            fold_change[empirical_bayes_id][i:i + 25] = fold_change.iloc[i - 1][empirical_bayes_id]

    # Get the variation estimate for the final bin and set the remaining values in the empirical bayes column
    # equal to this estimate
    results = fold_change.iloc[no_of_guides - (half_window_size + 1)][empirical_bayes_id]
    fold_change[empirical_bayes_id][no_of_guides - half_window_size:] = results

    # Calculate the z_score for each guide (fc/eb_std)
    log_.debug('Calculating guide-level z-scores')
    fold_change[fc_zscore_id] = fold_change[fc_replicate_id] / fold_change[empirical_bayes_id]

    return fold_change, fc_zscore_id


def calculate_gene_drugz_score(per_gene_score, min_observations):
    """
    Calculate per gene statistics for the z-scores aggregated across all comparisons
    The summed z-scores and the number of observations for each gene are first aggregated. These z-scores are then
    normalised and p-value estimates (assuming gaussian distribution), rank position, and FDR are calculated
    The statistics are first (with normalised z-scores ranked smallest to largest) to identify synthetic
    interactions and then (with normalised z-scores now ranked largest to smallest) to identify suppressor interactions

    :param per_gene_score: Data frame containing calculated z-scores per comparison
    :param min_observations: An integer value to act as a threshold for the minimum number observations to be included
    in the analysis (default=1)

    :return: per_gene_results: A dataframe of summary statistics for each gene
    """

    # Get a dataframe of values for genes where the number of observations is greater than the minimum threshold
    per_gene_results = per_gene_score.loc[per_gene_score.numObs >= min_observations, :]

    # Update the row number (number of genes) for this new dataframe.
    no_of_genes = per_gene_results.shape[0]

    # Calculate normalized gene z-score by:
    # 1. normalizing the sumZ values by number of observations,
    # 2. re-normalizing these values to fit uniform distribution of null p-vals
    normalised_z_scores = stats.zscore(per_gene_results['sumZ'] / np.sqrt(per_gene_results['numObs']))
    per_gene_results['normZ'] = normalised_z_scores

    # Sort the data frame by normZ (ascending) to identify synthetic interactions
    per_gene_results.sort_values('normZ', ascending=True, inplace=True)

    # Calculate p-values (from normal dist), and fdrs (by benjamini & hochberg correction)
    per_gene_results['pval_synth'] = stats.norm.sf(per_gene_results['normZ'] * -1)
    per_gene_results['rank_synth'] = np.arange(1, no_of_genes + 1)
    scale = per_gene_results['rank_synth'] / float(no_of_genes)
    per_gene_results['fdr_synth'] = per_gene_results['pval_synth'] / scale
    per_gene_results['fdr_synth'] = np.minimum.accumulate(per_gene_results['fdr_synth'][::-1])[::-1]

    # Resort by normZ (descending) and recalculate above values to identify suppressor interactions
    per_gene_results = per_gene_results.sort_values('normZ', ascending=False)
    per_gene_results['pval_supp'] = stats.norm.sf(per_gene_results['normZ'])
    per_gene_results['rank_supp'] = np.arange(1, no_of_genes + 1)
    scale = per_gene_results['rank_supp'] / float(no_of_genes)
    per_gene_results['fdr_supp'] = per_gene_results['pval_supp'] / scale
    per_gene_results['fdr_supp'] = np.minimum.accumulate(per_gene_results['fdr_supp'][::-1])[::-1]

    per_gene_results = per_gene_results.sort_values('normZ', ascending=True)

    return per_gene_results


def write_drugz_output(outfile, output):
    """
    Function that writes DrugZ results to a file.

    :param outfile: Dataframe with DrugZ results
    :param output: Per gene calculated statistics

    :return: fout: Output file
    """
    fout = outfile
    if not hasattr(fout, 'write'):
        fout = open(fout, 'w')
    fout.write('GENE')
    cols = output.columns.values
    for c in cols:
        fout.write('\t' + c)
    fout.write('\n')

    for i in output.index.values:
        fout.write('{0:s}\t{1:3.2f}\t{2:d}\t{3:4.2f}\t{4:.3g}\t{5:d}\t{6:.3g}\t{7:.3g}\t{8:d}\t{9:.3g}\n'.format(
            i, output.loc[i, 'sumZ'], int(output.loc[i, 'numObs']),
            output.loc[i, 'normZ'], output.loc[i, 'pval_synth'],
            int(output.loc[i, 'rank_synth']), output.loc[i, 'fdr_synth'],
            output.loc[i, 'pval_supp'], int(output.loc[i, 'rank_supp']),
            output.loc[i, 'fdr_supp']))
    fout.close()

    return fout


def get_args():
    """ Parse user giver arguments"""

    parser = argparse.ArgumentParser(description='DrugZ for chemogenetic interaction screens',
                                     epilog='dependencies: pylab, pandas')
    parser._optionals.title = "Options"
    parser.add_argument("-i", dest="infile", type=argparse.FileType('r'), metavar="sgRNA_count.txt",
                        help="sgRNA read count file", default=sys.stdin)
    parser.add_argument("-o", dest="drugz_output_file", type=argparse.FileType('w'), metavar="drugz-output.txt",
                        help="drugz output file", default=sys.stdout)
    parser.add_argument("-f", dest="fc_outfile", type=argparse.FileType('w'), metavar="drugz fold-change.txt",
                        help="drugz normalized fold-change file (optional")
    parser.add_argument("-c", dest="reference_sample", metavar="reference sample", required=True,
                        help="reference samples, comma delimited")
    parser.add_argument("-x", dest="target_sample", metavar="target sample", required=True,
                        help="target samples, comma delimited")
    parser.add_argument("-r", dest="remove_genes", metavar="remove genes", help="genes to remove, comma delimited")
    parser.add_argument("-p", dest="pseudocount", type=int, metavar="pseudocount", help="pseudocount (default=5)",
                        default=5)
    parser.add_argument("-I", dest="index_column", type=int,
                        help="Index column in the input file (default=0; GENE_CLONE column)", default=0)
    parser.add_argument("--minobs", dest="minObs", type=int, metavar="minObs", help="min number of obs (default=1)",
                        default=1)
    parser.add_argument("--half_window_size", dest="half_window_size", type=int, metavar="half_window_size",
                        help="width of variance-estimation window", default=500)
    parser.add_argument("-unpaired", dest="unpaired", action='store_true', default=False,
                        help='comparison status, paired (default) or unpaired')
    parser.add_argument("-replicates", dest="replicates", action='store_true', default=False,
                        help='replicates in the input file, True or False (default)')
    return parser.parse_args()


def calculate_unpaired_foldchange(reads, normalized_counts, target_sample, reference_sample, pseudocount):
    """
    Function that calculates the unpaired fold-change for each guide.

    :param reads: Dataframe with the raw read counts
    :param normalized_counts: Dataframe with the raw read counts
    :param target_sample: Name of the target sample column
    :param reference_sample: Name of the reference sample column
    :param pseudocount: Constant value added to all reads (default 5) - prevents log(0) problems

    :return: fold_change: Log2 fold-change between target and reference sample
    """

    fold_change = pd.DataFrame(index=reads.index.values)
    fold_change['GENE'] = reads[reads.columns.values[0]]

    # Add the control samples mean read counts column to the fold change dataframe and sort by this column
    fold_change['ctrl_mean_reads'] = reads[reference_sample].mean(axis=1)
    fold_change.sort_values('ctrl_mean_reads', ascending=False, inplace=True)

    # Add the column for unpaired fold-change (i.e. mean fold-change)
    fold_change['mean_fc'] = np.log2((normalized_counts[target_sample].mean(axis=1) + pseudocount) /
                                     (normalized_counts[reference_sample].mean(axis=1) + pseudocount))

    # set empty columns for eb variance and z-scores
    fold_change['eb_std'] = np.zeros(normalized_counts.shape[0])
    fold_change['zscore'] = np.zeros(normalized_counts.shape[0])

    return fold_change


def drugz_analysis(args):
    """ Call all functions and run DrugZ analysis

    :param args: User given arguments
    :return gene_normZ: normalized read counts
    """

    log_.debug("Initiating analysis")

    target_sample = str(args.target_sample)
    reference_sample = str(args.reference_sample)

    if args.remove_genes is None:
        remove_genes = []
    else:
        remove_genes = args.remove_genes.split(',')

    log_.debug("Target sample: " + target_sample)
    log_.debug("Reference sample: " + reference_sample)

    log_.debug("Loading the read count matrix")
    reads = load_reads(filepath=args.infile, index_column=0, genes_to_remove=remove_genes)
    no_of_guides = reads.shape[0]

    reference_sample, target_sample = check_input(reads, reference_sample=reference_sample, target_sample=target_sample,
                                                  replicates=args.replicates,
                                                  unpaired=args.unpaired)

    log_.debug("Normalizing read counts")
    normalized_counts = normalize_readcounts(reads=reads, target_sample=target_sample,
                                             reference_sample=reference_sample)

    num_replicates = len(reference_sample)

    if args.unpaired:
        log_.debug('Calculating raw fold change - unpaired approach')
        fold_change = calculate_unpaired_foldchange(reads, normalized_counts,
                                                    target_sample=target_sample,
                                                    reference_sample=reference_sample, pseudocount=args.pseudocount)

        log_.debug('Calculating smoothed Empirical Bayes estimates of stdev')
        fold_change = empirical_bayes(fold_change=fold_change, half_window_size=args.half_window_size,
                                      no_of_guides=no_of_guides, fc_replicate_id='mean_fc',
                                      empirical_bayes_id='eb_std', fc_zscore_id='zscore')[0]

        gRNA_stds = pd.DataFrame(
            fold_change.groupby('GENE')['zscore'].apply(lambda x: pd.Series([np.nanstd(x.values),
                                                                             np.count_nonzero(x.values)]))).unstack()
        gRNA_stds.columns = ['std sgRNAs/gene', 'numObs']
        gRNA_stds = gRNA_stds.reset_index()
        gRNA_names = fold_change.index
        fold_change = fold_change.merge(gRNA_stds, how="left", on="GENE")
        fold_change.set_index(gRNA_names, inplace=True)

        if args.gRNA_outfile:
            fold_change.to_csv(args.gRNA_outfile, sep='\t', float_format='%4.3f')

        # Sum up the fold-changes from sgRNAs for each gene
        per_gene_scores = pd.DataFrame(
            fold_change.groupby('GENE')['zscore'].apply(lambda x: pd.Series([np.nansum(x.values),
                                                                             np.count_nonzero(x.values)]))).unstack()
        per_gene_scores.columns = ['sumZ', 'numObs']

        gene_normZ = calculate_gene_drugz_score(per_gene_score=per_gene_scores, min_observations=1)

        log_.debug('Writing output file unpaired results')
        write_drugz_output(outfile=args.drugz_output_file, output=gene_normZ)

        log_.debug('Finished the analysis.')

        return gene_normZ

    else:
        fc_zscore_ids = list()
        fold_changes = list()
        for i in range(num_replicates):
            log_.debug('Calculating raw fold change for replicate {0}'.format(i + 1))
            fold_change = calculate_fold_change(reads, normalized_counts,
                                                target_sample=target_sample,
                                                reference_sample=reference_sample,
                                                pseudocount=args.pseudocount, replicate=i)

            log_.debug('Calculating smoothed Empirical Bayes estimates of stdev for replicate {0}'.format(i + 1))
            fold_change, fc_zscore_id = empirical_bayes(fold_change=fold_change, half_window_size=args.half_window_size,
                                                        no_of_guides=no_of_guides,
                                                        fc_replicate_id='fc_{replicate}'.format(replicate=i),
                                                        empirical_bayes_id='eb_std_{replicate}'.format(replicate=i),
                                                        fc_zscore_id='zscore_fc_{replicate}'.format(replicate=i))

            fold_changes.append(fold_change)

            fc_zscore_ids.append(fc_zscore_id)

        fold_changes = pd.concat(fold_changes, axis=1, sort=False)
        fold_changes = fold_changes.loc[:, ~fold_changes.columns.duplicated()]

        if args.gRNA_outfile:
            fold_changes.to_csv(args.gRNA_outfile, sep='\t', float_format='%4.3f')

        log_.debug('Calculating gene-level z-scores')
        # This is going to produce a per gene summation of the z-scores for each comparison. Missing values are
        # converted to zero. The aggregated zscore will be stored in a column named sumZ.
        # The second column will be the number of all non_zero observations (non-zero guides/gene) for each gene.
        per_gene_scores = fold_changes.groupby('GENE')[fc_zscore_ids].apply(lambda x: pd.Series([np.nansum(x.values),
                                                                                                 np.count_nonzero(
                                                                                                     x.values)]))
        per_gene_scores.columns = ['sumZ', 'numObs']
        gene_normZ = calculate_gene_drugz_score(per_gene_score=per_gene_scores, min_observations=1)

        log_.debug('Writing output file paired results')
        write_drugz_output(outfile=args.drugz_output_file, output=gene_normZ)

        log_.debug('Finished the analysis.')
        return gene_normZ
