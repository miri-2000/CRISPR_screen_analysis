# Steps for the CRISPR screen analysis

## Windows Setup

1. Download [Anaconda](https://www.anaconda.com/download) (conda: 4.13.0)

2. Open the Anaconda prompt and change into a folder location of your choice to store the program

3. Paste the following into the prompt

   ```cmd
    git clone https://github.com/miri-2000/CRISPR_screen_analysis.git
    cd CRISPR_screen_analysis/scripts
    conda env create -f environment.yml
    conda activate crispr_env
    Rscript install_r_packages.R
    python launch_app.py
   ```
   You should now see a window that allows you to type in the CRISPR screen parameters

4. Fill out the parameters prompted in the window.
   In case you would rather use a file to directly write in your parameters, adapt the "file_based_launch.py"
   with your parameters and instead of running:
   ```cmd
   python launch_app.py
   ```
   run:
   ```cmd
   python file_based_launch.py
   ```

5. If analysis went well, the run will end with "Analysis of the given dataset completed"

## Example run

* The example directory contains the data from a CRISPR screen performed at the Netherlands Cancer Institute (Publication can be found [here](https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fwww.biorxiv.org%2Fcontent%2F10.1101%2F2023.08.17.553681v1&data=05%7C02%7Cmiriam.mueller%40ru.nl%7C1e2c665e10cd45619ea908dce6f4554d%7C084578d9400d4a5aa7c7e76ca47af400%7C1%7C0%7C638639184354131836%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C0%7C%7C%7C&sdata=7Xeqvu5LtzmaMkibjqBgYJn3XG3exG5ulVol59X8fiA%3D&reserved=0))

* In order to run the example after performing the set up, run the command with the default parameters:
   run:
   ```cmd
   python file_based_launch.py
   ```
* The results will be stored in "example/results"

# Expected Input

The program is expected to run with different CRISPR screens, considering that the read counts and library
file do not vary a lot when using the same screening library and sequencing facility. In case you run into problems
while running the program, please verify that your input files fulfill the requirements stated below:

**How the required additional files should look like**

- Library file:
    - needs to have two columns named “Target Gene Symbol” and “sgRNA Target Sequence” containing the gene name and
      guide
      sequence respectively
    - Non targeting controls need to have the content "Non-Targeting_Control" or "Non-Targeting Control" in the column
      “Target Gene Symbol”
- Essential/Non-essential file:
    - Gene names cannot end with a dash ("-")

- read_counts input file:
    - The first column contains the sgRNA names, the second one the sgRNA sequences
        - both of these two columns must contain unique values
        - the sgRNA names are expected to contain the name of the target gene followed by a dash ("-") (e.g. HTT-1)
    - All annotation columns contain non-numerical data while data columns only have
      numerical values
    - The file must contain a row called "nohit_row" (specifying reads that could not be mapped to any guide)
    - Rows with non-targeting controls are expected to be called "Non-Targeting_Control" or "Non-Targeting Control"
    - The input file is expected to have columns that start with "guide_mm1_ ": columns with 1bp mismatch on the sgRNA (
      any possible 1bp mismatch, and if only possible from this sgRNA) -> required for creating counts_summary file

- Gene names are expected to have the same casing in all input files (essential,- non-essential, library and read_counts
  file)

**How to call your conditions in the input file**

In order to ensure that no naming issues can cause an error, a standard naming convention has been defined. This name
consists of: *Timepoint*\_*culture-details*\_*treatment\_replicate*

1) <u>Timepoint</u>: "t0", "t1",... (upper and lowercase possible)
2) <u>Details of the cell culture</u>: "3d2d", ”3d”,”2d”,"3d_to_2d", … (upper and lowercase possible)
3) <u>Treatment</u> (optional): "tr","ut", "untreated","treated",... (upper and lowercase possible)
4) <u>Replicate information</u>: "\_r" or "\_rep" followed by a digit (upper and lowercase possible)

In case the condition name contains a barcode after the condition name, separated by a colon (":"), it will be
automatically removed. Valid condition names are: "T1_3D_to_2D_Rep3:TGGTCA", "t1_2d_Treated_rep26",...

By default, the program will convert the given name to the shortest of these examples: eg "t1_3d2d_r3", "
t1_2d_tr_r26",...

Be aware that aberrant naming may cause problems when executing the program. The program tries to catch small changes (
e.g. the treatment and details of the cell culture are mixed) but, if possible please adhere to this naming convention.

---

# Default settings for the analysis

## Output file naming

By default, the program assumes that the dataset identifier is the first part of the name of the read_counts input file.
For instance, if your read counts file is called "1234_CRISPR_screen.txt", the program will extract the file name before
the first underscore and use it as a prefix for the output files (e.g. "1234_counts_summary.csv").

## Input parameters

**essential_genes**

- **Definition:** csv file that contains list of essential genes in the screen library
- will be denoted in output files and plots as type="p"

**non_essential_genes**

- **Definition:** csv file that contains list of non-essential genes in the screen library
- will be denoted in output files and plots as type="n"

**library_file**

- **Definition:** txt file that contains a list of all guides in the screen library

**input_file**

- **Definition:** absolute filepath to the txt file with the read counts from the screen
- the first part of the filename (before the first underscore) will be used as the prefix in all result files (e.g.
  7234_Brunello_screen_read_counts.txt will have result files called '7234_read_counts.txt',...)

**replicate_type**

- **Definition:** used as basis to create correlation plots (if type="biological", analysis is done assuming replicates
  are biological replicates; vice versa for "technical")

**threshold_reads**

- **Definition:** minimum sum of read counts for one gRNA over all conditions (default: 0)
- **Example:** if threshold=3 and there are two conditions that each have one read count for guide A, this guide will be
  excluded

**unwanted_columns**

- **Definition:** used to specify column names that should be removed

- **Default values for Brunello screens:** "guide_mm1_mismatch1", "mismatch1_", "nohit_cols", "guide_mm1_nohit"
    - The mentioned columns are columns that are created in case a read matches the reads from the same gRNA of a
      certain condition (e.g. t0_r1) with a 1 base pair mismatch

- **Example:** column “t1_r2_mismatch1_t0_r1” will be removed as it contains “mismatch1_”

**unwanted_rows**

- **Definition:** used to specify row names that should be removed

- **Default values for Brunello screens:** all row names are expected to be important for the screen analysis

**unwanted_row_substrings**

- **Definition:** used to specify any additional text that should be removed from the column with the sgRNA names
- **Default values for Brunello screens:** “:mismatch”
    - **Rationale:** Since the guides can be very similar with respect to the read sequences, sometimes a read can
      appear to be a perfect match for one sgRNA while also matching another sgRNA if 1 base pair is changed. Therefore,
      some rows can end up looking like “USP17L7-1:mismatch1_USP17L19-1”. This indicates that the read matches USP17L7-1
      perfectly while it has a 1 base pair mismatch for USP17L19-1. Therefore, as by default it is better to only
      consider perfect matches, everything after “:mismatch” will be removed (USP17L7-1:mismatch1_USP17L19-1->
      USP17L7-1).

**target_samples**

- **Definition:** the names of the target conditions

**reference_samples**

- **Definition:** the names of the reference/control samples

- ---

# Explanation of output files and figures

## Summary file

This file serves as a quality metric that provides information about the distribution of all reads over each condition.
The reads are grouped into three classes. The first class describes the number of reads that have been mapped to a
certain guide in a certain condition. The second class describes the reads that can be mapped to a certain condition if
one mismatch in the read sequence is allowed. The third class indicates the sum of all reads that could not be mapped to
any guide in any condition.

Be aware that a good screen has a mapping ratio of at least 65%, meaning that at least 65% of all reads are part of to
class 1. If this is not the case, either because many reads either contain one mismatch or cannot be mapped to a
condition at all, the screen should be repeated.

## Data_prep file

This file contains the cleaned version of the read counts file, where all unnecessary rows and columns have been
removed. The file contains the sgRNA and gene names as well as the type of gene (e.g. = target, essential,
non-essential, non-targeting) and any additional annotation columns that were given in the read counts file + the read
counts for each guide per condition.

## Sizefactors file

This file shows the R<sup>2</sup>-value of each condition to indicate how good the normalization method fits the data.
Two normalization methods are compared: the median of ratios sizefactor calculation (used by DeSeq2) and the calculation
of sizefactors using the geometric mean alone .

In the past DeSeq2 was a very good way to compute sizefactors that are used as part of the normalization. Nowadays, the
computation using the geometric mean alone has been found to work better. To compare which of the two techniques are
better for the current plot, the R<sup>2</sup>-values of the two sizefactors are compared for each condition (e.g.
t2_2d_tr,...).

As the default normalization method is based on the total, this file should be checked to see whether the total method
shows R<sup>2</sup>-values that are closer to 1 (optimal goodness of fit) than DeSeq2's median of ratios method. In case
the median ratio method has values that are closer to 1 than the total-based method, change *line 521* in *R_analysis_1*
accordingly.

```R
> updates=normalizeData(reads,Conditions,sizefactor_file,normalized_file,sizeFactorsBasedOnAll=F)
```

In case both methods perform well, it might be worth executing the screen twice with both normalization methods.

## Normalized file and DrugZ input file

These files contain the normalized read counts for each guide and each condition. The only difference between both of
them is that the normalized file contains all annotation columns that were given and created while running the program.
The DrugZ input file only contains the columns that are actually required to run DrugZ.

## Barplot

This figure shows the sum of all reads for each condition before and after the normalization step to verify that the
normalization lead to the expected result (every condition has the same sum of reads). In case the normalized plot shows
differences in the sum, the normalization step should be checked.

## Rankplot

This figure displays the distribution of the sum of reads per guide. The x-axis indicates the rank (the index of each
guide when sorting the guides by the sum of their reads (displayed in %)), while the y-axis shows the sum of the reads
per guide on a log10 scale (e.g. log10 of 2 = 100 reads). In the title, the percentage of guides that have a sum of reads
above a log10 of 2 is displayed to indicate whether there are big differences between different conditions.

## Correlation plot

This figure compares the distribution of read sums per guide between the replicates of each condition in a pairwise
manner. The read sums are displayed on a log10 scale (e.g. log10 of 2 = 100 reads) and the classes (target, essential,
non-essential, non-targeting) are shown by different colours. The R<sup>2</sup>-value of the pairwise comparison between
the replicates of each condition can be seen in above the plot.

## Correlation file

This file shows the R<sup>2</sup>-value and the pearson correlation coefficient of the pairwise comparison between the
replicates of each condition. This is an additional file for the correlation plot that shows the results from the
comparison plot (plus the pearson correlation) in tabular form.

## Heatmap

This figure shows the heatmap of the replicates of each condition. Darker fields represent a stronger correlation while
lighter fields represent a weaker correlation. On top of the heatmap, a dendrogram is visualised that shows the relation
between each replicate to all other replicates.

In case replicates from different conditions indicate a stronger correlation than replicates from the same condition (
possibly due to a batch effect), an additional cleaning step should be executed to remove this effect.

## Distribution plot

This figure displays the distribution of the guides that belong to the class essential genes (positive control) and the
class non-essential genes (negative control). On top of the plot, the f-measure is displayed which represents the
balance between precision (=true positives/(true + false positives)) and recall (=true positives/(true positives + false
negatives)) and is often used to evaluate the performance of a binary classification model.

In optimal conditions, this value should be as close to 1 as possible and the positive control distribution should be
shifted to the left compared to the negative one. In case no separation between positive and negative control can be
seen, this indicates that possibly something went wrong with the screen validation using these positive and negative
controls.

The reason why non-targeting controls are not used as a negative controls is because they do not account for the effect
of cutting itself, like the non-essentials do. Therefore, nowadays it is more common to use non-essential genes as
negative controls for a screen.

## Normal approximation plot

This figure allows to see whether the log<sub>2</sub> fold-change of the conditions that are being compared (stated in
target_samples and reference_samples) can be approximated by a normal distribution. This plot thereby allows to see
whether the assumption of the chosen hit identification model (here DrugZ; assumes read counts follow a normal
distribution ) fits the data. On the x-axis the log<sub>2</sub> fold-change is plotted together with the statistical
significance (p-value) of the observed log<sub>2</sub> fold-change.

In case the normal approximation does not fit the data well, it may be better to use a different model beside DrugZ (
e.g. MAGECK-VISPR) that potentially fits the data better.

## DrugZ: sgRNA output files

These files contain information about the performance of each guide including the mean z-score of the reads for each
guide of all replicates per condition as well as the mean fold-change, the mean standard deviation, the z-score per
guide and the number of guides per gene.

This file is the foundation of the sgRNA plots that are subsequently created and allow to compare the performance of
individual guides in terms of different parameters that were stated above.

## DrugZ: gene output files

These files contain information of each gene based on the guides targeting that gene including the sumZ (sum of z-scores
per gene), normZ (normalised z-score), pval_synth/supp (p-value of synthetic/suppressive interactions),
rank_synth/supp (index of each gene based on the sorted z-score) and fdr_synth/supp (false discovery rate of normZ score
per gene with respect to synthetic/suppressive interaction).

## File containing log<sub>2</sub> fold-change for each comparison

This file represents the concatenated version of the DrugZ gene output files. This file will be used to create the plots
below.

## Gene-level significance plots

These plots show the normZ scores (or log fold-change scores depending on defined value for *x-axis*) of all genes from
the screen in grey on the x-axis while the scores of the significant genes are displayed in colour. On the y-axis, the
false discovery rate of each gene is displayed. In case the current target and reference sample belong to the same
timepoint (e.g. t2_2d) but received a different treatment (e.g. t2_2d_tr vs t2_2d_ut), also a XY-plot is created which
shows the log<sub>2</sub> fold-changes between the initial state (t1_2d) and each of the conditions together with the
hits that were identified in the first plot. In case a *top* was defined in the input parameters, only the top x values
with the lowest false discovery rates will be displayed in the plot. The *fdr_threshold* can also be modified in the
input parameters if needed (default: 0,25).

## Guide-level significance plots

These plots visualize the performance (z-score) of all guides that target one of the genes that was considered
significant in the previous step. On the x-axis, the z-score of each guide is shown together with the -log<sub>10</sub>
of the p-value that corresponds to the z-score in a normal distribution. As a threshold, the -log<sub>10</sub> of the
p-value 0,05 was calculated which corresponds to approximately 1,3. This value was marked in the plots to infer which
sgRNAs have a significant normZ value.

---

# Further information

In case you want to know more about how the analysis is done, please refer to the report "*In silico* and *in vitro*
validation of a genome-wide CRISPR screen for the identification of genes that can confer resistance or sensitivity to
HER2-targeted therapy".

In case you want to know about DrugZ, please refer
to [the Publication describing DrugZ](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0665-3) and
to [the GitHub page](https://github.com/hart-lab/drugz) for the source code.





