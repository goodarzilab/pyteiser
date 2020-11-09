![Alt text](/images_github/pyteiser_logo_v2_hossein.png?raw=true)
# Pyteiser
A framework for identifying the structural motifs that are informative of whole-genome measurements across all the transcripts

| **Build Status** |
|:---:|
|[![][travis-img]][travis-url]|


### Introduction
pyteiser identifies structural motifs that could explain genome-wide changes in gene expression, stability or other quantitative transcriptomic measures.
pyteiser encodes structural motifs as context-free grammars (CFG) that represents short stem-loop structures along with a known primary sequence. <br>
First, pyteiser generates a comprehensive set of short seeds (of length 8-16) with a given information content: the seeds are generated in the way that they are not too specific but also not too general. An example seed would have a secondary structure of `<<<<......>>>>` and a sequence of `AAUNNGNGNUNAUU`. <br>
Then, the given RNA sequences (for example, 3'UTRs) are scanned for the occurrence of all the seeds. At that stage, each seed is assigned a binary representation vector, showing which transcripts contain matches for this seed and which ones do not. For example, if a user provided expression values for N sequences, the representation vector will be of length N and each element will be set to "1" if the corresponding fragment has a match to the given seed or "0" if it doesn't.  <br>
Then, each seed's binary vector is tested to assess whether it is informative of the input whole-genome measurements. For example, if the expression change values for transcripts [A, B, C and D] are [High, High, Low, Low] and the binary representation vector for seed X is [1, 1, 0, 0] we make a conclusion that the seed X might be informative of the expression changes. To capture such dependency we use Mutual Information (MI) (see Cover and Thomas, 2006). <br>
Then, pyteiser runs several non-parametric statistical tests to determine which seeds are significantly informative of expression changes and which ones are not. The seeds that did not pass the statistical tests get filtered out. <br>
Then, the seeds that passed the tests are classified by "families" - groups of seeds with very similar representation vectors. The initial set of seeds pyteiser works with is redundant, so several seeds can be very similar to each other and match the same sequences. Such redundant groups of seeds get collapsed into groups or "families". <br>
Then, a single representative seed is chosen for family. The representative seeds is then optimized by a greedy algorithm. The algorithm applies small consecutive changes to sequence and structure of the motif up until the changes don't make the seed more informative of genome-wide measurements anymore. For example, if changing the 1st nucleotide of the seed `AAUNNGNGNUNAUU` to 'G' makes the seed's representation vector more informative of the observed expression values, the algorithm keeps such change. The final oprimized seed is being tested for robustness: the statistical tests for seed significance are being repeated with down-sampled expression data. If the seed is not robust to down-sampling of expression data, it gets filtered out. <br>
Finally, the seeds that have passed all the tests are ranked by how informative are they in regard to the transcriptomic measures provided. Text and graphic reports are being printed, showing motif logos, patterns of their representation, their statistical significance etc. <br>

pyteiser is a successor of TEISER (see [link](https://tavazoielab.c2b2.columbia.edu/TEISER/)). pyteiser is built around the same concept. The overall pipeline is similar, however, several changes and improvements were made. pyteiser performs additional testing of seed matches using *in silico* RNA secondary structure prediction. The statistical tests and optimization algorithms were improved. pyteiser is also capable of handling SHAPE RNA secondary structure probing data; such data provide additional information about RNA secondary structures that do or do not exist in the cell *in vivo*.

pyteiser is intentionally designed in a way that allows storage of all the intermediate files in compressed binary format. The reason is that an academic lab is often focused on one model system; if using pyteiser for the same organism / the same transcriptome, there is no need to re-run time-consuming steps of seed occurence profile calculation over and over again; user can save computational resources by just storing intermediate files for a given transcriptome. We provide two different modes of compressing the seed occurence profiles: either store only indices of matching sequences (takes less space if running pyteiser on big number of short sequences); or store full profiles (takes less space if running pyteiser on small number of long sequences).

pyteiser is written in Python 3.7. The computationally heavy funcions are implemented efficiently through extensive usage on [numpy](https://numpy.org/) arrays and [numba](http://numba.pydata.org/) Python compiler.

## Getting started
### Requirements
pyteiser requires the following dependencies to be installed:
- Python modules:
	- numpy
	- numba
	- pandas
- Other tools:
	- ViennaRNA
It is highly recommended to install the dependencies through conda.
install miniconda


### Installation
Currently, pyteiser is distributed as a set of scripts rather than a pip-installable package. Release of pip-installable version is scheduled for February 2020. <br>
You can either run a command `git clone https://github.com/goodarzilab/pyteiser.git` in your terminal or click "Clone or download" in the top right corner of this page.

## Usage

***Small number of seeds, automatic pipeline that could be run on a PC:*** <br>
We provide a single wrapper script that runs the whole pipeline (`pyteiser/pyteiser/wrappers/pyteiser_pipeline.py`) starting with a set of sequences, corresponding measurements (expression or other ones) and a set of seeds, and reports a list of top candidate seeds and their matches in the input sequences. See the specifications of input files and description of parameters below.

***Large number of seeds, a set of scripts that must be run on HPC:*** <br>
Depending on the size of sequence set and on desired number of seeds to be analyzed, the pipeline might require a lot of computing resources; if searching among a large number of seeds (milliones), it might not be possible to run the pipeline on a PC, and a high-performance computing (HPC) machine might be required. Depending on the institution / company you are at the HPC you are using might have different job submission requirements. In particular, the keywords for memory or time requests might differ among individual HPC systems. Also, depending on cluster resources availability, it might be hard to reserve cores on HPC for long enough to be able to run the whole pipeline in a single run. In that case, we recommend running each step of the pipeline individually. We provide frameworks for either running the scripts on your own machine or submitting it to SGE-based HPCs through `qsub` command. For each computationally heavy step of the pipeline, we provide a script that runs on its own and also a script that is adjusted for submission through qsub. Most scripts can be submitted to HPC with a universal script for qsub submission (named `qsub_universal_submission.py`). It lets you (i) specify the keywords the HPC machine you're using requires, (ii) request how much time, memory and cores do you want to use and (iii) submit any of the computationally heavy scripts from the pipeline. A set of wrapper scripts that implement individual steps of the pipeline is provided in `pyteiser/pyteiser/wrappers`. <br>

Input files:
- Necessary:
	- `rna_fastafile`: a fasta file with RNA sequences of interest
	- `exp_values_file`: expression values in a csv format. It should have 2 or more columns, one column with names of of sequences, another column with values. The names of these two columns have to be specified with the `--anno_name_column` and `--measur_column` arguments. All the names of the sequences listed in this file must also have a corresponding sequence record in the fasta file provided with `rna_fastafile`. The measurement values might be integer or float numbers.
	- `seeds_file`: a binary file containing the seeds to search through. Such file can be created with `seed_generator.py` script
- Optional:
	- user can include a file with RNA structure probing data (SHAPE or DMS-seq) to guide the possible match selection. There is no commonly used standard format for SHAPE RNA reactivity data; therefore, we are using the two-column SHAPE file format used by RNAstructure package ([link](https://rna.urmc.rochester.edu/Text/File_Formats.html#SHAPE)). SHAPE file provided by user should contain SHAPE profiles for multiple sequences, separated with `>`, like in fasta file. SHAPE file can be provided to the `filter_profiles_by_folding.py` script with the `--shape_profile` argument

Output files:
The pipeline generates two files:
- `pyteiser_info.bin`: a table, each row corresponds to one motif that has been identified as significant. The columns include: the motif identifier, sequence and structure of the motif, mutual information of the seed occurence profile and the sequences expression data, pvalue for a permutation test, zscore, the result of robustness test and the number of sequences that contain this seed
- `pyteiser_matches.bin`: a table, each row corresponds to one sequence from the user-provided fasta file. The first two columns contain the sequence and its name; every one of the following columns corresponds to one of the motifs that have been identified as significant; in these columns, cells contain sub-sequences matching the given motif

Arguments for the automatic pipeline (parameters for all the individual steps included):
- input / output files:
	- `rna_fastafile`: fasta file with RNA sequences, see above
	- `exp_values_file`: expression values in a csv format, see above
	- `anno_name_column`: column name in exp_values file that contains annotations, see above
	- `measur_column`: column name in exp_values file that contains expression measurements, see above
	- `seeds_file`: file with seeds in binary format, see above
	- `temp_folder`: folder to write temporary files to
	- `out`: output folder
- optional arguments:
	- `nbins`: number of bins for discretization of expression profile
	- `min_occurences`: minimal number of seed occurence in the transcriptome for a seed to be considered at all
	- `maxfreq`: maximal seed frequency in the sequences analyzed for a seed to be considered
	- `n_permutations`: number of permutations for the rank test for a seed
	- `max_pvalue`: p-value threshold for filtering seeds
	- `min_zscore`: z-score threshold for filtering seeds
	- `step_1_jump`: step size at the 1st round of empirical threshold search
	- `step_2_min_interval`: resolution of empirical threshold search
	- `step_1_min_fraction`: minimal fraction of passing seeds for the 1st round of empirical threshold search
	- `step_2_min_fraction`: minimal fraction of passing seeds for the 2nd round of empirical threshold search
	- `step_3_min_fraction`: minimal fraction of passing seeds for the 3rd round of empirical threshold search
	- `min_ratio`: threshold on ratio of CMI/MI for the conditional information test for seed novelty
	- `are_input_seeds_degenerate`: do input seeds contain degenerate nucleotides (can be ignored in most cases)
	- `indices_mode`: choice of compression mode for storing seed occurence profiles: if True, store only indices of matching sequences (takes less space if running pyteiser on big number of short sequences); if False, store full profiles (takes less space if running pyteiser on small number of long sequences)
	- `index_bit_width`: number of bits per one index when compressing seed occurence profiles. Depends on the number of sequences provided. Default value 24, it's enough if number of sequences provided is smaller than 16.777.216
	- `random_noseed`: when choosing the order of positions to optimize, do not set the random number generator to a specific seed
	- `jackknife_n_permutations`: number of permutations for pvalue calculation in jackknife test
	- `jackknife_max_pvalue`: maximal pvalue for jackknife test
	- `jackknife_n_samples`: how many permutations to do in jackknife test
	- `jackknife_fraction_retain`: what fraction of the sample to retain for each test
	- `jackknife_min_fraction_passed`: what fraction of all iterations should
- arguments for `seed_generator.py` script
	- `outfolder`: output folder
	- `prefix`: prefix for naming the seed file
	- `num_motifs_per_file`: maximal number of seeds to write into a single file
	- `min_stem_length`: minimal stem length to consider
	- `max_stem_length`: maximal stem length to consider
	- `min_loop_length`: minimal loop length to consider
	- `max_loop_length`: maximal loop length to consider
	- `print_sequences`: print the sequences of generated seeds
	- `print_structures`: print the structures of generated seed
	- `print_first_motif`: print the first seed after each increase in stem or loop length
	- `min_inf_bases`: minimal number of informative (non-N) bases
	- `max_inf_bases`: maximal number of informative (non-N) bases
	- `minI`: minimal information content
	- `maxI`: maximal information content


You will have to specify the input and output folders you want to use. All the numeric parameters have preset default values; changing them is not recommended unless you have a very specific reason to do so. <br>
Below, the steps of the pipeline are listed along with the name of the corresponding script.

### Usage of the automatic pipeline
We provide example input files; the automatic pipeline can be launched with the command  once the package has been installed. 

### Steps of the pipeline:
#### 1. Generate seeds
	Use pyteiser/seeds_generator.py
#### 2. Convert sequences from fasta to binary format
	Use pyteiser/wrappers/binarize_sequences.py
#### 3. Precalculate seed occurence profiles
	Use pyteiser/wrappers/calculate_seed_profiles.py - run on HPC!
#### 4. Preprocess the expression file of interest
	Use either pyteiser/wrappers/preprocess_expression_profile_ensembl.py or pyteiser/wrappers/preprocess_custom_expression_profile.py
#### 5. Calculate MI values for all the seeds
	Use pyteiser/wrappers/calculate_MI_profiles.py - run on HPC!
#### 5a. (optional) Filter possible seed matches with *in silico* RNA folding algorithm
	Use pyteiser/wrappers/filter_profiles_by_folding.py
#### 6. Choose significance thresholds
	Use pyteiser/wrappers/choose_significant_seeds_v3.py - run on HPC!
#### 7. Combine seeds that passed
	Use pyteiser/wrappers/combine_passed_seeds.py
#### 8. Classify seeds by families
	Use pyteiser/wrappers/filter_passed_seeds_with_CMI.py
#### 9. Optimize seeds
	Use pyteiser/wrappers/optimize_seeds_single_chunk.py - run on HPC! You can submit it with pyteiser/wrappers/qsub_optimize_seeds.py
#### 10. Combine optimized seeds
	Use pyteiser/wrappers/combine_optimized_seeds.py

### License
MIT license

### Citing
See the paper

### About pyteiser
pyteiser has been developed in Goodarzi lab at UCSF by Matvei Khoroshkin and Hani Goodarzi

[travis-img]: https://travis-ci.com/goodarzilab/pyteiser.svg?branch=master
[travis-url]: https://travis-ci.com/github/goodarzilab/pyteiser

