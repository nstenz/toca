# toca.pl: Transcriptome-based Ortholog Concordance Analyzer
This script detects orthologous single-copy protein families present in the input transcriptome files, and then seeks to quantify the observed levels of discordance between all families.

## Dependencies
1. [ProteinOrtho.pl](https://www.bioinf.uni-leipzig.de/Software/proteinortho/)
	* Used to identify orthologous sequences shared across the given transcriptomes
2. [Blast](http://1.usa.gov/1zTP2u6)
	* Used by ProteinOrtho.pl to detect orthologs. Also used to reduce orthologous sequences to only their shared homologous sites
3. [MUSCLE](http://www.drive5.com/muscle/downloads.htm)
	* Used to align orthologous family sequences after reduction to homologous sites
4. [MrBayes](http://mrbayes.sourceforge.net/download.php)
	* Used to generate posterior distributions for each orthologous family
5. [mbsum and BUCKy](http://www.stat.wisc.edu/~ane/bucky/downloads.html)
	* mbsum: Summarizes the output from MrBayes into the required format for BUCKy
	* BUCKy: Quantifies observed discordance between genes

## Script Workflow
Potential orthologous families are first identified using ProteinOrtho. In order to obtain higher quality alignments, the sequences of the members of the detected families are reduced to only shared homologous sites. These reduced sequences are then aligned with MUSCLE. Once the families are aligned, they are then run through MrBayes in order to determine their posterior distributions. The results of mcmc chain are then summarized with mbsum. All resulting MrBayes summaries are then pooled together and run in BUCKy to quantify discordance between the input files.

## Script Usage & Settings
### Usage
At the mininmum, the script requires four FASTA files specfied with -i or --input to run:

```
toca.pl -i transcriptome1.fa transcriptome2.fa transcriptome3.fa transcriptome4.fa
```

### Command Line Options
For further fine-tuning of the script, the following options can also be specified:

| Option Flag(s)             | Option Descripton                                                                                    | Default |
|:---------------------------|:----------------------------------------------------------------------------------------------------:|:-------:|
| -i, --input                |file names of at least four transcriptomes (in FASTA format) to use for analyses                      | none    |
| -p, --polyploids           |file names of transcriptomes which should be treated as polyploids, this allows protein families with multiple copies from the polyploid to run| none |
| -o, --out_dir              |name of the directory to store output files in                                                        | "toca-" + Unix time of script invocation) |
| -l, --min_length           |the minimum sequence length (nucleotides) of each family member in order to be analyzed               | 300 nucleotides |
| -T, --n_threads            |the number of families to analyze concurrently                                                        | current number of free CPUs |
| -c, --alg_conn             |the minimum algebraic connectivity for ProteinOrtho                                                   | 0.25 |
| --mb_nruns                 |the number of runs to be used in the MrBayes mcmc                                                     | 4 |
| --mb_nchains               |the number of chains each run should use in the MrBayes mcmc                                          | 3 |
| --mb_ntemp                 |adjusts the swap rate between chains, lower temperature is less likely to swap                        | 0.45 |
| --mb_burnin                |the proportion of mcmc generations which should be discarded as burnin                                | 0.10 |
| --mb_ngen                  |the number of generations to run the MrBayes mcmc                                                     | 1000000 |
| --mb_sfreq                 |the frequency at which the MrBayes mcmc chain should be samples                                       | 40 |
| --bucky_alpha              |specifies potentially multiple values of alpha to run BUCKy with                                      | 1 |
| --bucky_ngen               |the number of generations to run the BUCKy mcmc                                                       | 1000000 |
| -h, --help                 |display help and exit                                                                                 | N/A |

## Output Files
The following files can be found in the output directory upon successful completion of the script:

* For each alpha specified via the "--bucky-alpha" command line setting, the following files will be output (where N = a particular value of alpha):
	* **BUCKy-alpha_N.concordance**: the most important file output by the script, contains the concordance factors for the quartets relevant to the analysis
	* **BUCKy-alpha_N.cluster**: provides information on clustering of input families
	* **BUCKy-alpha_N.gene**: summary of input family topologies and probabilies
	* **BUCKy-alpha_N.input**: lists file names input to BUCKy
	* **BUCKy-alpha_N.out**: log of STDOUT from BUCKy invocation
* **alignments.tar.gz**: gzipped tarball containing the alignments used in the analysis as well as the MrBayes commands that were used to analyze them
* **mb-mcmc-avgs.txt**: summary of MrBayes mcmc chains across all families used in the analysis, this file can be used primarily to:
	1. Determine if the mcmc chains appeared to have converged by checking the average value of the standard deviation of split frequencies
	2. Determine if the temperature of the mcmc chains should be adjusted by checking the average swap frequencies between chains
* **toca.family-members**: contains the names of the contigs present in each family, not all families listed in this file may be included in the concordance analysis
* **toca.proteinortho**: output by ProteinOrtho, contains final sequence clustering information of input files. This the file parsed by the toca.pl to determine clustering
* **toca.blast-graph**: output by ProteinOrtho, contains similarity scores between contigs as determined by blast
* **toca.proteinortho-graph**: output by ProteinOrtho, contains similarity scores between contigs as determined by blast
