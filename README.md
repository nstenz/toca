# toca.pl: Transcriptome-based Ortholog Concordance Analyzer
This script detects orthologous single-copy protein families present in the input transcriptome files, and then seeks to quantify the observed levels of discordance between all families.

## Dependencies
1. [ProteinOrtho.pl](https://www.bioinf.uni-leipzig.de/Software/proteinortho/)
	* Used to identify orthologous sequences shared across the given transcriptomes
2. [Blast](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
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

## Script Usage Settings
### Usage
At the mininmum, the script requires four FASTA files specfied with -i or --input to run:
    toca.pl -i transcriptome1.fa transcriptome2.fa transcriptome3.fa transcriptome4.fa

For further fine tuning of the script, the following options can also be used.

### Command Line Options
  -i, --input                 file names of at least four transcriptomes (in FASTA format) to use for analyses (REQUIRED)
  -p, --polyploids            file names of transcriptomes which should be treated as polyploids, treating a transcriptome 
                              as a polyploid allows protein families with multiple copies in the polyploid to run (default: none)
  -o, --output_directory      name of the directory to store output files in (default: "toca-" + Unix time of script invocation)
  -l, --min_seq_length        the minimum sequence length (nucleotides) of each family member in order to be analyzed (default: 300)
  -T, --num_threads           the number of families to analyze concurrently (default: current number of free CPUs)
  -c, --p_ortho_alg_conn      the minimum algebraic connectivity for ProteinOrtho (default: 0.25)
  --mb_nruns                  the number of runs to be used in the MrBayes mcmc (default: 4)
  --mb_nchains                the number of chains each run should use in the MrBayes mcmc (default: 3)
  --mb_ntemp                  adjusts the swap rate between chains, lower temperature is less likely to swap (default: 0.45)
  --mb_burnin                 the proportion of mcmc generations which should be discarded as burnin (default: 0.10)
  --mb_ngen                   the number of generations to run the MrBayes mcmc (default: 1000000)
  --mb_samplefreq             the frequency at which the MrBayes mcmc chain should be samples (default: 40)
  --bucky_alpha               specifies potentially multiple values of alpha to run BUCKy with (default: 1)
  --bucky_ngen                the number of generations to run the BUCKy mcmc (default: 1000000)
  -h, --help                  display this help and exit

## Output Files

