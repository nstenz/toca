#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use Time::HiRes qw(usleep);
use Fcntl qw(:flock SEEK_END);
use POSIX qw(ceil :sys_wait_h);
use File::Path qw(remove_tree);

##################################################################################
## Initialize Global variables, check for dependecies, parse command line options:
##################################################################################

# Set maximum number of forks to number of free CPUs:
#   I use forks instead of threads because the computer I run analyses on doesn't have 
#   perl compiled with threads and I've read that perl's threading has poor performance
my $max_forks = get_free_cpus();

# Fork PIDs 
my @pids;  
my $parent_pid = $$;

# Files which should be deleted upon completion or termination of an analysis
my @unlink;

# Minimum length (in nucleotides) a contig must be to be used in the analysis
my $min_contig_length = 300;

# Minimum algebraic-connectivity for an orthologous group from ProteinOrtho.pl to be used
my $alg_conn_threshold = 0.25;

# MrBayes settings:
#  - These will result in 250,000 samples per family with a 100,000 generation burnin
my $nruns = 4;
my $nchains = 3;
my $temp = 0.45;
my $burnin = 0.10;
my $ngen = 1000000;
my $samplefreq = 40;

# BUCKy settings:
my @alphas;
my $ngen_bucky = 1000000;

# Name of output directory
my $project_name = "toca-".time();

# Check that dependencies are present in user's PATH
my $mb = check_path_for_exec("mb");
my $mbsum = check_path_for_exec("mbsum");
my $bucky = check_path_for_exec("bucky");
my $muscle = check_path_for_exec("muscle");
my $blastn = check_path_for_exec("blastn"); # required by ProteinOrtho
my $makeblastdb = check_path_for_exec("makeblastdb"); # required by ProteinOrtho
my $protein_ortho = check_path_for_exec("proteinortho5.pl");

# Store script settings
my $settings = "@ARGV";

my %polyploids;
my @polyploids;
my @transcriptomes;
GetOptions(
	# General Settings
	"input|i=s{4,}"      => \@transcriptomes,
	"polyploids|p=s{0,}" => \@polyploids,
	"out_dir|o=s"        => \$project_name,
	"min_length|l=i"     => \$min_contig_length,
	"n_threads|T=i"      => \$max_forks,
	# ProteinOrtho
	"alg_conn|c=f"       => \$alg_conn_threshold,
	# MrBayes
	"mb_nruns=i"         => \$nruns,
	"mb_nchains=i"       => \$nchains,
	"mb_temp=f"          => \$temp,
	"mb_burnin=f"        => \$burnin,
	"mb_ngen=i"          => \$ngen,
	"mb_samplefreq=i"    => \$samplefreq,
	# BUCKy
	"bucky_alpha=s{0,}"  => \@alphas,
	"bucky_ngen=i"       => \$ngen_bucky,

	"help|h"             => \&help,
);
foreach my $polyploid (@polyploids) {
	$polyploids{$polyploid}++;
}

# Default to alpha = 1 if unspecified
push(@alphas, 1) if (!@alphas);

# Input error handling
die "You need to specify at least four fasta files containing transcriptomes.\n".&usage if (scalar(@transcriptomes) < 4);
foreach my $transcriptome (@transcriptomes) {
	die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);
}

# Initialize our working directory and create symlinks to transcriptomes
mkdir($project_name) if (!-e $project_name) || die "Could not create directory '$project_name': $!.\n";
foreach my $index (0 .. $#transcriptomes) {
	my $transcriptome = $transcriptomes[$index];
	my $transcriptome_abs_path = abs_path($transcriptome);
	(my $transcriptome_root = $transcriptome) =~ s/.*\///;

	run_cmd("ln -s $transcriptome_abs_path $ENV{PWD}/$project_name/$transcriptome_root");
	$transcriptomes[$index] = $transcriptome_root;
}

chdir($project_name);

####################################
# Run ProteinOrtho and parse output:
####################################

my $family_id = "toca.family-members";
my $final_mb_sum = "mb-mcmc-avgs.txt";
my $mb_swap_sum = "mb-avg-swap-freq.txt";
my $mb_stddev_sum = "mb-avg-split-std-dev.txt";

my $mb_sum_dir = "mb-sums";
my $align_dir = "alignments";

# Echo script invocation
logger('', "\nInvocation: perl toca.pl $settings");

logger('', "Running ProteinOrtho to identify orthologous proteins...\n");
run_cmd("$protein_ortho @transcriptomes -p=blastn+ -clean -conn=$alg_conn_threshold -project=toca");

mkdir($align_dir) if (!-e $align_dir) || die "Could not create directory '$align_dir': $!.\n";
mkdir($mb_sum_dir) if (!-e $mb_sum_dir) || die "Could not create directory '$mb_sum_dir': $!.\n";

# Reduce output only to families containing protein(s) in each transcriptome
logger('', "\nParsing ProteinOrtho output for gene families...");

my %families;
my $count = 0;
open(my $family_id_file, ">", $family_id);
open(my $ortho_output, "<", "toca.proteinortho") || die "Could not open ProteinOrtho output: $!.\n";
FAMILY: while (my $line = <$ortho_output>) {

	chomp($line); 
	my @line = split("\t", $line);

	# Reorder @transcriptomes based on ProteinOrtho output headers
	if ($. == 1) {
		@transcriptomes = @line[3 .. $#line]; 
		next;
	}

	my ($num_species, $num_genes) = ($line[0], $line[1]);

	# Check that family has at least 1 protein from each transcriptome
	if ($num_species == scalar(@transcriptomes) && $num_genes >= scalar(@transcriptomes)) {
		my @contigs = @line[3 .. $#line];	

		my $num_genes = 0;
		foreach my $index (0 .. $#contigs) {
			my $contig = $contigs[$index];

			my @split_contig = split(",", $contig);
			$num_genes += scalar(@split_contig);

			# Check that family only has multiple genes from a single transcriptome if it's labeled as a polyploid
			if (scalar(@split_contig) > 1) {
				if (!exists($polyploids{$transcriptomes[$index]})) {
					next FAMILY;
				}
			}
			$contigs[$index] = \@split_contig;
		}
		$families{$count} = \@contigs;

		print {$family_id_file} "$count:\n";
		foreach my $index (0 .. $#contigs) {
			my $contig = $contigs[$index];
			my $transcriptome = $transcriptomes[$index];
			print {$family_id_file} "  $transcriptome: @{$contig}\n";
		}
		print {$family_id_file} "\n";

		$count++;
	}
}
close($ortho_output);
close($family_id_file);

if ($count == 0) {
	logger('', "No orthologous families were found.\n");
	exit(0);
}
logger('', "$count families passed selection criteria.\n");

# Modify SIGINT (Ctrl+C) handler so we can clean up
$SIG{INT} = 'INT_handler';

########################################
# Analyze up to $max_forks concurrently:
########################################

logger('', "Analyzing each family using a maximum of $max_forks CPU(s)...");

# Run each family
foreach my $family (sort { $a <=> $b } keys %families) {
	my $members = $families{$family};
	my @quartets = reduce_family_to_quartets($members);

	# Run each quartet with a separate fork
	foreach my $index (0 .. $#quartets) {
		my $quartet = $quartets[$index];

		# Wait until a CPU is available
		until(okay_to_run()) {};

		my $pid = fork();

		# The child fork
		if ($pid == 0) {
			setpgrp();
			analyze_family({'ID' => $family."_$index", 'MEMBERS' => $quartet});
			exit(0);
		}
		else {
			push(@pids, $pid);
		}

	}
}

# Wait for all forks to finish
foreach my $pid (@pids) {
	waitpid($pid, 0);
}

logger('', "Analyses completed for all families.\n");
logger('', "Summarizing MrBayes MCMC quality...");

# Create summaries
summarize_mb_stddev();
summarize_mb_swap_freqs();

# Read in std deviation summary
open(my $mb_stddev_sum_file, "<", $mb_stddev_sum);
my @mb_summary_out = <$mb_stddev_sum_file>;
close($mb_stddev_sum_file);

# Add in swap frequency summary
open(my $mb_swap_sum_file, "<", $mb_swap_sum);
push(@mb_summary_out, <$mb_swap_sum_file>);
close($mb_swap_sum_file);

# Output to new file
open(my $final_mb_sum_file, ">", $final_mb_sum);
print {$final_mb_sum_file} @mb_summary_out;
close($final_mb_sum_file);

# Remove files used to create final summary
unlink($mb_swap_sum);
unlink($mb_stddev_sum);

logger('', "MrBayes MCMC quality summaries complete.");

########################################
# Run BUCKy on the MrBayes output files:
########################################

foreach my $alpha (@alphas) {
	logger('', "Running BUCKy with $ngen_bucky MCMC generations and alpha = $alpha...\n");
	if ($alpha !~ /^inf/i) {
		run_cmd("$bucky -a $alpha -n $ngen_bucky -o BUCKy-alpha_$alpha $mb_sum_dir/*.sum");
	}
	else {
		run_cmd("$bucky --use-independence-prior -n $ngen_bucky -o BUCKy-alpha_$alpha $mb_sum_dir/*.sum");
	}
}
logger('', "Cleaning up and organizing remaining files...");

# Clean up MrBayes summary files
remove_tree($mb_sum_dir);

# Tarball and gzip alignments
run_cmd("tar czf alignments.tar.gz $align_dir/");
remove_tree($align_dir);

logger('', "Script execution complete.\n");

sub reduce_family_to_quartets {
	my $members = shift;

	my @members = @{$members};

	# Map contig names to unique IDs

	my %member_map;
	my @all_members;
	my %nonvariable_members;
	foreach my $index (0 .. $#members) {
		my $member = $members[$index];

		foreach my $index2 (0 .. scalar(@{$member}) - 1) {
			my $member_id = "$index-$index2";
			$member_map{$member_id} = @{$member}[$index2];
			push(@all_members, $member_id);
		}
	}

	# TODO: perform combination over only variable members to save memory won't 
	# really matter though unless you use a ton of transcriptomes
	
	# Create all possible combinations of contigs
	my @possible = combine(\@all_members, scalar(@transcriptomes));

	my @return;
	foreach my $quartet (@possible) {
		my @members = @{$quartet};
		
		my %species_count;
		my $nonvar_count = 0;
		foreach my $member (@members) {
			(my $species = $member) =~ s/-\d+//;
			$species_count{$species}++;
		}

		# A valid quartet has all transcriptomes represented
		if (scalar(keys %species_count) == scalar(@transcriptomes)) {
			my @quartet;
			foreach my $contig (@{$quartet}) {
				push(@quartet, $member_map{$contig});
			}
			push(@return, \@quartet);
		}
	}
	return @return;
}

# I grabbed this from StackOverflow so that's why its style is different #DontFixWhatIsntBroken:
# https://stackoverflow.com/questions/10299961/in-perl-how-can-i-generate-all-possible-combinations-of-a-list
sub combine {
	my ($list, $n) = @_;
	die "Insufficient list members" if ($n > @$list);

	return map [$_], @$list if ($n <= 1);

	my @comb;

	for (my $i = 0; $i+$n <= @$list; ++$i) {
		my $val  = $list->[$i];
		my @rest = @$list[$i + 1 .. $#$list];
		push(@comb, [$val, @$_]) for combine(\@rest, $n - 1);
	}

	return @comb;
}

sub analyze_family {
	my $args = shift;

	my $id = $args->{ID};
	my @members = @{$args->{MEMBERS}};

	# Get the sequence of each contig in this family
	my %sequences;
	foreach my $index (0 .. $#members) {
		(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\..*/$2/;
		my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $members[$index]});
		$sequences{$base_name} = $sequence;
	}

	# The unaligned sequences of the family
	my $family_raw = "$id-raw.fasta";

	# The aligned sequences of the family
	my $family_aligned = "$align_dir/$id-aligned.nex";

	# Only the homologous sequence of the family
	my $family_reduced = "$id-reduced.fasta";

	# Stores STDOUT of MrBayes run
	my $mb_log_name = "$id-aligned.nex.mb.log";

	# Add all of the temporary files we don't care about into the deletion array
	push(@unlink, $family_raw, $family_reduced, $family_raw.".nhr", $family_raw.".nin", $family_raw.".nsq", 
		$family_raw.".out", $family_aligned.".ckp", $family_aligned.".ckp~", $family_aligned.".mcmc", $mb_log_name);

	foreach my $run (1 .. $nruns) {
		push(@unlink, "$family_aligned.run$run.t");
		push(@unlink, "$family_aligned.run$run.p");
	}

	$SIG{TERM} = sub { unlink(@unlink); kill -9, $$; exit(0) };

	# Create the raw sequence file
	open(my $raw, ">", $family_raw);
	foreach my $taxon (sort { $a cmp $b } keys %sequences) {
		print {$raw} ">$taxon\n";
		print {$raw} "$sequences{$taxon}\n";
	}
	close($raw);

	logger($id, "Reducing sequence to homologous sites.");

	# Reverse complement contigs (if needed) so they that are all in the same direction
	reorient_contigs($id, \%sequences);

	# Reload the sequences in their new orientation
	%sequences = parse_fasta($family_raw);

	# Open BLAST output created from final reorientation
	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

	# Remove any BLAST hits in the reverse direction as they are unwanted coincidences
	foreach my $index (reverse(0 .. $#lines)) {
		my @line = split("\t", $lines[$index]);

		my $s_strand = $line[4];
		if ($s_strand eq "minus") {
			splice(@lines, $index, 1);	
		}
	}

	#########################################################################
	# Determine indices of homologous sequence for each contig in the family:
	#########################################################################

	my %counts;
	my %matches;
	my %quartet;
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];

		my @line = split("\t", $line);
		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		next if ($query eq $match);

		$counts{$query}++;

		# Check if next line is a hit between the same contigs, if so, we want its homologous sequence 
		# indices to be the union of those sites

		my $next_line = " \t \t";
		if ($index + 1 <= $#lines) {
			$next_line = $lines[$index + 1];
		}
		my @next_line = split("\t", $next_line);
		my ($next_query, $next_match) = ($next_line[0], $next_line[1]);

		# BLAST uses 1-based indexing and we don't want that
		$q_start--; $q_end--; 

		my $q_range = "$q_start-$q_end";

		# Perform the union
		if (exists($matches{"$query-$match"})) {
			$q_range = union($q_range, $matches{"$query-$match"});
		}

		# Set the homologous sequence indices to its intersection with the other contigs
		if ("$query-$match" ne "$next_query-$next_match") {
			if (exists($quartet{$query})) {
				$quartet{$query} = intersection($q_range, $quartet{$query});	
			}
			else {
				$quartet{$query} = $q_range;
			}
		}
		$matches{"$query-$match"} = $q_range;
	}

	# Check that all contigs actually share sequence
	if (scalar(keys %quartet) < scalar(@transcriptomes)) {
		unlink(@unlink);
		logger($id, "Could not identify any homologous sites for a species in this family.");
		exit(0);
	}
	foreach my $taxon (keys %counts) {
		my $count = $counts{$taxon};
		if ($count < scalar(@transcriptomes) - 1) {
			unlink(@unlink);
			logger($id, "Some taxa did not share homology.");
			exit(0);
		}
	}

	# Output the new sequence based on only the homologous sites
	open(my $out, ">", $family_reduced);
	foreach my $contig (sort { $a cmp $b } keys %quartet) {
		my $seq = $sequences{$contig};
		my $range = $quartet{$contig};

		print {$out} ">$contig\n";

		# Get the desired subsequence(s) based on the homologous site indices
		my $reduced_length = 0;
		foreach my $segment (split(",", $range)) {
			my $sequence = get_partial_seq($segment, $seq);
			print {$out} "$sequence\n";

			$reduced_length += length($sequence);
		}
		
		# Check that the sequence meets the minimum length requirements
		if ($reduced_length < $min_contig_length) {
			unlink(@unlink);
			logger($id, "A sequence was too short.");
			exit(0);
		}
		print {$out} "\n";
	}
	close($out);

	######################################################
	# Run MUSCLE, MrBayes, then summarize MrBayes results:
	######################################################

	# Align with MUSCLE
	logger($id, "Aligning reduced sites with MUSCLE.");
	run_cmd("$muscle -in $family_reduced -out $family_aligned >/dev/null 2>&1");

	# Load alignment created by MUSCLE, and rewrite it in NEXUS format with MrBayes commands
	my %family_aligned = parse_fasta($family_aligned);
	write_nexus({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned});

	# Run MrBayes
	logger($id, "Running MrBayes and summarizing runs.");
	run_cmd("$mb $family_aligned > $mb_log_name");

	# Open MrBayes log, extract swap frequencies + final std dev of split frequencies
	open(my $mb_log, "<", $mb_log_name);
	chomp(my @data = <$mb_log>);
	close($mb_log);

	my @matrix_lines;
	my $final_std_dev;
	foreach my $line (@data) {

		$line =~ s/^\s+|\s+$//g;
		next if ($line eq "");

		# Get last standard deviation of split frequencies
		if ($line =~ /Average standard deviation of split frequencies: (\S+)/) {
			$final_std_dev = $1;
		}

		if ($line =~ /(\d+) \|/) {
			push(@matrix_lines, "$line\n");
		}
	}

	# Get file lock on swap sum so multiple processes don't write at the same time
	open(my $mb_swap_sum_file, ">>", $mb_swap_sum);
	flock($mb_swap_sum_file, LOCK_EX) || die "Could not lock '$mb_swap_sum_file': $!.\n";
	seek($mb_swap_sum_file, 0, SEEK_END) || die "Could not seek '$mb_swap_sum_file': $!.\n";

	print {$mb_swap_sum_file} @matrix_lines;

	flock($mb_swap_sum_file, LOCK_UN) || die "Could not unlock '$mb_swap_sum_file': $!.\n";
	close($mb_swap_sum_file);

	# Get file lock on std dev sum so multiple processes don't write at the same time
	open(my $mb_stddev_sum_file, ">>", $mb_stddev_sum);
	flock($mb_stddev_sum_file, LOCK_EX) || die "Could not lock '$mb_stddev_sum_file': $!.\n";
	seek($mb_stddev_sum_file, 0, SEEK_END) || die "Could not seek '$mb_stddev_sum_file': $!.\n";

	# TODO: ensure mcmc chain is diagnosed on final generation 
	print {$mb_stddev_sum_file} "$final_std_dev\n";

	flock($mb_stddev_sum_file, LOCK_UN) || die "Could not unlock '$mb_stddev_sum_file': $!.\n";
	close($mb_stddev_sum_file);

	# The number of generations from each run to exclude as burnin
	my $trim = ($ngen * $burnin * $burnin + $nruns) / $nruns;
	run_cmd("$mbsum $family_aligned.*.t -n $trim -o $mb_sum_dir/$id.sum >/dev/null 2>&1");

	logger($id, "All analyses successfully completed.");

	unlink(@unlink);
}

sub INT_handler {
	logger('', "\rKeyboard interrupt or pipeline error detected, stopping analyses and cleaning up.");

	# Send SIGTERM to forks, they will in turn sent SIGKILL to their pgroup
	kill 15, @pids; 

	sleep(1);

	# Move back into the directory script was called in
	chdir("..");

	# Try to delete directory five times, if it can't be deleted print an error message
	# I've found this method is necessary for analyses performed on AFS drives
	my $count = 0;
	until (!-e $project_name || $count == 5) {
		$count++;

		remove_tree($project_name, {error => \my $err});
		sleep(1);
	}
	logger('', "Could not clean all files in './$project_name/'.") if ($count == 5);

	exit(0);
}

sub reorient_contigs { 
	my ($id, $sequences) = (@_);

	my $family_raw = "$id-raw.fasta";

	my %sequences = %{$sequences};

	# Create BLAST database and perform self-BLAST, an evalue of 1e-05 is the default for ProteinOrtho
	run_cmd("$makeblastdb -in $family_raw -input_type fasta -dbtype nucl >/dev/null");
	run_cmd("$blastn -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstrand' -evalue 1e-05 >/dev/null");

	# Open resulting BLAST output
	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

	# Determine which strands should be reverse complemented

	my (%hits, %quartet, %strands);
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];

		my @line = split("\t", $line);
		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		next if ($query eq $match);

		$quartet{$query}++;

		# Subject strand is in reverse direction, hasn't been labeled as being reversed, and doesn't currently have a hit with query
		if ($s_strand eq "minus" && !exists($strands{$query}) && !exists($hits{"$query-$match"})) {
			$strands{$match}++;
		}
		$hits{"$query-$match"}++;
	}

	# Check that all contigs actually share sequence
	foreach my $taxon (keys %quartet) {
		my $count = $quartet{$taxon};
		if ($count < scalar(@transcriptomes) - 1) {
			unlink(@unlink);
			logger($id, "Some taxa did not share homology.");
			exit(0);
		}
	}

	# Reverse complement any sequences in the reverse direction then rerun this method
	if (scalar(keys %strands) > 0) {

		logger($id, "Reorienting contig(s).");

		open(my $raw, ">", $family_raw);
		foreach my $contig (keys %quartet) {

			my $seq;
			if ($strands{$contig}) {
				$seq = rev_comp($sequences{$contig});
			}
			else {
				$seq = $sequences{$contig};
			}

			print {$raw} ">$contig\n";
			print {$raw} "$seq\n";
		}
		close($raw);
		reorient_contigs($id, \%sequences);
	}
}

sub summarize_mb_stddev {

	# Read in all std devs
	open(my $mb_stddev_sum_file, "<", $mb_stddev_sum);
	chomp(my @data = <$mb_stddev_sum_file>);
	close($mb_stddev_sum_file);

	my $summary = summarize(\@data, "STDDEV");

	# Reopen for output
	open($mb_stddev_sum_file, ">", $mb_stddev_sum);
	print {$mb_stddev_sum_file} "Average standard deviation of split frequencies across all families: $summary\n\n";
	close($mb_stddev_sum_file);

	return;
}

sub summarize_mb_swap_freqs {

	# Read in all swap matrices
	open(my $mb_swap_sum_file, "<", $mb_swap_sum);
	chomp(my @data = <$mb_swap_sum_file>);
	close($mb_swap_sum_file);

	# Parse the individual matrices from the file
	my @swap_matrices;
	foreach my $index (0 .. $#data) {
		my $line = $data[$index];

		if ($line =~ /(\d+) \|/) {
			if ($1 == $nchains) {
				my $matrix = join("\n", @data[$index - $nchains + 1 .. $index]);
				push(@swap_matrices, $matrix);
			}
		}
	}

	# Join individual matrices into a single one
	my @swap_values;
	foreach my $matrix (@swap_matrices) {

		my @matrix = split("\n", $matrix);
		foreach my $y (0 .. $#matrix) {
			my $line = $matrix[$y];
			$line =~ s/\d+ \|\s+//;

			my @line = split(/\s+/, $line);
			splice(@line, $y, 0, "");

			foreach my $x (0 .. $#line) {
				my $number = $line[$x];
				push(@{$swap_values[$y]->[$x]}, $number);
			}
		}
	}

	# Summarize each matrix entry and determine length of longest one for table formatting
	my $longest_entry = 0;
	foreach my $y (0 .. $#swap_values) {
		foreach my $x (0 .. $#{$swap_values[$y]}) {

			my $summary = summarize($swap_values[$y]->[$x], "SWAP");
			if (length($summary) > $longest_entry) {
				$longest_entry = length($summary);
			}
		}
	}

	# Output matrix summary to file
	open($mb_swap_sum_file, ">", $mb_swap_sum);
	print {$mb_swap_sum_file} "Summary of ".scalar(@swap_matrices)." swap matrices:\n\n";

	# The code is ugly but at least the matrix is nice
	my $nchains_width = length($nchains);
	my $field_width = $longest_entry + 2;
	foreach my $y (0 .. $#swap_values) {

		# Add the header before the first row
		my $line;
		if ($y == 0) {
			$line .= (" " x ($nchains_width + 2));
			foreach my $run (1 .. $nchains) {
				$line .= sprintf("%".($field_width - 1)."s", $run);
			}
			$line .= "\n";
			
			my $header_width = length($line);

			$line .= (" " x ($nchains_width + 1)).("-" x ($header_width - ($nchains_width + 2)))."\n";
		}

		# Add chain number to the line
		$line .= sprintf("%${nchains_width}d |", $y + 1);

		# Add actual values from the matrix to the line
		foreach my $x (0 .. $#{$swap_values[$y]}) {
			my $summary = summarize($swap_values[$y]->[$x], "SWAP");
			if ($summary) {
				$line .= sprintf("%${field_width}s", $summary);
			}
			else {
				$line .= ' ' x ($field_width - 1);
			}
		}
		print {$mb_swap_sum_file} "$line\n";
	}

	print {$mb_swap_sum_file} "\nUpper diagonal: Mean proportion of successful state exchanges between chains for all genes\n";
	print {$mb_swap_sum_file} "Lower diagonal: Mean number of attempted state exchanges between chains for all genes\n\n";
	print {$mb_swap_sum_file} "Values in parentheses are standard deviations\n";

	return;
}

sub parse_fasta {
	my $filename = shift;

	my $taxon;
	my %align;
	open(my $alignment_file, '<', $filename) 
		or die "Could not open '$filename': $!\n";

	while (my $line = <$alignment_file>) {
		$line =~ s/^\s+|\s+$//g;

		# Taxon name
		#if ($line =~ /^>(.*)/) {
		if ($line =~ /^>(\S+)/) {
			$taxon = $1;
		}
		else {
			# Taxon sequence
			$taxon =~ s/-/_/g;
			$align{$taxon} .= $line;
		}
	}
	close($alignment_file);
	
	return %align;
}

sub write_nexus {
	my $settings = shift;

	my $out_name = $settings->{'OUT'};
	my %align = %{$settings->{'ALIGN'}};

	# Determine alignment info required in NEXUS format
	my $ntaxa = scalar(values %align);
	my $nchar = length((values %align)[0]);

	# Total number of generations MCMC should be run
	my $total_gens = $ngen * (1 + $burnin);

	open(my $out, ">", $out_name);
	print {$out} "#NEXUS\nbegin data;\n dimensions ntax=$ntaxa nchar=$nchar;\n ";
	print {$out} "format datatype=dna gap=- missing=?;\n matrix\n";
	foreach my $taxon (sort {$a cmp $b} keys %align) {
		print {$out} "  $taxon $align{$taxon}\n";
	}
	print {$out} "\n ;\nend;\n\n";
	print {$out} "begin mrbayes;\n";
	print {$out} "\tset nowarnings=yes;\n";
	print {$out} "\tset autoclose=yes;\n";
	print {$out} "\tlset nst=6 rates=invgamma;\n";
	print {$out} "\tmcmcp ngen=$total_gens burninfrac=$burnin samplefreq=$samplefreq printfreq=10000 diagnfreq=50000 nruns=$nruns nchains=$nchains temp=$temp swapfreq=10;\n";
	print {$out} "\tmcmc;\n";
	print {$out} "end;\n";
	close($out);
}

sub okay_to_run {

	# Free up a CPU by sleeping for 10 ms
	usleep(10000);

	#my $running_forks = 0;
	my $current_forks = scalar(@pids);
	foreach my $index (reverse(0 .. $#pids)) {
	#	my $pid = $pids[$index];
	#	#if (kill 0, $pid) {
	#	if (kill 0, $pid) {
	#		$running_forks++;
	#	}
	#	else {
	#		splice(@pids, $index, 1);
	#	}
		my $pid = $pids[$index];
		my $wait = waitpid($pid, WNOHANG);

		# Successfully reaped child
		if ($wait > 0) {
			$current_forks--;
			splice(@pids, $index, 1);
		}
	}

	#return ($running_forks < $max_forks);
	return ($current_forks < $max_forks);
}

sub check_path_for_exec {
	my $exec = shift;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec && -x $dir.$exec && !-d $dir.$exec);
	}

	die "Could not find the following executable: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path));
	return $exec_path;
}

sub get_seq_from_fasta {
	my $args = shift;

	my $file = $args->{'FILE'};
	my $taxon = $args->{'TAXON'};

	my $seq;
	my $is_desired_seq;

	# Traverse file until desired contig is found

	open(my $align, "<", $file);
	while (my $line = <$align>) {
		chomp($line);

		if ($line =~ /^>\Q$taxon\E/) {
			$is_desired_seq++;
		}
		elsif ($line =~ /^>/ && $is_desired_seq) {
			last;
		}
		elsif ($is_desired_seq) {
			$seq .= $line;
		}
	}
	close($align);

	return $seq;
}

sub union {
	my ($add, $old) = (@_);
	
	my @add = split(",", $add);
	my @old = split(",", $old);

	# Populate a hash with all members of @add and @old

	my %old;
	foreach my $segment (@old) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}
	foreach my $segment (@add) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}

	my @sites = sort { $a <=> $b } keys %old;

	# Convert sites to a string (i.e. 1,2,3,4,6,7,8 => 1-4,6-8)

	my @ends;
	my @starts = ($sites[0]);
	my $previous = $sites[0] - 1;
	foreach my $index (0 .. $#sites) {
		my $site = $sites[$index];
		if ($previous != $site - 1) {
			push(@ends, $sites[$index - 1]);
			push(@starts, $site);
		}
		$previous = $site;
	}
	push(@ends, $sites[$#sites]);

	my $union;
	foreach my $index (0 .. $#starts) {
		$union .= $starts[$index]."-".$ends[$index].",";
	}
	chop($union);

	return $union;
}

sub intersection {
	my ($new, $old) = (@_);

	my @new = split(",", $new);
	my @old = split(",", $old);

	# Populate a hash with all members of @new and @old

	my %old;
	foreach my $segment (@old) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}
	foreach my $segment (@new) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}

	# Remove site if both arrays don't have it
	foreach my $site (keys %old) {
		if ($old{$site} != 2) {
			delete($old{$site});
		}
	}

	my @sites = sort { $a <=> $b } keys %old;

	return "0-0" if (scalar(@sites) == 0);

	# Convert sites to a string (i.e. 1,2,3,4,6,7,8 => 1-4,6-8)

	my @ends;
	my @starts = ($sites[0]);
	my $previous = $sites[0] - 1;
	foreach my $index (0 .. $#sites) {
		my $site = $sites[$index];
		if ($previous != $site - 1) {
			push(@ends, $sites[$index - 1]);
			push(@starts, $site);
		}
		$previous = $site;
	}
	push(@ends, $sites[$#sites]);

	my $intersection;
	foreach my $index (0 .. $#starts) {
		$intersection .= $starts[$index]."-".$ends[$index].",";
	}
	chop($intersection);

	return $intersection;
}

sub get_partial_seq {
	my ($range, $seq) = (@_);

	my $partial_seq;

	#TODO: check if split is necessary (I don't think it is)

	my @range = split(",", $range);
	foreach my $segment (@range) {
		my ($start, $end) = split("-", $segment);
		$partial_seq .= substr($seq, $start, $end - $start + 1);	
	}
	
	return $partial_seq;
}

sub rev_comp {
	my $seq = shift;
	
	my %comp = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C');

	my $rev_comp;
	my $rev = reverse($seq);
	foreach my $index (0 .. length($rev) - 1) {
		my $char = substr($rev, $index, 1);
		$rev_comp .= $comp{$char};
	}

	return $rev_comp;
}

sub get_free_cpus {

	my $os_name = $^O;

	# Returns a two-member array containing CPU usage observed by top,
	# top is run twice as its first output is usually inaccurate
	my @percent_free_cpu;
	if ($os_name eq "darwin") {
		# Mac OS
		chomp(@percent_free_cpu = `top -i 1 -l 2 | grep "CPU usage"`);
	}
	else {
		# Linux
		chomp(@percent_free_cpu = `top -b -n2 -d0.05 | grep "Cpu(s)"`);
	}

	my $percent_free_cpu = pop(@percent_free_cpu);

	if ($os_name eq "darwin") {
		# Mac OS
		$percent_free_cpu =~ s/.*?(\d+\.\d+)%\s+id.*/$1/;
	}
	else {
		# linux 
		$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
		$percent_free_cpu = eval($percent_free_cpu);
	}

	my $total_cpus;
	if ($os_name eq "darwin") {
		# Mac OS
		$total_cpus = `sysctl -n hw.ncpu`;
	}
	else {
		# linux
		$total_cpus = `grep --count 'cpu' /proc/stat` - 1;
	}

	my $free_cpus = ceil($total_cpus * $percent_free_cpu / 100);

	if ($free_cpus == 0 || $free_cpus !~ /^\d+$/) {
		$free_cpus = 1; # assume that at least one cpu can be used
	}
	
	return $free_cpus;
}

sub summarize {
	my ($array, $type) = @_;
	my @array = @{$array};

	# Empty space in matrix
	return '' if ($array[0] eq '');

	# Calculate mean
	my $sum = 0;
	foreach my $num (@array) {
		$sum += $num;	
	}

	my $mean = $sum / scalar(@array);

	# Calculate standard deviation
	my $deviation_square_sum;
	foreach my $num (@array) {
		$deviation_square_sum += ($mean - $num)**2;
	}

	# Divide by zero if sample size is one
	return $mean if (scalar(@array) - 1 == 0);

	my $sd = sqrt($deviation_square_sum / (scalar(@array) - 1));

	# Reformat number of decimal places to show based on input type
	if ($type eq "SWAP") {
		if ($mean > 1) {
			return sprintf("%.0f (± %.1f)", $mean, $sd);
		}
		else {
			return sprintf("%.2f (± %.2f)", $mean, $sd);
		}
	}
	else {
		return sprintf("%.6f (± %.6f)", $mean, $sd);
	}
}

sub logger {
	my ($id, $msg) = @_;

	my $time = "[".localtime(time())."]";

	# Allow for new lines and carriage returns before message
	if ($msg =~ s/^\n//) {
		$time = "\n$time";
	}
	elsif ($msg =~ s/^\r//) {
		$time = "\r$time";
	}

	my $longest_possible_id = length(scalar(keys %families)) + 2;

	if (length($id) > 0) {
		$id = "0".$id while (length($id) < $longest_possible_id);

		printf("$time [%${longest_possible_id}s] $msg\n", $id); 
	}
	else {
		print "$time $msg\n"; 
	}
}

sub run_cmd {
	my $command = shift;

	my $return = system($command);

	if ($return) {
		logger('', "'$command' died with error: '$return'.\n");
		kill(2, $parent_pid);
		exit(0);
	}
}

sub help {
print <<EOF; 
@{[usage()]}
Identify orthologous protein families shared between the given transcriptomes

  -i, --input                 file names of at least four transcriptomes (in FASTA format) to use for analyses (REQUIRED)
  -p, --polyploids            file names of transcriptomes which should be treated as polyploids, treating a transcriptome 
                              as a polyploid allows protein families with multiple copies in the polyploid to run (default: none)
  -o, --out_dir               name of the directory to store output files in (default: "toca-" + Unix time of script invocation)
  -l, --min_length            the minimum sequence length (nucleotides) of each family member in order to be analyzed (default: 300)
  -T, --n_threads             the number of families to analyze concurrently (default: current number of free CPUs)
  -c, --alg_conn              the minimum algebraic connectivity for ProteinOrtho (default: 0.25)
  --mb_nruns                  the number of runs to be used in the MrBayes mcmc (default: 4)
  --mb_nchains                the number of chains each run should use in the MrBayes mcmc (default: 3)
  --mb_ntemp                  adjusts the swap rate between chains, lower temperature is less likely to swap (default: 0.45)
  --mb_burnin                 the proportion of mcmc generations which should be discarded as burnin (default: 0.10)
  --mb_ngen                   the number of generations to run the MrBayes mcmc (default: 1000000)
  --mb_samplefreq             the frequency at which the MrBayes mcmc chain should be samples (default: 40)
  --bucky_alpha               specifies potentially multiple values of alpha to run BUCKy with (default: 1)
  --bucky_ngen                the number of generations to run the BUCKy mcmc (default: 1000000)
  -h, --help                  display this help and exit

Examples:
  perl toca.pl -i t1.fa t2.fa t3.fa t4.fa     runs the script with default settings using the t1.fa, t2.fa, t3.fa, and t4.fa as input

Mail bug reports and suggestions to <noah.stenz.github\@gmail.com>
EOF
exit(0);
}

sub usage {
	return "Usage: perl toca.pl -i FILES [OPTIONS]...\n";
}
