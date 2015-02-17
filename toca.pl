#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Time::HiRes qw(usleep);
use Fcntl qw(:flock SEEK_END);

##################################################################################
## Initialize Global variables, check for dependecies, parse command line options:
##################################################################################

# Set maximum number of forks to number of free CPUs:
#   I use forks instead of threads because the computer I run analyses on doesn't have 
#   perl compiled with threads and I've read that perl's threading has poor performance
my $max_forks = get_free_cpus();

# Fork PIDs 
my @pids;  

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
# TODO: allow user to specify more mb settings during invocation

# BUCKy settings:
my $alpha = 1;
my $ngen_bucky = 1000000;
# TODO: allow user to specify more BUCKy settings during invocation

# Check that dependencies are present in user's PATH
my $mb = check_path_for_exec("mb");
my $mbsum = check_path_for_exec("mbsum");
my $bucky = check_path_for_exec("bucky");
my $muscle = check_path_for_exec("muscle");
my $blastn = check_path_for_exec("blastn"); # required by ProteinOrtho
my $makeblastdb = check_path_for_exec("makeblastdb"); # required by ProteinOrtho
my $protein_ortho = check_path_for_exec("proteinortho5.pl");

my %polyploids;
my @polyploids;
my @transcriptomes;
GetOptions(
	"transcriptomes|t:s{4,}" => \@transcriptomes,
	"polyploids|p:s{0,}"     => \@polyploids,
	"help|h"                 => \&help,
);
foreach my $polyploid (@polyploids) {
	$polyploids{$polyploid}++;
}

# Input error handling
die "You need to specify at least four fasta files containing transcriptomes.\n".&usage if (scalar(@transcriptomes) < 4);
foreach my $transcriptome (@transcriptomes) {
	die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);
}

####################################
# Run ProteinOrtho and parse output:
####################################

my $project_name = "toca-".time();
my $mb_swap_sum = $project_name."-mb-swap.txt";
my $mb_stddev_sum = $project_name."-mb-stddev.txt";

system("$protein_ortho @transcriptomes -p=blastn+ -clean -conn=$alg_conn_threshold -project=$project_name");

# Reduce output only to families containing protein(s) in each transcriptome
my %families;
my $count = 0;
open(my $ortho_output, "<", "$project_name.proteinortho") || die "Could not open ProteinOrtho output: $!.\n";
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
		$count++;

		print "$count: ";
		foreach my $contig (@contigs) {
			print "@{$contig} ";
		}
		print "\n";
	}
}
close($ortho_output);

die "No orthologous families were found.\n" if ($count == 0);
print "$count families passed selection criteria.\n";

# Don't let forks become zombies
$SIG{CHLD} = 'IGNORE';

# Modify SIG_INT (ctrl+c) handler so we can clean up
$SIG{INT} = 'INT_handler';

########################################
# Analyze up to $max_forks concurrently:
########################################

# Run each family
foreach my $family (sort { $a <=> $b } keys %families) {
	my $members = $families{$family};
	my @quartets = reduce_family_to_quartets($members);

	# Run each quartet separately 
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

# Create summaries
summarize_mb_stddev();
summarize_mb_swap_freqs();

########################################
# Run BUCKy on the MrBayes output files:
########################################

$project_name .= "-BUCKy-alpha_$alpha";

if ($alpha !~ /^inf/i) {
	system("$bucky -a $alpha -n $ngen_bucky -o $project_name *.sum");
}
else {
	system("$bucky --use-independence-prior -n $ngen_bucky -o $project_name *.sum");
}

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
	my $family_aligned = "$id-aligned.nex";

	# Only the homologous sequence of the family
	my $family_reduced = "$id-reduced.fasta";

	# Stores STDOUT of MrBayes run
	my $mb_log_name = "$family_aligned.mb.log";

	# Add all of the temporary files we don't care about into the deletion array
	push(@unlink, $family_raw, $family_reduced, $family_raw.".nhr", $family_raw.".nin", $family_raw.".nsq", 
		$family_raw.".out", $family_aligned.".ckp", $family_aligned.".ckp~", $family_aligned.".mcmc", $mb_log_name);

	foreach my $run (1 .. $nruns) {
		push(@unlink, "$family_aligned.run$run.t");
		push(@unlink, "$family_aligned.run$run.p");
	}

	$SIG{INT} = sub { unlink(@unlink); exit(0) };
	$SIG{KILL} = sub { unlink(@unlink); exit(0) };

	# Create the raw sequence file
	open(my $raw, ">", $family_raw);
	foreach my $taxon (sort { $a cmp $b } keys %sequences) {
		print {$raw} ">$taxon\n";
		print {$raw} "$sequences{$taxon}\n";
	}
	close($raw);

	print "[$id] Reducing sequence to homologous sites.\n";

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
	#if (keys %quartet < 4) {
	if (scalar(keys %quartet) < scalar(@transcriptomes)) {
		unlink(@unlink);
		die "[$id] Could not identify any homologous sites for a species in this family.\n";
	}
	foreach my $taxon (keys %counts) {
		my $count = $counts{$taxon};
		if ($count < scalar(@transcriptomes) - 1) {
			unlink(@unlink);
			die "[$id] Some taxa did not share homology.\n";
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
			die "[$id] A sequence did not meet the minimum length requirements.\n";
		}
		print {$out} "\n";
	}
	close($out);

	######################################################
	# Run MUSCLE, MrBayes, then summarize MrBayes results:
	######################################################

	# Align with MUSCLE
	print "[$id] Aligning reduced sites with MUSCLE.\n";
	system("$muscle -in $family_reduced -out $family_aligned >/dev/null 2>&1") || die;

	# Load alignment created by MUSCLE, and rewrite it in NEXUS format with MrBayes commands
	my %family_aligned = parse_fasta($family_aligned);
	write_nexus({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned});

	# Run MrBayes
	print "[$id] Running MrBayes and summarizing runs.\n";
	#system("$mb $family_aligned >/dev/null") || die;
	system("$mb $family_aligned > $mb_log_name") || die;

	# Open MrBayes log, extract swap frequencies + final std dev of split frequencies
	open(my $mb_log, "<", $mb_log_name);
	chomp(my @data = <$mb_log>);
	close($mb_log);

	my @matrix_lines;
	my $final_std_dev;
	foreach my $line (@data) {

		$line =~ s/^\s+|\s+$//g;
		next if ($line eq "");

		if ($line =~ /Average standard deviation of split frequencies: (\S+)/) {
			$final_std_dev = $1;
		}

		if ($line =~ /(\d+) \|/) {
			#print {$storage} $line,"\n";
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
	system("$mbsum $family_aligned.*.t -n $trim -o $id.sum >/dev/null 2>&1") || die;

	unlink(@unlink);
}

sub INT_handler {
	unlink(@unlink); 
	kill -9, @pids; 
	exit(0);
}

sub reorient_contigs { 
	my ($id, $sequences) = (@_);

	my $family_raw = "$id-raw.fasta";

	my %sequences = %{$sequences};

	# Create BLAST database and perform self-BLAST, an evalue of 1e-05 is the default for ProteinOrtho
	system("$makeblastdb -in $family_raw -input_type fasta -dbtype nucl >/dev/null") || die;
	system("$blastn -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstrand' -evalue 1e-05 >/dev/null") || die;

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
			die "[$id] Some taxa did not share homology.\n";
			exit(0);
		}
	}

	# Reverse complement any sequences in the reverse direction then rerun this method
	if (scalar(keys %strands) > 0) {

		print "[$id] Reorienting contig(s).\n";

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
	print {$mb_stddev_sum_file} "Average standard deviation of split frequencies across all families: $summary\n";
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
			my $summary = summarize($swap_values[$y]->[$x]);
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
	print {$mb_swap_sum_file} "Standard deviations are in parentheses\n";

	return;
}

sub parse_fasta {
	my $filename = shift;

	my %align;
	open(my $alignment_file, '<', $filename) 
		or die "Could not open '$filename': $!\n";
	chomp(my @data = <$alignment_file>);
	close($alignment_file);
	
	my $taxon;
	foreach my $line (@data) {
		if ($line =~ /^>(\S+)/) {
		#if ($line =~ /^>(.*)/) {
			$taxon = $1;
		}
		else {
			$taxon =~ s/-/_/g;
			$align{$taxon} .= $line;
		}
	}
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

	my $running_forks = 0;
	foreach my $index (reverse(0 .. $#pids)) {
		my $pid = $pids[$index];
		if (kill 0, $pid) {
			$running_forks++;
		}
		else {
			splice(@pids, $index, 1);
		}
	}

	return ($running_forks < $max_forks);
}

sub check_path_for_exec {
	my $exec = shift;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec);
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

	# Returns a two-member array containing CPU usage observed by top,
	# top is run twice as its first output is usually inaccurate
	chomp(my @percent_free_cpu = `top -bn2d0.05 | grep "Cpu(s)"`);

	my $percent_free_cpu = pop(@percent_free_cpu);
	my $test = $percent_free_cpu;
	$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
	$percent_free_cpu = eval($percent_free_cpu);

	my $total_cpus = `grep 'cpu' /proc/stat | wc -l` - 1;
	die "$test\n" if (!defined($percent_free_cpu));

	my $free_cpus = ceil($total_cpus * $percent_free_cpu / 100);

	if ($free_cpus == 0) {
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
		return sprintf("%.5f (± %.6f)", $mean, $sd);
	}
}

sub help {
	return "";
}

sub usage {
	return "";
}

