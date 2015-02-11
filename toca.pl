#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

##################################################################################
## Initialize Global variables, check for dependecies, parse command line options:
##################################################################################

# Set maximum number of forks to number of free CPUs:
#   - I use forks instead of threads because the computer I run analyses on doesn't have 
#   perl compiled with threads plus I've read that perl's threading has poor performance
my $max_forks = get_free_cpus();

# Fork PIDs 
my @pids;  

# Files which should be deleted upon completion or termination
my @unlink;

# Minimum length (in nucleotides) a contig must be to be used in the analysis
my $min_contig_length = 300;

# Minimum algebraic-connectivity for an orthologous group from ProteinOrtho.pl to be used
my $alg_conn_threshold = 0.25;

# MrBayes settings:
#  - These will result in 250,000 samples and a 100,000 generation burnin
my $nruns = 4;
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
my $blastn = check_path_for_exec("blastn"); # required by proteinortho
my $makeblastdb = check_path_for_exec("makeblastdb"); # required by proteinortho
my $protein_ortho = check_path_for_exec("proteinortho5.pl");

my %polyploids;
my @polyploids;
my @transcriptomes;
GetOptions(
	"transcriptomes|t:s{3,}" => \@transcriptomes,
	"polyploids|p:s{0,}"     => \@polyploids,
	"help|h"                 => \&help,
);
foreach my $polyploid (@polyploids) {
	$polyploids{$polyploid}++;
}

# Input error handling
die "You need to specify at least three fasta files containing transcriptomes.\n".&usage if (scalar(@transcriptomes) < 3);
foreach my $transcriptome (@transcriptomes) {
	die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);
}

####################################
# Run ProteinOrtho and parse output:
####################################

my $project_name = "toca-".time();
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
	#my ($num_species, $num_genes, $alg_conn) = ($line[0], $line[1], $line[2]);
	#next if ($alg_conn < $alg_conn_threshold);

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
					#print "$transcriptomes[$index] has multiple genes but is not labeled as a polyploid, skipping.\n";
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

die "No orthologus families were found.\n" if ($count == 0);
print "$count families passed selection criteria.\n";

# Don't let forks become zombies
$SIG{CHLD} = 'IGNORE';
# Modify SIG_INT (ctrl+c) handler so we can clean up
$SIG{INT} = 'INT_handler';

##############################################
# Analyze up to $max_forks families at a time:
##############################################

foreach my $family (sort { $a <=> $b } keys %families) {
	my $members = $families{$family};
	my @quartets = reduce_family_to_quartets($members);

	foreach my $index (0 .. $#quartets) {
		my $quartet = $quartets[$index];

		# Wait until a CPU is available
		until(okay_to_run()) {};

		my $pid = fork();
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
foreach my $pid (@pids) {
	waitpid($pid, 0);
}

my $out_name;
foreach my $transcriptome (@transcriptomes) {
	(my $base_name = $transcriptome) =~ s/(.*\/)?(.*)\..*/$2/;
	$out_name .= $base_name."-";
}
$out_name .= "a_$alpha";

if ($alpha !~ /^inf/i) {
	system("$bucky -a $alpha -n $ngen_bucky -o $out_name *.sum");
}
else {
	system("$bucky --use-independence-prior -n $ngen_bucky -o $out_name *.sum");
}

sub reduce_family_to_quartets {
	my $members = shift;

	my @members = @{$members};

	my %member_map;
	my @all_members;
	my %nonvariable_members;
#	#foreach my $member (@members) {
	foreach my $index (0 .. $#members) {
		my $member = $members[$index];
		#push(@all_members, @{$member});

		foreach my $index2 (0 .. scalar(@{$member}) - 1) {
			my $member_id = "$index-$index2";
			$member_map{$member_id} = @{$member}[$index2];
#			print "$member_id is @{$member}[$index2]\n";
			push(@all_members, $member_id);
		}

#		if (scalar(@{$member} == 1)) {
##			my $member_name = @{$member}[0];
#			print "$index-0 is non variable.\n";
##			$nonvariable_members{$member_name}++;
#			$nonvariable_members{"$index-0"}++;
#		}
#		else {
##			my $members = @{$members};
##			foreach my $index1 (0 .. $#members - 1) {
##				my $member1 = $members[$index1];
##				foreach my $index2 (1 .. $#members) {
##					my $member2 = $members[$index2];
##					$
##				}
##			}
##			print "variable: @{$member}\n";
#		}
	}
#	print "ALL: @all_members\n";
	
	#my @possible = combine(\@all_members, 4);
	my @possible = combine(\@all_members, scalar(@transcriptomes));

	#print "@$_\n" for @possible;

	my @return;
	foreach my $quartet (@possible) {
		my @members = @{$quartet};
		
		my %species_count;
		my $nonvar_count = 0;
		foreach my $member (@members) {
#			if (exists($nonvariable_members{$member})) {
#				$nonvar_count++;
#				#$nonvar_count += $nonvariable_members{$member};
#			}
			(my $species = $member) =~ s/-\d+//;
			$species_count{$species}++;
		}

#		print "$nonvar_count (".scalar(keys %nonvariable_members).") ".scalar(keys %species_count)." (".scalar(@transcriptomes).")\n";
#		print "has: ".scalar(keys %species_count)." needs: ".scalar(@transcriptomes)."\n";

		#if ($count == scalar(keys %nonvariable_members)) {
		#if ($nonvar_count == scalar(keys %nonvariable_members) && scalar(keys %species_count) == scalar(@transcriptomes)) {
		if (scalar(keys %species_count) == scalar(@transcriptomes)) {
#			print "@{$quartet}\n";

			my @quartet;
			foreach my $contig (@{$quartet}) {
				push(@quartet, $member_map{$contig});
			}
#			print "@quartet\n";
			push(@return, \@quartet);
		}
	}
	return @return;
}

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
	#my $members = @{$args->{MEMBERS}};

	my %sequences;
	foreach my $index (0 .. $#members) {
		(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\..*/$2/;
		my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $members[$index]});
		$sequences{$base_name} = $sequence;
		#print ">$base_name\n$sequence\n\n";
	}

#	my %sequences;
#	#foreach my $index (0 .. scalar(@{$members}) - 1) {
#	foreach my $index (0 .. $#members) {
#		#my $member = @{$members}[$index];
#		my $member = $members[$index];
#		(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\..*/$2/;
#		foreach my $contig (@{$member}) {
#			my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $contig});
#			$sequences{$base_name} = $sequence;
#			print ">$base_name\n$sequence\n\n";
#		}
#	}
#	die;

#	foreach my $index (0 .. $#members) {
#		(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\..*/$2/;
#		my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $members[$index]});
#		$sequences{$base_name} = $sequence;
#	}

	my $family_raw = "$id-raw.fasta";
	my $family_aligned = "$id-aligned.nex";
	my $family_reduced = "$id-reduced.fasta";

	push(@unlink, $family_raw, $family_reduced, $family_raw.".nhr", $family_raw.".nin", $family_raw.".nsq", $family_raw.".out", $family_aligned.".ckp", $family_aligned.".ckp~", $family_aligned.".mcmc");
	foreach my $run (1 .. $nruns) {
		push(@unlink, "$family_aligned.run$run.t");
		push(@unlink, "$family_aligned.run$run.p");
	}
	$SIG{INT} = sub { unlink(@unlink); exit(0) };
	$SIG{KILL} = sub { unlink(@unlink); exit(0) };

	open(my $raw, ">", $family_raw);
	foreach my $taxon (sort { $a cmp $b } keys %sequences) {
		print {$raw} ">$taxon\n";
		print {$raw} "$sequences{$taxon}\n";
	}
	close($raw);

	#print "[$id] Reducing sequence to homologous sites.\n";

	reorient_contigs($id, \%sequences);
	%sequences = parse_fasta($family_raw);

	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

	my $has_minus;

	# remove matches in reverse direction
	foreach my $index (reverse(0 .. $#lines)) {
		my @line = split("\t", $lines[$index]);
		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		if ($s_strand eq "minus") {
			splice(@lines, $index, 1);	
			$has_minus++;
		}
	}
	#exit(0) if (!defined($has_minus));

	my %counts;
	#my %strands;
	my %matches;
	my %quartet;
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];

		my @line = split("\t", $line);
		#my ($query, $match, $q_start, $q_end) = ($line[0], $line[1], $line[2], $line[3]);
		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		next if ($query eq $match);

#		next if ($s_strand eq "minus"); ###################

		$counts{$query}++;

	#	if ($s_strand eq "minus" && !exists($strands{$query})) {
	#		$strands{$match}++;
	#	}

		my $next_line = " \t \t";
		if ($index + 1 <= $#lines) {
			$next_line = $lines[$index + 1];
		}
		my @next_line = split("\t", $next_line);
		my ($next_query, $next_match) = ($next_line[0], $next_line[1]);

		$q_start--; $q_end--; # blast uses 1-based indexing and we don't want that

		my $q_range = "$q_start-$q_end";

		if (exists($matches{"$query-$match"})) {
			$q_range = union($q_range, $matches{"$query-$match"});
		}

	#	print "$query-$match: $q_range ($s_strand)\n";
	#	print "  next : $next_query $next_match\n";
		if ("$query-$match" ne "$next_query-$next_match") {
			if (exists($quartet{$query})) {
				$quartet{$query} = intersection($q_range, $quartet{$query});	
			}
			else {
				$quartet{$query} = $q_range;
			}
		}
		#print "\n" if $next_query ne $query;
		$matches{"$query-$match"} = $q_range;
	}
	if (keys %quartet < 4) {
		unlink(@unlink);
		die "[$id] Could not identify homologous sites for a species in this family.\n";
	}
	foreach my $taxon (keys %counts) {
		my $count = $counts{$taxon};
		#if ($count != scalar(@transcriptomes) - 1) {
		if ($count < scalar(@transcriptomes) - 1) {
			#die "[$id] Some taxa did not share homology.\n";
			unlink(@unlink);
			exit(0);
		}
	}

#	foreach my $contig (keys %strands) {
#		print "$contig: $strands{$contig}\n";
#	}

	open(my $out, ">", $family_reduced);
	#foreach my $contig (keys %quartet) {
	foreach my $contig (sort { $a cmp $b } keys %quartet) {
		#print "frame: $frames{$contig}\n";

#		my $seq;
#		if ($strands{$contig}) {
#			$seq = rev_comp($sequences{$contig});
#		}
#		else {
#			$seq = $sequences{$contig};
#		}
		my $seq = $sequences{$contig};
		my $range = $quartet{$contig};

		print {$out} ">$contig\n";
#		print "[$id] $contig : $range\n";

		my $reduced_length = 0;
		foreach my $segment (split(",", $range)) {
			#print {$out} get_partial_seq($segment, $seq),"\n";
			my $sequence = get_partial_seq($segment, $seq);
			print {$out} "$sequence\n";

			$reduced_length += length($sequence);
		}
		if ($reduced_length < $min_contig_length) {
			unlink(@unlink);
			die "[$id] A sequence did not meet the minimum length requirements.\n";
		}
		print {$out} "\n";
	}
	close($out);

	print "[$id] Aligning reduced sites with muscle.\n";
	system("$muscle -in $family_reduced -out $family_aligned >/dev/null 2>&1") || die;
	#system("$clustalw2 $family_reduced -output=nexus -outfile=$family_aligned >/dev/null 2>&1") || die;

	my %family_aligned = parse_fasta($family_aligned);
	write_nexus({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned});
	#die;

	print "[$id] Running MrBayes and summarizing runs.\n";

	#my $trim = $ngen * $burnin * $burnin + $nruns;
	my $trim = ($ngen * $burnin * $burnin + $nruns) / $nruns;
	system("$mb $family_aligned >/dev/null") || die;
	#system("$mbsum $family_aligned.*.t -n $trim -o $id.sum >/dev/null") || die;
	system("$mbsum $family_aligned.*.t -n $trim -o $id.sum >/dev/null 2>&1") || die;

	unlink(@unlink);
}

sub INT_handler {
	unlink(@unlink); 
	kill -9, @pids; 
	exit(0);
}

sub reorient_contigs { 
	#my ($family_raw, $sequences) = (@_);
	my ($id, $sequences) = (@_);

	my $family_raw = "$id-raw.fasta";
	#my $family_aligned = "$id-aligned.nex";
	#my $family_reduced = "$id-reduced.fasta";

	my %sequences = %{$sequences};
	#(my $id = $family_raw) =~ s/^(\d+)-.*/$1/;

	system("$makeblastdb -in $family_raw -input_type fasta -dbtype nucl >/dev/null") || die;
	#system("$blastn -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstrand' -evalue 1e-05 -perc_identity=50 >/dev/null") || die;
	system("$blastn -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstrand' -evalue 1e-05 >/dev/null") || die;

	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

	#my %matches;
	#my %counts;
	my %hits;
	my %quartet;
	my %strands;
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];

		my @line = split("\t", $line);
		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		next if ($query eq $match);

		$quartet{$query}++;

		if ($s_strand eq "minus" && !exists($strands{$query}) && !exists($hits{"$query-$match"})) {
		#if ($s_strand eq "minus" && !exists($strands{$query})) {
			$strands{$match}++;
		}
		$hits{"$query-$match"}++;

	#	my $next_line = " \t \t";
	#	if ($index + 1 <= $#lines) {
	#		$next_line = $lines[$index + 1];
	#	}
	#	my @next_line = split("\t", $next_line);
	#	my ($next_query, $next_match) = ($next_line[0], $next_line[1]);

	#	$q_start--; $q_end--; # blast uses 1-based indexing and we don't want that

	#	my $q_range = "$q_start-$q_end";

	#	if (exists($matches{"$query-$match"})) {
	#		$q_range = union($q_range, $matches{"$query-$match"});
	#	}

	#	print "$query-$match: $q_range ($s_strand)\n";
	#	print "$query-$match: ($s_strand)\n";
	#	if ("$query-$match" ne "$next_query-$next_match") {
	#		if (exists($quartet{$query})) {
	#			$quartet{$query} = intersection($q_range, $quartet{$query});	
	#		}
	#		else {
	#			$quartet{$query} = $q_range;
	#		}
	#	}
	#	$matches{"$query-$match"} = $q_range;
	}
#	if (keys %quartet < 4) {
#		#unlink(@unlink);
#		die "[$id] Could not identify homologous sites for a species in this family.\n";
#	}
	foreach my $taxon (keys %quartet) {
		my $count = $quartet{$taxon};
		#if ($count != scalar(@transcriptomes) - 1) {
		if ($count < scalar(@transcriptomes) - 1) {
			unlink(@unlink);
			#die "[$id] Some taxa did not share homology.\n";
			exit(0);
		}
	}

	if (keys %strands > 0) {
#		print "[$id] Reorienting contig(s).\n";
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
			#$align->{$taxon} .= $line;
			$align{$taxon} .= $line;
		}
	}
	return %align;
}

sub write_nexus {
	my $settings = shift;

	my $out_name = $settings->{'OUT'};
	my %align = %{$settings->{'ALIGN'}};

	my $ntaxa = scalar(values %align);
	my $nchar = length((values %align)[0]);

	my $total_gens = $ngen * (1 + $burnin);

	open(my $out, ">", $out_name);
	print {$out} "#NEXUS\nbegin data;\n dimensions ntax=$ntaxa nchar=$nchar;\n ";
	print {$out} "format datatype=dna gap=- missing=?;\n matrix\n";
	foreach my $taxon (sort {$a cmp $b} keys %align) {
		print {$out} "$taxon $align{$taxon}\n";
	}
	print {$out} "\n ;\nend;\n\n";
	print {$out} "begin mrbayes;\n";
	print {$out} "\tset nowarnings=yes;\n";
	print {$out} "\tset autoclose=yes;\n";
	print {$out} "\tlset nst=6 rates=invgamma;\n";
	#print {$out} "\tmcmcp ngen=$ngen burninfrac=$burnin samplefreq=$samplefreq printfreq=10000 diagnfreq=50000 nruns=$nruns nchains=3 temp=0.45 swapfreq=10;\n";
	print {$out} "\tmcmcp ngen=$total_gens burninfrac=$burnin samplefreq=$samplefreq printfreq=10000 diagnfreq=50000 nruns=$nruns nchains=3 temp=0.45 swapfreq=10;\n";
	print {$out} "\tmcmc;\n";
	print {$out} "end;\n";
	close($out);
}

sub okay_to_run {
	my $running_forks = 0;
	foreach  my $index (reverse(0 .. $#pids)) {
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
	
	my $path = $ENV{PATH}.":."; #include current directory as well
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
	#$taxon =~ s/-/_/g;

	my $seq;
	my $is_desired_seq;

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
	
	my $union;

	my @add = split(",", $add);
	my @old = split(",", $old);

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

	my @ends;
	my @starts = ($sites[0]);
	my $previous = $sites[0] - 1;
	#foreach my $site (@sites) {
	foreach my $index (0 .. $#sites) {
		my $site = $sites[$index];
		if ($previous != $site - 1) {
			push(@ends, $sites[$index - 1]);
			push(@starts, $site);
		}
		$previous = $site;
	}
	push(@ends, $sites[$#sites]);

	#print "starts: @starts\n";
	#print "ends:   @ends\n";
	foreach my $index (0 .. $#starts) {
		$union .= $starts[$index]."-".$ends[$index].",";
	}
	chop($union);

	#print "$add + $old = $union\n";

	return $union;
}

sub intersection {
	my ($new, $old) = (@_);

	my @new = split(",", $new);
	my @old = split(",", $old);

	my $intersection;

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

	foreach my $site (keys %old) {
		if ($old{$site} != 2) {
			delete($old{$site});
		}
	}

	my @sites = sort { $a <=> $b } keys %old;

	return "0-0" if (scalar(@sites) == 0);

	my @ends;
	my @starts = ($sites[0]);
	my $previous = $sites[0] - 1;
	#foreach my $site (@sites) {
	foreach my $index (0 .. $#sites) {
		my $site = $sites[$index];
		if ($previous != $site - 1) {
			push(@ends, $sites[$index - 1]);
			push(@starts, $site);
		}
		$previous = $site;
	}
	push(@ends, $sites[$#sites]);

	foreach my $index (0 .. $#starts) {
		$intersection .= $starts[$index]."-".$ends[$index].",";
	}
	chop($intersection);

	#print "$new intersect $old = $intersection\n";

	return $intersection;
}

sub get_partial_seq {
	my ($range, $seq) = (@_);

	my $partial_seq;

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

#	if ($no_forks) {
#		return 1; # assume that at least one cpu is free
#	}
#	else {
#
		# Returns a two-member array containing cpu usage observed by the program top,
		# command is run twice as top's first output is usually inaccurate
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
#	}
}

sub help {
	return "";
}

sub usage {
	return "";
}

