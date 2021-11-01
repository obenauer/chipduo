#!/usr/bin/perl
# Compare two categories of ChIP-seq BAM files. Requires samtools.
# Usage: ChipDuo file1.bam file2.bam file3.bam ...
# Parameters:
#		-w window size (number of bases)
#		-y window sliding factor (slides per window)
#		-s read size (bases)
#		-b background level (read depth)
#		-T difference threshold
#		-p difference percent threshold
#		-r ratio threshold
#		-m maximum region for samtools
#		-G maximum gene distance from peaks
#		-t number of treatment files
#		-c number of control files
#		-g gene file (filename)
#		-B match only best gene to each peak
#		-e exclude a selected chromosome
#		-O analyze only a selected chromosome
#		-H write out histogram of scores
#		-a avoid normalization
#		-n file of external (spike-in) normalization values
#		-R regions of genome to analyze (BED format)
#		-S suppress header in results file
#		-d input directory of BAMs
#		-D output directory
#		-i input file of previous ChipDuo results
#		-o output file prefix
#		-h help (show command syntax)

# Copyright 2011 John Obenauer.  This program is free software. 
# You can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software 
# Foundation, either version 2 of the License, or (at your option) 
# any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.

use strict;
use Getopt::Long qw(:config no_ignore_case);
use IPC::System::Simple qw(capturex);
my ($i, $j, $k);

# Command syntax for help option
my $command_syntax = "Usage: ChipDuo file1.bam file2.bam file3.bam ...\n" . 
	"    -w window size (number of bases)\n". 
	"    -y window sliding factor (slides per window)\n" . 
	"    -s read size (bases)\n" . 
	"    -b background level (read depth)\n" . 
	"    -T difference threshold\n" . 
	"    -p difference percent threshold\n" . 
	"    -r ratio threshold\n" . 
	"    -m maximum region for samtools\n" . 
	"    -G maximum gene distance from peaks\n" . 
	"    -t number of treatment files\n" . 
	"    -c number of control files\n" . 
	"    -g gene file (filename)\n" . 
	"    -B match only best gene to each peak\n" . 
	"    -e exclude a selected chromosome\n" . 
	"    -O analyze only a selected chromosome\n" . 
	"    -H write out histogram of scores\n" . 
	"    -a avoid normalization\n" . 
	"    -n file of external (spike-in) normalization values\n" . 
	"    -R regions of genome to analyze (BED format)\n" . 
	"    -S suppress header in results file\n" . 
	"    -d input directory of BAMs\n" . 
	"    -D output directory\n" . 
	"    -i input file of previous ChipDuo results\n" . 
	"    -o output file prefix\n" . 
	"    -h help (show command syntax)\n";

# Check that required arguments were passed
if ($#ARGV < 1) {
	die $command_syntax;
}

# Read command-line options
my ($window, $slide, $readlength, $background, $dthresh, $pthresh, $rthresh,
	$maxregion, $maxdist, $tfiles, $cfiles, $genefile, $matchbest, $skipchrom, 
	$onlychrom, $savehist, $avoidnorm, $normfile, $regionfile, $suppresshdr, $inpeakfile, 
	$outfile, $help);
our ($indir, $outdir);
get_options(\@ARGV, $window, $slide, $readlength, $background, $dthresh,
	$pthresh, $rthresh, $maxregion, $maxdist, $tfiles, $cfiles, $genefile, $matchbest, 
	$skipchrom, $onlychrom, $savehist, $avoidnorm, $normfile, $regionfile, 
	$suppresshdr, $indir, $outdir, $inpeakfile, $outfile, $help);
# DEBUG
# print "w = $window\ny = $slide\ns = $readlength\nb = $background\n" .
# 	"T = $dthresh\np = $pthresh\nr = $rthresh\nm = $maxregion\nG = $maxdist\nt = $tfiles\n" .
# 	"c = $cfiles\ng = $genefile\nB = $matchbest\ne = $skipchrom\n" . 
# 	"O = $onlychrom\nH = $savehist\na = $avoidnorm\nn = $normfile\n" . 
#	"R = $regionfile\nS = $suppresshdr\nd = $indir\nD=$outdir\ti = $inpeakfile\no = $outfile\n" . 
#	"h = $help\n";
if ($help) {
	die $command_syntax;
}

# Categorize input files as treatments or controls
my (@files, @categories, $nfiles);
get_files(\@ARGV, \@files, \@categories, $nfiles, $tfiles, $cfiles);
print "Input files:\n";
for ($i = 0; $i <= $#ARGV; $i++) {
    print "$i: $files[$i] ($categories[$i])\n";
}

# Remove ".bam" from input files to make sample names
my @sample_names = ();
my $pos = -1;
for ($i = 0; $i < $nfiles; $i++) {
	$pos = rindex($files[$i], ".");
	if ($pos != -1) {
		$sample_names[$i] = substr($files[$i], 0, $pos);
	} else {
		$sample_names[$i] = $files[$i];
	}
}

# Normalize read counts
my @norm = ();
my $outnormfile = $outfile . "_duo_norm.txt";
if ($avoidnorm) {
	# User doesn't want normalization to be done
	for ($i = 0; $i < $nfiles; $i++) {
		$norm[$i] = 1;
	}
} elsif ($normfile) {
	# Calculate normalization constants using total mapped reads
	get_external_normalization($normfile, \@files, \@norm);
	print "External normalizations:\n";
	for ($i = 0; $i < $nfiles; $i++) {
		print "$files[$i]: $norm[$i]\n";
	}
} else {
	# Calculate normalization constants using total mapped reads
	get_normalization(\@files, \@norm, $outnormfile);
}

# Read genes file if one was specified
my (@genes, $ngenes);
get_genes($genefile, \@genes, $ngenes);
if ($ngenes) {
	print "$ngenes genes read from $genefile.\n";
} else {
	print "No gene file specified.\n";
}
# DEBUG
# print "First 5 genes:\n";
# for ($i = 0; $i < 5; $i++) {
# 	print $genes[$i][3] . " " . $genes[$i][0] . ":" . $genes[$i][1] . "-" . $genes[$i][2] . "\n";
# }

# Define chromosome names and sizes
my ($nchr, @chrnames, @chrsizes);
get_chromosomes($files[0], $nchr, \@chrnames, \@chrsizes);
print "$nchr chromosomes defined.\n";
# DEBUG
#for ($i = 0; $i < $nchr; $i++) {
#	print "$i $chrnames[$i] $chrsizes[$i]\n";
#}

# If option to analyze only a single chromosome was 
# selected, remove other chromosomes from the list
# my ($found, $onlychrom_idx);
# if ($onlychrom ne "") {
# 	$found = 0;
# 	for ($i = 0; $i < $nchr; $i++) {
# 		if ($chrnames[$i] eq $onlychrom) {
# 			$found = 1;
# 			$onlychrom_idx = $i;
# 			last;
# 		}
# 	}
# 	if ($found) {
# 		print "Restricting analysis to chromosome $onlychrom.\n";
# 		$nchr = 1;
# 		$chrnames[0] = $chrnames[$onlychrom_idx];
# 		$chrsizes[0] = $chrsizes[$onlychrom_idx];
# 	} else {
# 		die "Aborting: Selected chromosome $onlychrom not found in genome.\n";
# 	}
# }

my @onlylist = ();
my $nonlylist = 0;
my @old_chrnames = @chrnames;
my @old_chrsizes = @chrsizes;
my $old_nchr = $nchr;
my $found;
if ($onlychrom =~ /\,/) {
	@onlylist = split(/\,/, $onlychrom);
	$nonlylist = @onlylist;
	print "Restricting analysis to chromosomes $onlychrom.\n";
} elsif ($onlychrom ne "") {
	$onlylist[0] = $onlychrom;
	$nonlylist = 1;
	print "Restricting analysis to chromosome $onlychrom.\n";
}
if ($nonlylist > 0) {
	@chrnames = ();
	@chrsizes = ();
	$nchr = 0;
	for ($i = 0; $i < $nonlylist; $i++) {
		$found = 0;
		for ($j = 0; $j < $old_nchr; $j++) {
			if ($onlylist[$i] eq $old_chrnames[$j]) {
				$chrnames[$nchr] = $old_chrnames[$j];
				$chrsizes[$nchr] = $old_chrsizes[$j];
				$nchr++;
				$found = 1;
				last;
			}
		}
		if (!$found) {
			die "Aborting: Selected chromosome $onlylist[$i] not found in genome.\n";
		}
	}
}

# Indicate if a chromosome is being skipped
my @skiplist = ();
my $nskiplist = 0;
if ($skipchrom =~ /\,/) {
	@skiplist = split(/\,/, $skipchrom);
	$nskiplist = @skiplist;
	print "Chromosomes $skipchrom will be ignored.\n";
} elsif ($skipchrom ne "") {
	$skiplist[0] = $skipchrom;
	$nskiplist = 1;
	print "Chromosome $skipchrom will be ignored.\n";
}
# Check whether skipped chromosomes are present in BAM files
for ($i = 0; $i < $nskiplist; $i++) {
	$found = 0;
	for ($j = 0; $j < $nchr; $j++) {
		if ($chrnames[$j] eq $skiplist[$i]) {
			$found = 1;
			last;
		}
	}
	if (!$found) {
		print "Warning: Skipped chromosome $skiplist[$i] not found in BAM files.\n";
	}
}

# Calculate average read length 
my ($first1000, @bamlines, $nbamlines, @parts);
my $readlength_sum = 0;
if ($readlength == -1) {
	# Read in first 1000 reads from first BAM file
	$first1000 = `samtools view $indir$files[0] | head -n 1000`;
	@bamlines = split(/\n/, $first1000);
	$nbamlines = scalar @bamlines;
	for ($i = 0; $i < $nbamlines; $i++) {
		@parts = split(/\t/, $bamlines[$i]);
		$readlength_sum += length($parts[9]);
	}
	if ($nbamlines == 0) {
		die "Aborting: No lines could be read in from first BAM file.\n";
	} else {
		$readlength = int($readlength_sum / $nbamlines);
	}
}
print "Read length = $readlength\n";

# Initialize histogram of all scored windows
# Histogram is 0 to 1 million with bin sizes of 1
my @scorehist = ();
my $nscorehist = 1000000;
for ($i = 0; $i < $nscorehist; $i++) {
	$scorehist[$i][0] = $i; # Every bin size is 1
	$scorehist[$i][1] = 0; # All start out at zero counts
}
my $nscored = 0;

# Read in and compare the data in genome regions
my $margin = $maxregion * 0.01;
if ($margin < 1000) {
	$margin = 1000;
}
my ($reads, $nreads, $m, $p, $r, $w, $command);
my @lines = ();
my ($nwindows, @windows);
my (@counts, @depth, $maxdepth);
my @coverage = ();
my ($start, $index);
my ($ngenelist, @genelist, @sorted);
my ($genelist_string, $distance_string);
my ($dist1, $dist2, $closest, $closestdist);
my $q;
my (@peaks, $npeaks);
my ($genome_chr, $compare_genome_chr, $genes_chr, $compare_genes_chr);
my ($outnorm);
my (@regions, $nregions);
my $total_windows = 0;
my $skipped_windows = 0;
my $analyzed_windows = 0;
my @good_windows = ();
my $maxdepth_region;
my $discard_hdr;
my ($found1, $found2);

# Define array of genome regions that will be analyzed
my @regions = ();
my $nregions = 0;
if ($regionfile) {
	open(REGIONS, $regionfile);
	while(<REGIONS>) {
		chomp;
		@parts = split(/\t/);
		$regions[$nregions][0] = $parts[0]; # Chr 
		$regions[$nregions][1] = $parts[1]; # Start
		$regions[$nregions][2] = $parts[2]; # End
		$regions[$nregions][3] = $parts[2] - $parts[1]; # Size
		$nregions++;
	}
	close(REGIONS);
	print "Restricting analysis to genome regions in $regionfile.\n";
} else {
	for ($i = 0; $i < $nchr; $i++) {

		# Skip if this is an excluded chromosome
		$found = 0;
		for ($j = 0; $j < $nskiplist; $j++) {
			if ($skiplist[$j] eq $chrnames[$i]) {
				$found = 1; 
				last;
			}
		}
		next if $found;

		for ($j = 1; $j <= $chrsizes[$i]; $j = $j + $maxregion) {
			$regions[$nregions][0] = $chrnames[$i];
			$regions[$nregions][1] = $j;
			$regions[$nregions][2] = $j + $maxregion + $margin - 1;
			# Don't let region go past end of chromosome
			if ($regions[$nregions][2] > $chrsizes[$i]) {
				$regions[$nregions][2] = $chrsizes[$i];
			}
			$regions[$nregions][3] = $regions[$nregions][2] - $regions[$nregions][1] + 1;
			$nregions++;
		}
	}

	# Write out regions
	open(REGIONS, ">$outdir$outfile" . "_duo_regions.bed");
	for ($i = 0; $i < $nregions; $i++) {
		print REGIONS "$regions[$i][0]\t$regions[$i][1]\t$regions[$i][2]\n";
	}
	close(REGIONS);
	
}

# Examine windows within each region
for ($i = 0; $i < $nregions; $i++) {

	# Define windows in this region
	$nwindows = 0;
	@windows = ();
	for ($w = 0; $w < $regions[$i][3]; $w = $w + $window / $slide) {
		$windows[$nwindows][0] = $w;
		$windows[$nwindows][1] = $w + $window - 1;
		# Make sure window doesn't go past end of this region
		if ($windows[$nwindows][1] > $regions[$i][2]) {
			$windows[$nwindows][1] = $regions[$i][2];
		}
		$windows[$nwindows][2] = $windows[$nwindows][1] 
			- $windows[$nwindows][0] + 1;
		$nwindows++;
	}
	$total_windows += $nwindows;
	@good_windows = ();

	print "Scanning region $regions[$i][0]:$regions[$i][1]-$regions[$i][2] ($nwindows windows)...\n";
	
	# Fill coverage array
	&get_coverage($regions[$i][0], $regions[$i][1], $regions[$i][2], $regions[$i][3], 
		$readlength, $nfiles, \@files, \@coverage);
	$maxdepth_region = 0;

	# Loop through windows to find treatment vs. control differences
	for ($w = 0; $w < $nwindows; $w++) {

		# Create depth array from coverage, using windows
		@depth = ();
		$maxdepth = 0;
		for ($k = 0; $k < $nfiles; $k++) {
			$depth[$k] = 0;	# Initialize
			for ($r = $windows[$w][0]; $r <= $windows[$w][1]; $r++) {
				#print "  coverage[$k][$r] = $coverage[$k][$r]\n";
				if ($depth[$k] < $coverage[$k][$r]) {
					$depth[$k] = $coverage[$k][$r];
				}
			}
			# Normalize depth
			$depth[$k] = $depth[$k] * $norm[$k];
			if ($maxdepth < $depth[$k]) {
				$maxdepth = $depth[$k];
			}
		}
		
		if ($maxdepth_region < $maxdepth) {
			$maxdepth_region = $maxdepth;
		}
		
		# Skip this window if the depth in all files is below background
		if ($maxdepth <= $background) {
			$scorehist[0][1]++;
			$nscored++;
			$skipped_windows++;
			$good_windows[$w] = 0;
			next;
		} else {
			$good_windows[$w] = 1;
		}
		
		# Make sure at least one of the categories is above background
		my ($bkgflag, $tbkgflag, $cbkgflag, $tbkgcount, $cbkgcount);
		$bkgflag = $tbkgflag = $cbkgflag = 0;
		$tbkgcount = $cbkgcount = 0;
		# Treatment files
		for ($k = 0; $k < $tfiles; $k++) {
			if ($depth[$k] > $background) {
				$tbkgcount++;
			}
		}
		if ($tbkgcount == $tfiles) {
			$tbkgflag = 1;
		}
		# Control files
		for ($k = $tfiles; $k < ($tfiles + $cfiles); $k++) {
			if ($depth[$k] > $background) {
				$cbkgcount++;
			}
		}
		if ($cbkgcount == $cfiles) {
			$cbkgflag = 1;
		}
		# If either whole category is above background (or both are), set bkgflag
		if ($tbkgcount || $cbkgcount) {
			$bkgflag = 1;
			$good_windows[$w] = 1;
		} else {
			$scorehist[0][1]++;
			$nscored++;
			$skipped_windows++;
			$good_windows[$w] = 0;
			next; # If no full category is above background, skip this window
		}

		# Create counts array from coverage, using windows
		@counts = ();
		for ($k = 0; $k < $nfiles; $k++) {
			$counts[$k] = 0;	# Initialize
			for ($r = $windows[$w][0]; $r <= $windows[$w][1]; $r++) {
				$counts[$k] += $coverage[$k][$r];
			}
			# Count number of reads instead of individual bases
			$counts[$k] = $counts[$k] / $readlength;
			# Normalize counts
			$counts[$k] = $counts[$k] * $norm[$k];
		}

		# Calculate average counts of each category
		my ($tsum, $csum, $tavg, $cavg, $tavgd, $cavgd);
		$tsum = 0;
		for ($k = 0; $k < $tfiles; $k++) {
			$tsum += $counts[$k];
		}
		$tavg = $tsum / $tfiles;
		$csum = 0;
		for ($k = $tfiles; $k < ($tfiles + $cfiles); $k++) {
			$csum += $counts[$k];
		}
		$cavg = $csum / $cfiles;

		# Calculate average depth of each category
		$tsum = 0;
		for ($k = 0; $k < $tfiles; $k++) {
			$tsum += $depth[$k];
		}
		$tavgd = $tsum / $tfiles;
		$csum = 0;
		for ($k = $tfiles; $k < ($tfiles + $cfiles); $k++) {
			$csum += $depth[$k];
		}
		$cavgd = $csum / $cfiles;

		# Check if treatment / control ratio is higher than rthresh
		my ($rcount, $ratio, $ratioflag, $dcount, $difference, $differenceflag, $product, $score, $rwinner, $dwinner, $pwinner, $winner, $log2fc);
		$rwinner = "-";
		$rcount = 0;
		for ($k = 0; $k < $tfiles; $k++) {
			for ($m = 0; $m < $cfiles; $m++) {
				if ($depth[$m + $tfiles] == 0) {	# Protect against division by zero
					if ($depth[$k] > 0) {
						$rcount++;	
						# positive / zero = infinity, which meets any ratio threshold!
					}
				# Compare each treatment to each control
				} elsif (($depth[$k] / $depth[$m + $tfiles]) > $rthresh) {
					$rcount++;
				}
			}
		}
		if ($rcount == ($tfiles * $cfiles)) {
			$ratioflag = 1;
			$rwinner = "Treatment";
			if ($cavgd) {
				$ratio = $tavgd / $cavgd;
			} else {
				$ratio = $tavgd;
			}
			if ($ratio > 0) {
				$log2fc = log($ratio) / log(2);
			} else {
				$log2fc = 0;
			}
		} else {

			# If treatment / control ratio wasn't high enough, try control / treatment ratio
			$rcount = 0;
			for ($k = 0; $k < $tfiles; $k++) {
				for ($m = 0; $m < $cfiles; $m++) {
					if ($depth[$k] == 0) {
						if ($depth[$m + $tfiles] > 0) {
							$rcount++;
						}
					} elsif (($depth[$m + $tfiles] / $depth[$k]) > $rthresh) {
						$rcount++;
					}
				}
			}
			if ($rcount == ($tfiles * $cfiles)) {
				$ratioflag = 1;
				$rwinner = "Control";
				if ($tavgd) {
					$ratio = $cavgd / $tavgd;
				} else {
					$ratio = $cavgd;
				}
				if ($ratio > 0) {
					$log2fc = -1.0 * log($ratio) / log(2);
				} else {
					$log2fc = 0;
				}
			} else {
				# Neither T/C nor C/T ratio worked
				$ratioflag = 0;
				$ratio = 1;
				if ($tavgd > $cavgd) {
					if ($cavgd) {
						$ratio = $tavgd / $cavgd;
					} else {
						$ratio = $tavgd;
					}
					if ($ratio > 0) {
						$log2fc = log($ratio) / log(2);
					} else {
						$log2fc = 0;
					}
				} else {
					if ($tavgd) {
						$ratio = $cavgd / $tavgd;
					} else {
						$ratio = $tavgd;
					}
					if ($ratio > 0) {
						$log2fc = -1.0 * log($ratio) / log(2);
					} else {
						$log2fc = 0;
					}
				}
			}
		}

		# Check if treatment / control difference percent is higher than pthresh
		my ($pcount, $diffpct, $diffpctflag);
		$pwinner = "-";
		$pcount = 0;
		for ($k = 0; $k < $tfiles; $k++) {
			for ($m = 0; $m < $cfiles; $m++) {
				if ($counts[$m + $tfiles] == 0) {	# Protect against division by zero
					if (($counts[$k] - $counts[$m + $tfiles]) > 0) {
						$pcount++;	
						# positive / zero = infinity, which meets any ratio threshold!
					}
				# Compare each treatment to each control
				} elsif ((($counts[$k] - $counts[$m + $tfiles]) / $counts[$m + $tfiles]) > $pthresh) {
					$pcount++;
				}
			}
		}
		if ($pcount == ($tfiles * $cfiles)) {
			$diffpctflag = 1;
			$pwinner = "Treatment";
			if ($cavg) {
				$diffpct = ($tavg - $cavg) / $cavgd;
			} else {
				$diffpct = ($tavg - $cavg);
			}
		} else {

			# If treatment / control diffpct wasn't high enough, try control / treatment diffpct
			$pcount = 0;
			for ($k = 0; $k < $tfiles; $k++) {
				for ($m = 0; $m < $cfiles; $m++) {
					if ($counts[$k] == 0) {
						if (($counts[$m + $tfiles] - $counts[$k]) > 0) {
							$pcount++;
						}
					} elsif ((($counts[$m + $tfiles] - $counts[$k]) / $counts[$k]) > $pthresh) {
						$pcount++;
					}
				}
			}
			if ($pcount == ($tfiles * $cfiles)) {
				$diffpctflag = 1;
				$pwinner = "Control";
				if ($tavg) {
					$diffpct = ($cavg - $tavg) / $tavg;
				} else {
					$diffpct = ($cavg - $tavg);
				}
			} else {
				# Neither T/C nor C/T ratio worked
				$diffpctflag = 0;
				if ($tavg > $cavg) {
					if ($cavg) {
						$diffpct = ($tavg - $cavg) / $cavg;
					} else {
						$diffpct = $tavg - $cavg;
					}
				} else {
					if ($tavg) {
						$diffpct = ($cavg - $tavg) / $tavg;
					} else {
						$diffpct = ($cavg - $tavg);
					}
				}
			}
		}

		# Check if treatment - control difference is higher than dthresh
		$dwinner = "-";
		$dcount = 0;
		for ($k = 0; $k < $tfiles; $k++) {
			for ($m = 0; $m < $cfiles; $m++) {
				if (($counts[$k] - $counts[$m + $tfiles]) > $dthresh) {
					$dcount++;
				}
			}
		}
		if ($dcount == ($tfiles * $cfiles)) {
			$differenceflag = 1;
			$dwinner = "Treatment";
			$difference = $tavg - $cavg;
		} else {

			# If treatment - control difference wasn't high enough, try control - treatment
			$dcount = 0;
			for ($k = 0; $k < $tfiles; $k++) {
				for ($m = 0; $m < $cfiles; $m++) {
					if (($counts[$m + $tfiles] - $counts[$k]) > $dthresh) {
						$dcount++;
					}
				}
			}
			if ($dcount == ($tfiles * $cfiles)) {
				$differenceflag = 1;
				$dwinner = "Control";
				$difference = ($cavg - $tavg);
			} else {
				$differenceflag = 0;
				$difference = abs($tavg - $cavg);
			}
		}

		# Calculate a score for this window
		$score = $ratio * $difference * $diffpct;
		$analyzed_windows++;
		
		# Compare winner samples
		if (($rwinner eq "Treatment") && ($dwinner eq "Treatment") 
			&& ($pwinner eq "Treatment")) {
			$winner = "Treatment";
		} elsif (($rwinner eq "Control") && ($dwinner eq "Control") 
			&& ($pwinner eq "Control")) {
			$winner = "Control";
		} else {
			$winner = "-";
		}
		
		# Check if significance criteria are met
		my $significant;
		if ($bkgflag && $differenceflag && $diffpctflag && $ratioflag) {
			$significant = 1;
		} else {
			$significant = 0;
		}
		
		# Add this score to histogram
		if (int($score) > 999999) {
			$scorehist[999999][1]++;
		} else {
			$scorehist[int($score)][1]++;
		}
		$nscored++;
		
		# Check if background, ratio, and difference thresholds were all met
		if ($significant) {
			# DEBUG
			#printf("Found significant difference: %s:%d-%d, ratio %.2f, 
			#	difference %d, diffpct %.2f, score %.2f\n",
			#	$chrnames[$windows[$w][0]],
			#	($windows[$w][1] + $regionstart),
			#	($windows[$w][2] + $regionstart),
			#	$ratio, $difference, $diffpct, $score);

			# Record this peak
			$peaks[$npeaks][0] = $regions[$i][0];
			$peaks[$npeaks][1] = $windows[$w][0] + $regions[$i][1];
			$peaks[$npeaks][2] = $windows[$w][1] + $regions[$i][1];
			$peaks[$npeaks][3] = $winner;
			$peaks[$npeaks][4] = $log2fc;
			$peaks[$npeaks][5] = $score;
			$peaks[$npeaks][6] = 0; # PValue will be added later
			$peaks[$npeaks][7] = 0; # AdjP will be added later
			$peaks[$npeaks][8] = ""; # Gene will be added later
			$peaks[$npeaks][9] = ""; # Distance will be added later
			$npeaks++;
			
		}
	}
}
# Close output file
close(OUT);

# Write out scores histogram if requested
if ($savehist) {
	open(SCORES, ">$outdir$outfile" . "_duo_hist.txt");
	print SCORES "Bin\tCount\n";
	for ($i = 0; $i < $nscorehist; $i++) {
		print SCORES "$scorehist[$i][0]\t$scorehist[$i][1]\n";
	}
	close(SCORES);
}

# Use scores histogram to calculate p-values
my $scorebin = 0;
my $p_numerator;
for ($i = 0; $i < $npeaks; $i++) {
	$scorebin = int($peaks[$i][5]);
	# Add up counts in this peak's bin and in all higher ones
	$p_numerator = 0;
	for ($j = $scorebin; $j < $nscorehist; $j++) {
		$p_numerator += $scorehist[$j][1];
	}
	# P = p_numerator / nscored
	$peaks[$i][6] = 1.0 * $p_numerator / $nscored;
}

# Sort by winner=T or winner=C and by peak start positions
@peaks = sort { $b->[3] cmp $a->[3] || $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @peaks;
my $ncol = 10;

# Make a second pass through the peak list to join overlapping windows
print "Joining peaks with overlapping regions...\n";
my ($A, $B, $C, $D, $AB, $CD, $AC, $AD, $BC, $BD, $sum);
for ($i = 0; $i < $npeaks; $i++) {
	for ($j = $i + 1; $j < $npeaks; $j++) {
	
		# Make sure conditions are met: don't compare a peak to itself; 
		# chromosome names have to match; HigherGroup (Treatment or Control) 
		# has to match; and neither peak being compared was already deleted
		if (($i != $j) && ($peaks[$i][0] eq $peaks[$j][0]) 
			&& ($peaks[$i][3] eq $peaks[$j][3]) && ($peaks[$i][0] ne "DELETED") 
			&& ($peaks[$j][0] ne "DELETED")) {
			
			# Calculate overlap, if any
			$A = $peaks[$i][1] - 1; # -1/+1 to combine adjacent but non-overlapping peaks
			$B = $peaks[$i][2] + 1;
			$C = $peaks[$j][1] - 1;
			$D = $peaks[$j][2] + 1;
			$AC = $AD = $BC = $BD = 0;
			$AC = 1 if (($A - $C) >= 0);
			$AD = 1 if (($A - $D) >= 0);
			$BC = 1 if (($B - $C) >= 0);
			$BD = 1 if (($B - $D) >= 0);
			$sum = $AC + $AD + $BC + $BD;
			
			if (($sum > 0) && ($sum < 4)) {
				# Overlap detected
				# Use peak with better score, but merge start and end regions
				if ($peaks[$i][5] > $peaks[$j][5]) {
					# Keep first (i) peak
					if ((!$AC) && (!$BD)) {
						# Right side of i peak overlaps left side of j peak; 
						# change i right boundary to equal j right boundary
						$peaks[$i][2] = $peaks[$j][2];
					} elsif (($AC) && (!$BD)) {
						# i peak is contained completely within j peak;
						# change both i boundaries
						$peaks[$i][1] = $peaks[$j][1];
						$peaks[$i][2] = $peaks[$j][2];					
					} elsif ((!$AC) && ($BD)) {
						# j peak is contained completely within i peak;
						# keep i boundaries (do nothing)
					} elsif (($AC) && ($BD)) {
						# Left side of i peak overlaps right side of j peak; 
						# change i left boundary to equal j left boundary
						$peaks[$i][2] = $peaks[$j][2];									
					}
					# Mark j peak for later deletion
					$peaks[$j][0] = "DELETED";
				} else {
					# Keep second (j) peak
					if ((!$AC) && (!$BD)) {
						# Right side of i peak overlaps left side of j peak; 
						# change j left boundary to equal i left boundary
						$peaks[$j][1] = $peaks[$i][1];
					} elsif (($AC) && (!$BD)) {
						# i peak is contained completely within j peak;
						# keep j boundaries (do nothing)
					} elsif ((!$AC) && ($BD)) {
						# j peak is contained completely within i peak;
						# change both j boundaries
						$peaks[$j][1] = $peaks[$i][1];
						$peaks[$j][2] = $peaks[$i][2];
					} elsif (($AC) && ($BD)) {
						# Left side of i peak overlaps right side of j peak; 
						# change j right boundary to equal i right boundary
						$peaks[$j][2] = $peaks[$i][2];
					}
					# Mark i peak for later deletion
					$peaks[$i][0] = "DELETED";
				}
			}
		} 
	}
}

# Remove deleted peaks
my @origpeaks = @peaks;
my $norigpeaks = $npeaks;
@peaks = ();
$npeaks = 0;
for ($i = 0; $i < $norigpeaks; $i++) {
	if ($origpeaks[$i][0] ne "DELETED") {
		for ($j = 0; $j < $ncol; $j++) {
			$peaks[$npeaks][$j] = $origpeaks[$i][$j];
		}
		$npeaks++;
	}
}

# Remove HigherGroup (winner) and score columns from output
@origpeaks = @peaks;
$norigpeaks = $npeaks;
@peaks = ();
$npeaks = 0;
for ($i = 0; $i < $norigpeaks; $i++) {
	$peaks[$npeaks][0] = $origpeaks[$i][0]; # Chr
	$peaks[$npeaks][1] = $origpeaks[$i][1]; # Start
	$peaks[$npeaks][2] = $origpeaks[$i][2]; # End 
	$peaks[$npeaks][3] = $origpeaks[$i][4]; # Log2FC
	$peaks[$npeaks][4] = $origpeaks[$i][6]; # PValue
	$peaks[$npeaks][5] = $origpeaks[$i][7]; # AdjPValue
	$peaks[$npeaks][6] = $origpeaks[$i][8]; # Gene
	$peaks[$npeaks][7] = $origpeaks[$i][9]; # Distance
	$npeaks++;
}
my $outcol = 8;

# If an input file of previous ChipDuo results 
# was specified on the command line, combine the 
# significant results from this run with the 
# input list and remove duplicates
my @input = ();
my $ninput = 0;
my $margin = 10;
my ($nparts, $ab, $cd, $pointsum, $overlap);
if ($inpeakfile) {

	# Read in input peak list
	print "Combining new peaks with input peaks...\n";
	open(IN, "$outdir$inpeakfile");
	$discard_hdr = <IN>;
	while(<IN>) {
		
		# If a header line is present, skip it
		next if (/^Chrom/);
		
		chomp;
		@parts = split(/\t/);
		$nparts = scalar @parts;
		for ($i = 0; $i < $outcol; $i++) {
			$input[$ninput][$i] = $parts[$i];
		}
		$ninput++;
	}
	close(IN);
	
	# Find and remove duplicates
	for ($i = 0; $i < $npeaks; $i++) {

		# Look for overlapping peaks
		$found = 0;
		for ($j = 0; $j < $ninput; $j++) {
			
			# First check if chromosomes match
			if ($peaks[$i][0] eq $input[$j][0]) {
			
				# Next check coordinates
				$A = $peaks[$i][1] - 1;
				$B = $peaks[$i][2] + 1;
				$C = $input[$j][1] - 1;
				$D = $input[$j][2] + 1;
				$AB = abs($A - $B);
				$CD = abs($C - $D);
				if ($A >= $C) { $AC = 1; } else { $AC = 0; }
				if ($A >= $D) { $AD = 1; } else { $AD = 0; }
				if ($B >= $C) { $BC = 1; } else { $BC = 0; }
				if ($B >= $D) { $BD = 1; } else { $BD = 0; }
				$pointsum = $AC + $AD + $BC + $BD;
				
				if (($pointsum == 1) || ($pointsum == 2) || ($pointsum == 3)) {
				
					# Peaks overlap, calculate overlap size
					$overlap = 0;
					if ($pointsum == 1) {
						# Right side of hits peak overlaps 
						# left side of input peak
						$overlap = abs($B - $C);
					} elsif ($pointsum == 2) {
						# One peak is contained inside the 
						# other one; overlap is the size 
						# of the narrower peak
						if ($AB <= $CD) {
							$overlap = $AB;
						} else {
							$overlap = $CD;
						}
					} elsif ($pointsum == 3) {
						# Left side of hits peak overlaps 
						# right side of input peak
						$overlap = abs($D - $A);
					}
				
					# Test if overlap is large enough
					if ($overlap >= $margin) {
					
						# Overlap found, so ignore this peak
						$found = 1; 
						last;
					}
				}
			}
		}
		# If this peak wasn't found, add it to the 
		# input peak list
		if (!$found) {
			for ($j = 0; $j < $outcol; $j++) {
				$input[$ninput][$j] = $peaks[$i][$j];
			}
			$ninput++;
		}
	}
	
	# Copy input peaks plus added peaks to new peaks array
	@peaks = @input;
	$npeaks = $ninput;
}

# Sort combined peaks by increasing p-value
@peaks = sort { $a->[4] <=> $b->[4] } @peaks;

# Use raw p-values to calculate adjusted p-values
my @rawp = ();
my @reverse_idx = ();
for ($i = 0; $i < $npeaks; $i++) {
	$rawp[$i][0] = $i;
	$rawp[$i][1] = $peaks[$i][4];
	$reverse_idx[$i] = $npeaks - $i;
}
my @nipi = ();
for ($i = 0; $i < $npeaks; $i++) {
	$nipi[$i] = $rawp[$reverse_idx[$i] - 1][1] * $npeaks / $reverse_idx[$i]
}
# Calculate cumulative minimum
my @cum_min = ();
$cum_min[0] = $nipi[0];
if ($cum_min[0] > 1) {
	$cum_min[0] = 1;
}
for ($i = 1; $i < $npeaks; $i++) {
	if ($nipi[$i] < $cum_min[$i - 1]) {
		$cum_min[$i] = $nipi[$i];
	} else {
		$cum_min[$i] = $cum_min[$i - 1];
	}
	if ($cum_min[$i] > 1) {
		$cum_min[$i] = 1;
	}
}
# Calculate adjusted p-value and add it to peaks array
for ($i = 0; $i < $npeaks; $i++) {
	$peaks[$i][5] = $cum_min[$reverse_idx[$i] - 1];
}

# Add gene annotations to peaks if a gene file was included
my ($peakpos, $genelist, $distances, $dist, $g);
if ($genefile) {
	for ($i = 0; $i < $npeaks; $i++) {

		$peakpos = int(($peaks[$i][1] + $peaks[$i][2]) / 2);
		@genelist = ();
		$ngenelist = 0;
		$closest = "";
		$closestdist = "NA";
		for ($g = 0; $g < $ngenes; $g++) {
		
			# Check whether gene annotations use "chr" in chromosome names
			if (substr($peaks[$i][0], 0, 3) eq "chr") {
				$genome_chr = 1;
			} else {
				$genome_chr = 0;
			}
			if (substr($genes[$g][0], 0, 3) eq "chr") {
				$genes_chr = 1;
			} else {
				$genes_chr = 0;
			}
			if ($genome_chr) {
				$compare_genome_chr = $peaks[$i][0];
			} else {
				$compare_genome_chr = "chr" . $peaks[$i][0];
			}
			if ($genes_chr) {
				$compare_genes_chr = $genes[$g][0];
			} else {
				$compare_genes_chr = "chr" . $genes[$g][0];
			}

			# Look only at entries with matching chromosome name
			if ($compare_genome_chr eq $compare_genes_chr) {
				$dist1 = $dist2 = $dist = 0;
				$dist1 = $peakpos - $genes[$g][1];
				$dist2 = $peakpos - $genes[$g][2];
				if (abs($dist1) <= abs($dist2)) {
					$dist = $dist1;
				} else {
					$dist = $dist2;
				}
				# If peak is somewhere within the gene, make distance 0
				if ((($dist1 > 0) && ($dist2 < 0)) || (($dist1 < 0) && ($dist2 > 0))) {
					$dist = 0;
				}
				# Keep track of closest gene, even if it's not 
				# within the maxdist distance
				if (($closestdist eq "NA") || (abs($dist) < $closestdist)) {
					$closest = $genes[$g][3];
					$closestdist = $dist;
				}
				
				if (abs($dist) <= $maxdist) {
					# Make sure this gene isn't already in the list
					$found = 0;
					for ($q = 0; $q < $ngenelist; $q++) {
						if ($genelist[$q][0] eq $genes[$g][3]) {
							$found = 1;

							# Gene is already in list; if new distance 
							# is closer, substitute the closer one
							if (abs($dist) < $genelist[$q][1]) {
								$genelist[$q][0] = $genes[$g][3];
								$genelist[$q][1] = $dist;
								# ngenelist is not incremented because it's a substitution
							}
						}
					}

					# Gene wasn't in list already, so add it
					if (!$found) {
						$genelist[$ngenelist][0] = $genes[$g][3];
						$genelist[$ngenelist][1] = $dist;
						$ngenelist++;
					}
				}
			}
		}

		# Collapse genelist array into strings for output
		if ($ngenelist == 0) {
			$genelist_string = "($closest)";
			$distance_string = "($closestdist)";
		} elsif ($ngenelist == 1) {
			$genelist_string = $genelist[0][0];
			$distance_string = $genelist[0][1];
		} else {

			# If 2 or more genes were found, sort by distance
			@sorted = sort { abs($a->[1]) <=> abs($b->[1]) } @genelist;
			
			# Output single closest gene or all nearby genes
			if ($matchbest) {
				$genelist_string = $sorted[0][0];
				$distance_string = $sorted[0][1];
			} else {
				for ($g = 1; $g < $ngenelist; $g++) {
					$genelist_string .= "," . $sorted[$g][0];
					$distance_string .= "," . $sorted[$g][1];
				}
			}
		}
		
		# Add gene and distance information to this peak
		$peaks[$i][6] = $genelist_string;
		$peaks[$i][7] = $distance_string;
		
	}
}

# Write out peaks
open(PEAKS, ">$outdir$outfile" . "_duo_results.txt");
if (!$suppresshdr) {
	print PEAKS "Chrom\tStart\tEnd\tLog2FC\tPValue\t" . 
		"AdjPValue\tGene\tDistance\n";
}
for ($i = 0; $i < $npeaks; $i++) {
	for ($j = 0; $j < $outcol; $j++) {
		print PEAKS $peaks[$i][$j];
		if ($j == ($outcol - 1)) {
			print PEAKS "\n";
		} else {
			print PEAKS "\t";
		}
	}
}
close(PEAKS);

# Initialize counts matrix
my @counts = ();
my $counts_max = $npeaks;
if ($counts_max > 500) {
	$counts_max = 500;
}
for ($i = 0; $i < $counts_max; $i++) {
	for ($j = 0; $j < $nfiles; $j++) {
		$counts[$i][$j] = 0;
	}
}

# Calculate read counts for significant peaks in all samples
print "Calculating read counts for heatmap...\n";
my $regionsize;
for ($i = 0; $i < $counts_max; $i++) {

	# Fill coverage array
	$regionsize = $peaks[$i][2] - $peaks[$i][1] + 1;
	&get_coverage($peaks[$i][0], $peaks[$i][1], $peaks[$i][2], $regionsize, 
		$readlength, $nfiles, \@files, \@coverage);

	# Measure counts from coverage
	for ($j = 0; $j < $nfiles; $j++) {
		for ($k = 0; $k < $regionsize; $k++) {
			$counts[$i][$j] += $coverage[$j][$k];
		}
		# Count number of reads instead of individual bases
		$counts[$i][$j] = $counts[$i][$j] / $readlength;
		# Normalize
		$counts[$i][$j] = $counts[$i][$j] * $norm[$j];
	}
}

# Write out counts for heatmaps
open(COUNTS, ">$outdir$outfile" . "_duo_counts.txt");
print COUNTS "Peak";
for ($i = 0; $i < $nfiles; $i++) {
	print COUNTS "\t" . $sample_names[$i];
}
print COUNTS "\n";
for ($i = 0; $i < $counts_max; $i++) {
	print COUNTS $peaks[$i][0] . ":" . $peaks[$i][1] . "-" . $peaks[$i][2];
	for ($j = 0; $j < $nfiles; $j++) {
		print COUNTS "\t" . $counts[$i][$j];
	}
	print COUNTS "\n";
}
close(COUNTS);

# Write out treatment and control groups for heatmaps
open(SAMPLES, ">$outdir$outfile" . "_duo_samples.txt");
print SAMPLES "Sample\tGroup\n";
for ($i = 0; $i < $nfiles; $i++) {
	print SAMPLES $sample_names[$i] . "\t";
	if ($categories[$i] eq "t") {
		print SAMPLES "Treatment\n";
	} elsif ($categories[$i] eq "c") {
		print SAMPLES "Control\n";
	}
}
close(SAMPLES);

# Write out R script to plot heatmap
&write_heatmap_script($outfile);
print "Total windows defined = $total_windows.\n";
print "Windows analyzed = $analyzed_windows.\n";
print "Windows skipped = $skipped_windows.\n";
print "Finished analysis for window size = $window.\n";

# Read command line options
sub get_options(\@ARGV, $window, $slide, $readlength, $background, $dthresh, $rthresh, $maxregion, $maxdist, $tfiles, $cfiles, $genefile, $matchbest, $skipchrom, $onlychrom, $savehist,  $avoidnorm, $normfile, $suppresshdr, $regionfile, $indir, $outdir, $inpeakfile, $outfile, $help) {

	# Set default values
	$window = 400;
	$slide = 2;
	$readlength = -1; # Will be estimated if not specified on command line
	$background = 8;
	$dthresh = 10;
	$pthresh = 1.2;
	$rthresh = 1.2;
	$maxregion = 10000000;
	$maxdist = 10000;
	$tfiles = 0;
	$cfiles = 0;
	$genefile = "";
	$matchbest = 0;
	$skipchrom = "";
	$onlychrom = "";
	$savehist = 0; 
	$avoidnorm = 0;
	$normfile = "";
	$suppresshdr = 0;
	$regionfile = "";
	$indir = "";
	$outdir = "";
	$inpeakfile = "";
	$outfile = "duo";
	$help = 0;

	# Replace defaults if values were specified
	GetOptions(
		'w=i' => \$window,
		'y=i' => \$slide,
		's=i' => \$readlength,
		'b=i' => \$background,
		'T=i' => \$dthresh,
		'p=f' => \$pthresh,
		'r=f' => \$rthresh,
		'm=i' => \$maxregion,
		'G=i' => \$maxdist, 
		't=i' => \$tfiles,
		'c=i' => \$cfiles,
		'g=s' => \$genefile,
		'B' => \$matchbest, 
		'e=s' => \$skipchrom, 
		'O=s' => \$onlychrom,
		'H' => \$savehist, 
		'a' => \$avoidnorm,
		'n=s' => \$normfile,
		'S' => \$suppresshdr, 
		'R=s' => \$regionfile, 
		'd=s' => \$indir, 
		'D=s' => \$outdir, 
		'i=s' => \$inpeakfile, 
		'o=s' => \$outfile, 
		'h' => \$help
	);
	if (($indir ne "") && (substr($indir, length($indir) - 1, 1) ne "/")) {
		$indir .= "/";
	}
	if (($outdir ne "") && (substr($outdir, length($outdir) - 1, 1) ne "/")) {
		$outdir .= "/";
	}
	
	return;
}

# Categorize treatment and control input files
sub get_files(\@ARGV, \@files, \@categories, $outfile, $nfiles, $tfiles, $cfiles) {

	# First check if command line options gave the file counts
	$nfiles = $#ARGV + 1;
	my $i;
	if (($tfiles) && ($cfiles)) {
		# User specified number of files for each category
		if (($tfiles + $cfiles) == $nfiles) {
			@files = ();
			@categories = ();
			for ($i = 0; $i < $tfiles; $i++) {
				$files[$i] = $ARGV[$i];
				$categories[$i] = "t";
			}
			for ($i = 0; $i < $cfiles; $i++) {
				$files[$tfiles + $i] = $ARGV[$tfiles + $i];
				$categories[$tfiles + $i] = "c";
			}

		} else {
			# File counts don't add up correctly
			print "Aborting:  Treatment and control file counts do not add up to " .
				"the input files specified.\n";
			exit;
		}
	} else {
		# User didn't specify file counts; assume equal numbers of each
		my $replicates = int($nfiles / 2);
		if (($replicates * 2) == $nfiles) {
			print "Assuming $replicates files per T/C category.\n";
			$tfiles = $replicates;
			$cfiles = $replicates;
			@files = ();
			@categories = ();
			for ($i = 0; $i < $tfiles; $i++) {
				$files[$i] = $ARGV[$i];
				$categories[$i] = "t";
			}
			for ($i = 0; $i < $cfiles; $i++) {
				$files[$tfiles + $i] = $ARGV[$tfiles + $i];
				$categories[$tfiles + $i] = "c";
			}

		} else {
			# Unequal numbers of T/C files are present
			print "Aborting: Unequal numbers of treatment and control files are present.\n";
			print "Unequal numbers are allowed if the -t and -c options are used.\n";
			exit;
		}
	}

	return;
}

sub get_normalization(\@files, \@norm, $outnormfile) {

	# Count number of reads in each file for normalization
	my @nreads = ();
	my $minreads = -1;
	my ($i, $j, @stats, @parts);
	for ($i = 0; $i <= $#files; $i++) {
		print "Counting reads in $files[$i]...\n";
		@stats = `samtools flagstat $indir$files[$i] 2>/dev/null`;
		for ($j = 0; $j <= $#stats; $j++) {
			# Allow for both old and new formats of samtools flagstat output
			if ($stats[$j] =~ /(\d+) \+ (\d+) in total/) {
				$nreads[$i] = $1 + $2;
				last;
			} elsif ($stats[$j] =~ /(\d+) in total/) {
				$nreads[$i] = $1;
				last;
			}
		}
		if ($nreads[$i] == 0) {
			print "Warning: No reads found.  This may mean " . 
				"samtools is not installed or the output is not " .
				"being parsed correctly.\n";
		}
		if ($minreads == -1) {
			$minreads = $nreads[$i];
		}
		if ($minreads > $nreads[$i]) {
			$minreads = $nreads[$i];
		}
	}
	for ($i = 0; $i <= $#files; $i++) {
		$norm[$i] = $minreads / $nreads[$i];
	}

	if ($outnormfile) {
		open(OUTNORM, ">$outdir$outnormfile");
		print OUTNORM "File\tNormalization\n";
		for ($i = 0; $i <= $#files; $i++) {
			print OUTNORM $files[$i] . "\t" . $norm[$i] . "\n";
		}
		close(OUTNORM);
	}
	
	return;
}

sub get_external_normalization($normfile, \@files, \@norm) {

	my $normfile = $_[0];
	my ($i, $j);
	my ($found, $index);
	my @normarray;
	my $nnorm = 0;
	open(NORM, "$outdir$normfile");
	my $discard_hdr = <NORM>;
	while(<NORM>) {
		chomp;
		@parts = split(/\t/);
		$normarray[$nnorm][0] = $parts[0];	# Filename
		$normarray[$nnorm][1] = $parts[1];	# Normalization value
		$nnorm++;
	}
	close(NORM);
	
	for ($i = 0; $i <= $#files; $i++) {
		# Find this file's normalization value
		$found = 0;
		for ($j = 0; $j < $nnorm; $j++) {
			if ($files[$i] eq $normarray[$j][0]) {
				$norm[$i] = $normarray[$j][1];
				$found = 1;
				last;
			}
		}
		if (!$found) {
			die "No external normalization values found for $files[$index] in $normfile.\n";		
		}
	}
	
	return;
}

# Read genes file if one was specified
sub get_genes(my $genefile, \@genes, $ngenes) {

	# If a gene file was specified, read in gene positions
	# Format is BED4. If BED formats with more columns are used 
	# (like BED12 or BED14), the additional columns are ignored.
	@genes = ();
	$ngenes = 0;
	my (@parts, $chr, $strand, $start, $end, $gene);
	if ($genefile) {
		open(GENE, $genefile);
		my $discard_hdr = <GENE>;
		while(<GENE>) {
			chomp;
			@parts = split(/\t/);

			$chr = $parts[0];
			$start = $parts[1];
			$end = $parts[2];
			$gene = $parts[3];
			$genes[$ngenes][0] = $chr;
			$genes[$ngenes][1] = $start;
			$genes[$ngenes][2] = $end;
			$genes[$ngenes][3] = $gene;

			$ngenes++;
		}
		close(GENE);
	}

	return;
}

# Define chromosome names and sizes from header of first BAM
sub get_chromosomes($file1, $nchr, \@chrnames, \@chrsizes) {

	my $file1 = $_[0];
	my $command = "samtools view -H $indir$file1";
	my $header = `$command`;
	my @lines = split(/\n/, $header);
	my $nlines = @lines;
	$nchr = 0;
	my $i;
	if ($nlines) {
		for ($i = 0; $i < $nlines; $i++) {
			if ($lines[$i] =~ /^\@SQ\s+SN\:(\S+)\s+LN\:(\d+)\s*$/) {
				$chrnames[$nchr] = $1;
				$chrsizes[$nchr] = $2;
				$nchr++;
			} elsif ($lines[$i] =~ /^\@SQ\s+LN\:(\d+)\s+SN\:(\S+)\s*$/) {
				$chrnames[$nchr] = $2;
				$chrsizes[$nchr] = $1;
				$nchr++;
			}
		}
	} else {
		die "Aborting: Chromosome names and sizes missing from BAM header.\n";
	}
	return;
}

sub get_coverage() {

	my $targetchr = $_[0];
	my $regionstart = $_[1];
	my $regionstop = $_[2];
	my $regionsize = $_[3];
	my $readlength = $_[4];
	my $nfiles = $_[5];
	my $files_ref = $_[6];
	my @files = @$files_ref;
	my $coverage_ref = $_[7];
	my $command;
	my @arguments;
	
	# Initialize coverage array
	my ($i, $j);
	for ($i = 0; $i < $nfiles; $i++) {
		for ($j = 0; $j < $regionsize; $j++) {
			#$coverage[$i][$j] = 0;
			$$coverage_ref[$i][$j] = 0;
		}
	}

	# Fill coverage array for this region
	my ($reads, $nreads, @lines, @parts, $start, $index, $j, $k);
	for ($i = 0; $i < $nfiles; $i++) {
	
		# Old method using backticks; sometimes causes segfaults
		#$command = "samtools view $indir$files[$i] $targetchr:" .
		#	"$regionstart-$regionstop 2>/dev/null";
		#$reads = `$command`;
		
		$command = "samtools";
		@arguments = ("view", "$indir$files[$i]", 
			"$targetchr:$regionstart-$regionstop");
		$reads = capturex($command, @arguments);
		
		@lines = split(/\n/, $reads);
		$nreads = @lines;
		if ($nreads) {
			for ($j = 0; $j < $nreads; $j++) {
				@parts = split(/\t/, $lines[$j]);
				$start = $parts[3];
				for ($k = 0; $k < $readlength; $k++) {
					$index = $start - $regionstart + $k;
					if ($index >= 0) {
						$$coverage_ref[$i][$index]++;
					}
				}
			}
		}
	}
	
	return;
}

sub write_heatmap_script(my $outfile) {
	open(SCRIPT, ">$outfile" . "_duo.R");
	print SCRIPT "library(DESeq2)\n";
	print SCRIPT "library(gplots)\n";
	print SCRIPT "counts_file = \"$outdir$outfile" . "_duo_counts.txt\"\n";
	print SCRIPT "sample_info_file = \"$outdir$outfile" . "_duo_samples.txt\"\n";
	print SCRIPT "comparison_title = \"$outfile ChipDuo\"\n";
	print SCRIPT "heatmap_file = \"$outdir$outfile" . "_duo_heatmap.png\"\n";
	print SCRIPT "count.data = read.csv(counts_file, sep=\"\\t\")\n";
	print SCRIPT "rownames(count.data) = make.names(count.data\$Peak, unique=TRUE)\n";
	print SCRIPT "count.data = count.data[, -1] # Remove Peak column\n";
	print SCRIPT "pheno.data = read.csv(sample_info_file, sep=\"\\t\", row.names=1)\n";
	print SCRIPT "rownames(pheno.data) = make.names(rownames(pheno.data), unique=TRUE)\n";
	print SCRIPT "count.data = round(count.data) # Round to integer\n";
	print SCRIPT "eset = DESeqDataSetFromMatrix(countData = count.data, \n" . 
		"	colData = pheno.data, design = ~ Group)\n";
	print SCRIPT "png(heatmap_file)\n";
	print SCRIPT "par(oma=c(6,0,0,0))\n";
	print SCRIPT "hmcols = colorRampPalette(c(\"blue\", \"white\", \"red\"))(256)\n";
	print SCRIPT "color.map <- function(Group) { if (Group == \"Control\") \n" . 
		"	\"#0000FF\" else \"#FF0000\" }\n";
	print SCRIPT "samplecolors <- unlist(lapply(eset\$Group, color.map))\n";
	print SCRIPT "heatmap.2(counts(eset, normalized=FALSE), \n" . 
		"	main=comparison_title, cexRow = 0.1, cexCol = 1.2, \n" . 
		"	trace=\"none\", scale=\"row\", ColSideColors = samplecolors, \n" . 
		"	key=TRUE, density.info=\"none\", symkey=FALSE, col=hmcols, \n" . 
		"	labRow = rep(\"\", nrow(eset)))\n";
	print SCRIPT "dev.off()\n";
	close(SCRIPT);
	return;
}

