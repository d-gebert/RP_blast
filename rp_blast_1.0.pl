#!/usr/bin/perl
use strict;
use warnings;

# Options
my $max_eval = 3;		# Maximal e-value
my $min_qcov = 50;		# Query coverage threshold
my $min_hit_ident = 0; 	# Minimal hit identity
my $min_grp_ident = 0;	# Minimal group hit identity
# Global constants
my $tbnopt = '-seg no -word_size 2';
# Global variables
my %orthologs = ();
$|=1; #Autoflush

# Species identifiers
my @species = ('dm6','Dmel','Dsec','Dsim','Dyak','Dere','Dana','Dpse','Dper','Dwil','Dmoj','Dvir','Dgri');
my %spec_id = ();
for my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <genes_list.txt> <dmel_cds.fas> <path_to_genomes>\n";
unless ($ARGV[0]&&$ARGV[1]&&$ARGV[2]) {
	die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $genelist_file = $ARGV[0];
my $cds_seqs_file = $ARGV[1];
my $genomes_dir   = $ARGV[2];
unlink('.log');

# Make directory for sequences
my $seqs_dir = "CDS_seqs";
mkdir($seqs_dir);
my $blast_dir = "$seqs_dir/Raw_blast";
mkdir($blast_dir);
my $orthos_dir = "$seqs_dir/Orthologs";
mkdir($orthos_dir);
my $alignm_dir = "$seqs_dir/Alignments";
mkdir($alignm_dir);

# Get genes list
my $genes_list = get_tab_fields($genelist_file,1);
# Create a dict with gene accessions and names
my %gene_names = ();
# Go through each gene in genes list
for my $gene_i (sort keys %{$genes_list}) {
	my $nam = $genes_list->{$gene_i}->[0];
	my $acc = $genes_list->{$gene_i}->[1];
	$gene_names{$acc} = $nam;
}

# Get cds sequences
my $cds_seqs = get_fasta_seqs($cds_seqs_file);
# Exchange long fasta header with gene accession number
for my $cds (sort keys %{$cds_seqs}) {
    my($acc) = ($cds =~ /FlyBase_Annotation_IDs\:(\w+)/);
    $cds_seqs->{$acc} = $cds_seqs->{$cds};
    delete($cds_seqs->{$cds});
}

# Get genome files list
opendir my $dir, $genomes_dir or die "Cannot open directory: $!";
my @genome_files = readdir $dir;
@genome_files = grep {$_ !~ /^\./} @genome_files;
closedir($dir);

# Get genome sequences
my %genome_seqs = ();
# Go through each genome in genomes list
for my $genome_file (@genome_files) {
	# Get genome species name
	my($g_sp) = ($genome_file =~ /^([^\.]+)\./);
	# Get genome sequence
	$genome_seqs{$g_sp} = get_fasta_seqs("$genomes_dir/$genome_file", 1);
}

# # Blast each gene on each genome
# Go through each gene in genes list
for my $gene_i (sort keys %{$genes_list}) {
    # # Extract gene of interest
	# Gene name and accession number
	my $nam = $genes_list->{$gene_i}->[0];
    my $acc = $genes_list->{$gene_i}->[1];
	# Name dna and protein sequence files
    my $gene_file = "$seqs_dir/$acc.cds.fas";
	my $prot_file = "$seqs_dir/$acc.pep.fas";
	# Print dna sequence to file if it does not yet exist
    unless (-s $gene_file) {
        my $out = open_outfile($gene_file);
        print($out ">$acc\n$cds_seqs->{$acc}\n");
        close($out);
    }
	# Print protein sequence to file if it does not yet exist
	unless (-s $prot_file) {
        my $out = open_outfile($prot_file);
		my $pep_seq = dna2peptide($cds_seqs->{$acc});
		$pep_seq =~ s/\*//;
        print($out ">$acc\n$pep_seq\n");
        close($out);
    }
    # # Blast gene on genomes
	# Go through each genome in genomes list
    for my $genome_file (@genome_files) {
		# Name tblastn output files: table and alignment
		my $blast_out = "$blast_dir/$acc.$genome_file.tbn.aln";
		# Run tblastn if the output files do not yet exist
		unless (-s $blast_out) {
        	system("tblastn -query $prot_file -subject $genomes_dir/$genome_file -out $blast_out $tbnopt >>.log 2>&1");
        }
        # Get genome species name
        my($g_sp) = ($genome_file =~ /^([^\.]+)\./);
		# Extract filtered blast hits
        my $blast_hits = tblastn_alignment_hits($blast_out,$max_eval,$min_hit_ident);
		# Get blast hit groups
		my $hit_groups = group_blast_hits($blast_hits,$min_qcov,$min_grp_ident);
		# Split hits with gaps in query
		$hit_groups = split_gapped_hits($hit_groups);
		# Adjust exon boundaries of hit groups and hit alignments
		$hit_groups = adjust_exon_bounds($hit_groups);
		# Remove short hits with much lower identity
		$hit_groups = remove_odd_hits($hit_groups);
		# Get total identities of hit groups
		my $grp_idents = hit_group_identities($hit_groups);
		# # Save groups as orthologs
		# Initialize ortholog id number
		my $i = 0;
		# Go through each group
		for my $gid (sort {$grp_idents->{$b} <=> $grp_idents->{$a}} keys %{$hit_groups}) {
			# Increment ortholog number
			$i++;
			# Set ortholog id
			my $ort_id = "$acc.$g_sp.$i";
			# Intron length variable
			my $intron_len = 0;
			my $int_beg;
			my $int_end;
			my $prev_send = 0;
			# Go through each hit in group
			for my $hit (sort {$a->[6] <=> $b->[6]} @{$hit_groups->{$gid}}) {
				# Blast hit information
				my $strd = $hit->[11] eq 'plus' ? '+' : '-';
				my $sbeg = $hit->[9] < $hit->[10] ? $hit->[9] : $hit->[10];
				my $send = $hit->[9] < $hit->[10] ? $hit->[10] : $hit->[9];
				my $schr = $hit->[1];
				my $iden = $hit->[2];
				my $eval = $hit->[3];
				my $qlen = $hit->[5];
				my $qbeg = $hit->[6];
				my $qend = $hit->[7];
				my $scor = $hit->[13];
				my $gide = $grp_idents->{$gid};
				# Get intron length
				if ($prev_send) {
					if ($strd eq '+') {
						$int_beg = $prev_send + 1;
						$int_end = $hit->[9] - 1;
					}
					else {
						$int_beg = $hit->[9] + 1;
						$int_end = $prev_send - 1;
					}
					$intron_len = $int_end-$int_beg+1;
				}
				# Get alignment for hit
				my $qaa_seq = $hit->[16];
				# Get subject amino acid sequence for hit
				my $saa_seq = $hit->[17];
				my $saa_stp = $saa_seq =~ tr/\*//;
				# Save identified ortholog in species of gene
				push(
					@{$orthologs{$acc}{$g_sp}{$ort_id}},
					[$schr,$sbeg,$send,$strd,$qbeg,$qend,$qlen,$iden,$eval,$scor,$saa_stp,$intron_len,$qaa_seq,$saa_seq,$gide]
				);
				# Save previous subject hit end position
				$prev_send = $hit->[10];
			}
		}
    }
}

# Count stop codon (including due to small gap inclusion)
my %stop_counts = ();
# # Print each ortholog protein sequence to file
# Go through each original gene accession
for my $acc (sort keys %orthologs) {
	# Go through each species
    for my $gsp (@species) {
		# Go through each ortholog in this species
        for my $ort (sort keys %{$orthologs{$acc}{$gsp}}) {
			# Concatenate amino acid sequence with potential gaps
			my $prot_seq = '';
			# Last position variable for gaps
			my $last_pos = 0;
			# Go through each exon
            for my $ex (@{$orthologs{$acc}{$gsp}{$ort}}) {
				# Amino acid positions
				my $qbeg = $ex->[4];
				my $qend = $ex->[5];
				# Leading gap length
				my $gap_len = $qbeg - $last_pos - 1;
				# Add gap
				$prot_seq .= '-' x $gap_len;
				# Add exon piece to sequence
				$prot_seq .= $ex->[-1];
				# Get new last pos
				$last_pos = $qend;
            }
			# Add potential gap at end of sequence
			my $tot_qlen = $orthologs{$acc}{$gsp}{$ort}[0][6];
			my $gap_len = $tot_qlen - $last_pos;
			# Add final gap
			$prot_seq .= '-' x $gap_len;
			# Open output file
			my $outfile_aa = "$orthos_dir/$ort.pep.fas";
			my $out_aa = open_outfile($outfile_aa);
			# Fasta header
			print($out_aa ">$ort [$gene_names{$acc}]\n");
			# Amino acid sequence
			print($out_aa "$prot_seq\n");
			# Close file
			close($out_aa);
			# Concatenate amino acid sequence with potential gaps
			my $dna_seq = '';
			# Last position variable for gaps
			$last_pos = 0;
			# Go through each exon
            for my $ex (@{$orthologs{$acc}{$gsp}{$ort}}) {
				# CDS positions
				my $sctg = $ex->[0];
				my $sbeg = $ex->[1];
				my $send = $ex->[2];
				# Amino acid positions
				my $qbeg = $ex->[4];
				my $qend = $ex->[5];
				# Intron length
				my $ilen = $ex->[11];
				# Intron length too small for real intron
				if ($ilen > 0 && $ilen < 40) {
					# Get fake intron DNA sequence
					my $int_seq = substr($genome_seqs{$gsp}->{$sctg}, $sbeg-$ilen-1, $ilen);
					if ($acc eq 'CG5920' && $gsp eq 'Dana') {
						print("$ilen $int_seq\n");
					}
					# Add fake intron piece to sequence
					$dna_seq .= $int_seq;
				}
				# Leading gap length
				my $gap_len = $qbeg - $last_pos - 1;
				# Add gap
				$dna_seq .= '-' x ($gap_len*3);
				# Get exon DNA sequence
				my $exn_seq = substr($genome_seqs{$gsp}->{$sctg}, $sbeg-1, $send-$sbeg+1);
				# Get reverse complement if on minus strand
				if ($ex->[3] eq '-') {
					$exn_seq = rev_com($exn_seq);
				}
				# Add exon piece to sequence
				$dna_seq .= $exn_seq;
				# Get new last pos
				$last_pos = $qend;
            }
			# Open output file
			my $outfile_nt = "$orthos_dir/$ort.dna.fas";
			my $out_nt = open_outfile($outfile_nt);
			# Fasta header
			print($out_nt ">$ort [$gene_names{$acc}]\n");
			# DNA sequence
			print($out_nt "$dna_seq\n");
			# Close file
			close($out_nt);
			# Translate final dna sequence, which might contain additional gap sequences
			my $pep_seq = dna2peptide($dna_seq);
			$pep_seq =~ s/\?/-/g;
			# Open output file
			my $outfile_aan = "$orthos_dir/$ort.pepn.fas";
			my $out_aan = open_outfile($outfile_aan);
			# Fasta header
			print($out_aan ">$ort [$gene_names{$acc}]\n");
			# Amino acid sequence
			print($out_aan "$pep_seq\n");
			# Close file
			close($out_aan);
			# Count stop codons
			$stop_counts{$acc}{$gsp}{$ort} = ($pep_seq =~ tr/*/*/);
        }
    }
}

# # Build alignments between original sequence and each best hit
# Go through each original gene accession
for my $acc (sort keys %orthologs) {
	# Name dna and protein sequence files
	my @gene_files = ("$seqs_dir/$acc.cds.fas");
	my @prot_files = ("$seqs_dir/$acc.pep.fas");
	# Go through each species
    for my $gsp (@species) {
		# Go through each ortholog in this species
        for my $ort (sort keys %{$orthologs{$acc}{$gsp}}) {
			# Skip if not best hit
			next unless $ort =~ /\.1$/;
			# Push ortholog files to lists
			push(@gene_files,"$orthos_dir/$ort.dna.fas");
			push(@prot_files,"$orthos_dir/$ort.pep.fas");
		}
	}
	# Output files
	my $outfile_nt = "$alignm_dir/$acc.1.dna.fas";
	my $outfile_aa = "$alignm_dir/$acc.1.pep.fas";
	# Concatenate files
	system("cat @gene_files > $outfile_nt");
	system("cat @prot_files > $outfile_aa");
}

# # Main output
# Open output file
my $outfile1 = "homologous_genes.txt";
my $out1 = open_outfile($outfile1);
# Print title line of output file
print($out1 "#acc\tspec\tgr_id\tchr\tmul\tbeg\tend\tstr\taa_beg\taa_end\taa_len\tident\teval\tscore\tn_stop\tn_stop_n\tin_len\ttot_id\tgt|ag\tn_exons\n\n");
# Go through each original gene accession
for my $acc (sort keys %orthologs) {
	# Title with gene accession, name and protein length
	my $prot_len = (length($cds_seqs->{$acc})/3)-1;
	print($out1 "### $acc - $gene_names{$acc} - ${prot_len}aa\n");
	# Go through each species
    for my $gsp (@species) {
		# Print species name
		print($out1 "# $gsp\n");
		# Go through each ortholog in this species
        for my $ort (sort keys %{$orthologs{$acc}{$gsp}}) {
			# Number of exons
            my $n_exs = @{$orthologs{$acc}{$gsp}{$ort}};
			# Intron bounds
			my %sp_sites = ();
			# Go through each exon
            for my $ex_i (0..$n_exs-1) {
				# Up-/downstream exons
				my $us_ex = $orthologs{$acc}{$gsp}{$ort}[$ex_i];
				my $ds_ex = $orthologs{$acc}{$gsp}{$ort}[$ex_i+1];
				# First exon
				$sp_sites{$us_ex} = "-" unless $sp_sites{$us_ex};
				last if $ex_i == $n_exs-1;
				# Intron coordinates
				my $i_ctg = $us_ex->[0];
				my $i_beg = $us_ex->[3] eq '+' ? $us_ex->[2]+1 : $ds_ex->[2]+1;
				my $i_end = $us_ex->[3] eq '+' ? $ds_ex->[1]-1 : $us_ex->[1]-1;
				# Get Up-/downstream splice sites
				my $us_ss = substr($genome_seqs{$gsp}->{$i_ctg}, $i_beg-1, 2);
				my $ds_ss = substr($genome_seqs{$gsp}->{$i_ctg}, $i_end-2, 2);
				my $sp_site = $us_ex->[3] eq '+' ? "$us_ss|$ds_ss" : rev_com("$us_ss|$ds_ss");
				# Splice sites are not canonical
				if ($sp_site ne 'GT|AG') {
					# Shift intron border coordinates
					$i_beg += 1;
					$i_end += 1;
					# Get Up-/downstream splice sites
					$us_ss = substr($genome_seqs{$gsp}->{$i_ctg}, $i_beg-1, 2);
					$ds_ss = substr($genome_seqs{$gsp}->{$i_ctg}, $i_end-2, 2);
					my $sp_site_p1 = $us_ex->[3] eq '+' ? "$us_ss|$ds_ss" : rev_com("$us_ss|$ds_ss");
					# Shift intron border coordinates
					$i_beg -= 2;
					$i_end -= 2;
					# Get Up-/downstream splice sites
					$us_ss = substr($genome_seqs{$gsp}->{$i_ctg}, $i_beg-1, 2);
					$ds_ss = substr($genome_seqs{$gsp}->{$i_ctg}, $i_end-2, 2);
					my $sp_site_m1 = $us_ex->[3] eq '+' ? "$us_ss|$ds_ss" : rev_com("$us_ss|$ds_ss");
					# Check if any 1-off shifted intron has GT-AG borders
					if ($sp_site_p1 eq 'GT|AG') {
						$sp_site = $sp_site_p1;
						if ($us_ex->[3] eq '+') {
							$us_ex->[2] += 1;
							$ds_ex->[1] += 1;
						}
						else {
							$ds_ex->[2] += 1;
							$us_ex->[1] += 1;
						}
					}
					elsif ($sp_site_m1 eq 'GT|AG') {
						$sp_site = $sp_site_m1;
						if ($us_ex->[3] eq '+') {
							$us_ex->[2] -= 1;
							$ds_ex->[1] -= 1;
						}
						else {
							$ds_ex->[2] -= 1;
							$us_ex->[1] -= 1;
						}
					}
				}
				$sp_sites{$ds_ex} = $sp_site;
			}
			# Go through each exon
            for my $ex (@{$orthologs{$acc}{$gsp}{$ort}}) {
				# Print information about exon
                print($out1 "$acc\t$gsp\t$ort\t");
				for my $ex_fi (0 .. scalar(@{$ex})-1) {
					if ($ex_fi == 7) {
						printf($out1 "%0.1f%%\t", $ex->[$ex_fi]);
					} 
					elsif ($ex_fi == 10) {
						print($out1 "$ex->[$ex_fi]\t$stop_counts{$acc}{$gsp}{$ort}\t");
					}
					elsif ($ex_fi == 12 || $ex_fi == 13) {
						next;
					}
					elsif ($ex_fi == 14) {
						printf($out1 "%0.1f%%\t", $ex->[$ex_fi]);
					}
					else {
						print($out1 "$ex->[$ex_fi]\t");
					}
				}
				print($out1 "$sp_sites{$ex}\t$n_exs\n");
            }
            print($out1 "\n");
        }
    }
	print($out1 "\n");
}
# Add information on muller elements to main output
system("perl add_muller_info.pl homologous_genes.txt muller_shares");

exit;

################################# subroutines #################################

sub get_tab_fields {
	# Take name of tab file
	my($infile,$skip_header) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	if ($skip_header) { shift(@in_data); }
	# Global tree hashes for each species
	my %data_fields = ();
	# Set 0 as start index
	my $id_i = 0;
	# Go through file data
	for my $line (@in_data) {
		# Get line data
        my @d = split(/\t/,$line);
        # Save data fields
        @{$data_fields{$id_i}} = @d;
		$id_i++;
    }
	# Return data fields
	return \%data_fields;
}

sub tblastn_alignment_hits {
	# Take infile name
	my($file,$max_e_val,$prc_ident) = @_;
	# Get file data
	my @file_data = get_file_data_array($file);
    # Initialize variables
    my $hid = 0;
    my %begs_q = ();
    my %ends_q = ();
    my %begs_s = ();
    my %ends_s = ();
	my %nams_s = ();
	my %strs_s = ();
	my %bscore = ();
	my %evalue = ();
	my %idents = ();
	my %a_gaps = ();
	my %s_lens = ();
	my %a_lens = ();
	my %mismas = ();
	my %aln_seq_q = ();
    my %aln_seq_s = ();
    my %blast_hits = ();
	my $qry = '';
    my $sbj = '';
	my $qln = 0;
	my $sln = 0;
    # Parse blast alignment file
    for my $line (@file_data) {
		# Query sequence name
		if ($line =~ /^Query= /) {
			# Get strand from frame orientation
			($qry) = ($line =~ /^Query= (.*)/);
		}
		# Query/Subject length
        elsif ($line =~ /^Length=/) {
			# First occurrence: query
			if (!$qln) {
				($qln) = ($line =~ /^Length=(\d+)/);
			}
			# Later occurrences: subjects
			else {
				($sln) = ($line =~ /^Length=(\d+)/);
			}
        }
        # Subject sequence name
        elsif ($line =~ /^> /) {
			($sbj) = ($line =~ /^> (\S+)/);
        }
		# Score and evalue
		elsif ($line =~ /^ Score = /) {
			# Increment hit id
			$hid++;
			# Save subject name
			$nams_s{$hid} = $sbj;
			# Get Score
            my($bit) = ($line =~ /^ Score = (\S+) bits/);
			# Save bit score
			$bscore{$hid} = $bit;
			# Get evalue
            my($eva) = ($line =~ / Expect.* = (\S+), Method/);
			# Save evalue
			$evalue{$hid} = $eva;
			# Save subject length
			$s_lens{$hid} = $sln;
		}
		# Identity and gaps
		elsif ($line =~ /^ Identities = /) {
			# Get identity percentage
            my($ide) = ($line =~ /^ Identities = \d+\/\d+ \((\d+)\%\), Pos/);
			# Save identity
			$idents{$hid} = $ide;
			# Get matches and alignment length
            my($mas) = ($line =~ /^ Identities = (\d+)\/\d+ \(\d+\%\), Pos/);
			my($aln) = ($line =~ /^ Identities = \d+\/(\d+) \(\d+\%\), Pos/);
			# Save alignment length
			$a_lens{$hid} = $aln;
			# Get gaps
            my($gap) = ($line =~ /, Gaps = (\d+)\/\d+ /);
			# Save identitygaps
			$a_gaps{$hid} = $gap;
			# Save number of mismatches
			$mismas{$hid} = $aln - $mas - $gap;
		}
		# Frame orientation
		elsif ($line =~ /^ Frame/) {
			# Get strand from frame orientation
            my($str) = ($line =~ /^ Frame = (\D)\d/);
			# Save hit strand
			$strs_s{$hid} = $str eq '+' ? 'plus' : 'minus';
        }
		# Query sequence
        if ($line =~ /^Query  /) {
			# Split the line: beg, seq, end
            my @d = split(/\s+/, $line);
            my $beg = $d[1];
            my $seq = $d[2];
            my $end = $d[3];
			# Add sequence to query alignment hash
            $aln_seq_q{$hid} .= uc($seq);
			# Add coordinates to query coordinates hash
            push(@{$begs_q{$hid}},$beg);
            push(@{$ends_q{$hid}},$end);
        }
		# Subject sequence
        if ($line =~ /^Sbjct  /) {
			# Split the line: beg, seq, end
            my @d = split(/\s+/, $line);
            my $beg = $d[1];
            my $seq = $d[2];
            my $end = $d[3];
			# Add sequence to subject alignment hash
            $aln_seq_s{$hid} .= uc($seq);
			# Add coordinates to subject coordinates hash
            push(@{$begs_s{$hid}},$beg);
            push(@{$ends_s{$hid}},$end);
        }
    }
	# Go through each hit
	for my $hit (sort {$a <=> $b} keys %aln_seq_q) {
		# Filter hits
		if ($idents{$hit} < $prc_ident) { next; }
		if ($evalue{$hit} > $max_e_val) { next; }
		# Get hit coordinates on query
		my $qbeg = $begs_q{$hit}[0];
		my $qend = $ends_q{$hit}[-1];
		# Get query coverage per hit
		my $qcov = int((($qend - $qbeg + 1) / $qln * 100)+0.5);
		# Save query/subject coordinates, sequence and strand
		push (@{$blast_hits{$nams_s{$hit}}}, [
			$qry,				# Query id
			$nams_s{$hit},		# Subject id
			$idents{$hit},		# % identity
			$evalue{$hit},		# evalue
			$qcov,				# % query coverage per subject (qcov per hit)
			$qln,				# query length
			$begs_q{$hit}[0],	# q. start
			$ends_q{$hit}[-1],	# q. end
			$s_lens{$hit},		# subject length
			$begs_s{$hit}[0],	# s. start
			$ends_s{$hit}[-1],	# s. end
			$strs_s{$hit},		# subject strand
			$a_lens{$hit},		# alignment length
			$bscore{$hit},		# bit score
			$mismas{$hit},		# mismatches
			$a_gaps{$hit},		# gap opens (total length of gaps)
			$aln_seq_q{$hit},	# q. alignment sequence
			$aln_seq_s{$hit}	# s. alignment sequence
		]);
	}
	# Return blast_hits hash
	return \%blast_hits;
}

sub group_blast_hits {
	# Take gene blast hits
	my($blast_hits,$min_qcov,$min_grp_ident) = @_;
	# Option: maximum distance between hits for initial grouping
	my $max_distance = 5000;
	# # Group hits on same strand on chromosome with small distance
	# Variable for initial grouped hits
	my %init_hit_groups = ();
	# Initial group id
	my $init_group = 0;
	# Go through each contig
	for my $chr (sort keys %{$blast_hits}) {
		# Initialize variables
		my $prev_hit_end = 0;
		# Go through each hit sorted by start position on subject
		for my $hit (sort {$a->[9] <=> $b->[9]} @{$blast_hits->{$chr}}) {
			# Hit coordinates
			my $hit_beg = $hit->[9] < $hit->[10] ? $hit->[9] : $hit->[10];
			my $hit_end = $hit->[9] < $hit->[10] ? $hit->[10] : $hit->[9];
			my $sstrand = $hit->[11];
			# Get distance
			my $distance = $hit_beg-$prev_hit_end;
			# If distance above threshold open new group
			if ($distance > $max_distance) { $init_group++ }
			# Allocate hit to group
			push(@{$init_hit_groups{$chr}{$init_group}{$sstrand}},$hit);
			# Save end position for next hit
			$prev_hit_end = $hit_end;
		}
	}
	# # Split grouped hits if query is covered multiple times
	# Final groups
	my %groups = ();
	# Init group id
	my $gid = 0;
	# Go through each contig
	for my $chr (sort keys %init_hit_groups) {
		# Go through hit group
        for my $igr (sort keys %{$init_hit_groups{$chr}}) {
            for my $str (sort keys %{$init_hit_groups{$chr}{$igr}}) {
				# New group
				$gid++;
				# Sort grouped hits by subject pos according to strand
				if ($str eq 'plus') {
					# Sort plus strand hits by subject pos in ascending order
					@{$init_hit_groups{$chr}{$igr}{$str}} = sort {$a->[9] <=> $b->[9]} @{$init_hit_groups{$chr}{$igr}{$str}};
				}
				else {
					# Sort minus strand hits by subject pos in descending order
					@{$init_hit_groups{$chr}{$igr}{$str}} = sort {$b->[9] <=> $a->[9]} @{$init_hit_groups{$chr}{$igr}{$str}};
				}
				# Initialize group query positional count
				my %qcov_pos = ();
				# Go through each hit of current group
                for my $hit (@{$init_hit_groups{$chr}{$igr}{$str}}) {
					# Get properties
                    my $qbeg = $hit->[6];
                    my $qend = $hit->[7];
                    my $qlen = $qend-$qbeg+1;
					my $total_qlen = $hit->[5];
					# Group query coverage length
					my $gr_qcov_len = keys %{$qcov_pos{$gid}};
					# Check overlap length between hit and grouped hits
					my $overlap_len = 0;
					for my $pos ($qbeg..$qend) {
						if ($qcov_pos{$gid}{$pos}) {
							$overlap_len++;
						}
					}
					# If overlap above threshold open new group
					if ($gr_qcov_len) {
						if ($overlap_len/$gr_qcov_len > 0.95 || $overlap_len/$qlen > 0.95) {
							$gid++;
						}
					}
					# Allocate hit to group
					push(@{$groups{$gid}},$hit);
					# Add to group query positional count
					for my $pos ($qbeg..$qend) { $qcov_pos{$gid}{$pos}++; }
				}
			}
		}
	}
	# # Filter groups with below thershold query coverage
	my @delete_groups = ();
	# Go through each final group
	for $gid (sort keys %groups) {
		# Initialize group query positional count
		my %qcov_pos = ();
		my $total_qlen = 0;
		my $acc = '';
		# Go through each hit of current group
		for my $hit (@{$groups{$gid}}) {
			# Get coordinates
			my $qbeg = $hit->[6];
			my $qend = $hit->[7];
			# Add to group query positional count
			for my $pos ($qbeg..$qend) { $qcov_pos{$pos}++; }
			# Save total query length
			$total_qlen = $hit->[5];
			$acc = $hit->[0];
		}
		# Group query coverage length
		my $gr_qcov_len = keys %qcov_pos;
		# Group query coverage (%)
		my $gr_qcov = $gr_qcov_len/$total_qlen*100;
		# Mark groups with below threshold query coverage for deletion
		if ($gr_qcov < $min_qcov) {
			push(@delete_groups, $gid);
		}
	}
	# Delete marked groups
	for $gid (@delete_groups) {
		delete($groups{$gid});
	}
	# Return final hit groups
	return \%groups;
}

# Input: two strings (i.e. DNA sequences)
# Output: number of mismatches between strings
sub hamming_distance {
	# Get two sequences
	my($p, $q) = @_;
	# Get bitwise distance between sequences
	my $hamming_dist = ($p ^ $q) =~ tr/\0//c;
	# Return the hamming distance
	return $hamming_dist;
}

# Input: hit groups hash
# Output: modified hit groups hash
sub split_gapped_hits {
	# Arguments
	my($hit_groups) = @_;
	# Go through each group
	for my $gid (sort {$a <=> $b} keys %{$hit_groups}) {
		# Array of hits to be deleted after split
		my @hits_to_delete = ();
		# Go through each hit in group
		for my $hit (sort {$a->[6] <=> $b->[6]} @{$hit_groups->{$gid}}) {
			# Blast hit information
			my $strd = $hit->[11] eq 'plus' ? '+' : '-';
			my $sbeg = $hit->[9] < $hit->[10] ? $hit->[9] : $hit->[10];
			my $send = $hit->[9] < $hit->[10] ? $hit->[10] : $hit->[9];
			my $qbeg = $hit->[6];
			my $qend = $hit->[7];
			my $qseq = $hit->[16];
			my $sseq = $hit->[17];
			# Skip if query sequence contains no gaps in alignment
			if ($qseq !~ /-/) { next; }
			# Split alignment at query gaps
			# Initialize sequence and coordinate variables for split hits
			my @qseqs = ();
			my @sseqs = ();
			my @qbegs = ();
			my @sbegs = ();
			$qbegs[0] = $qbeg;
			$sbegs[0] = $hit->[9];
			my @qends = ();
			my @sends = ();
			$qends[0] = $qbeg-1;
			$sends[0] = $strd eq '+' ? $hit->[9]-1 : $hit->[9]+1;
			# Split index
			my $i = 0;
			my $in_gap = 0;
			# Go through positions in query alignment sequence
			for (my $pos=0; $pos<length($qseq); $pos++) {
				# Get character at current alignment position
				my $qaa = substr($qseq, $pos, 1);
				my $saa = substr($sseq, $pos, 1);
				# Non-gap character in query: extend current sub-hit
				if ($qaa ne '-') {
					# Add character to string in seqs array
					$qseqs[$i] .= $qaa;
					$sseqs[$i] .= $saa;
					# Update query hit end position
					$qends[$i] = $qends[$i-1] if !$qends[$i];
					$qends[$i] += 1;
					# Update query hit beg position
					$qbegs[$i] = $qends[$i] if !$qbegs[$i];
					# Update subject hit beg position
					if ($strd eq '+') {
						$sbegs[$i] = $sends[$i]+1 if !$sbegs[$i];
					} else {
						$sbegs[$i] = $sends[$i]-1 if !$sbegs[$i];
					}
					# In gap: false
					$in_gap = 0;
				}
				# First gap character in query: break and create next sub-hit
				elsif ($qaa eq '-' && !$in_gap) {
					# Just entered gap, increment array index
					$i++;
					# In gap: true
					$in_gap = 1;
				}
				# Non-gap character in subject
				if ($saa ne '-') {
					# Update subject hit end position
					$sends[$i] = $sends[$i-1] if !$sends[$i];
					$sends[$i] += 3 if $strd eq '+';
					$sends[$i] -= 3 if $strd eq '-';
					#print("$sbegs[$i] $sends[$i]\n");
				}
			}
			# Adjust last position if on minus strand
			#$sends[-1]++ if $strd eq '-';
			# Go through each index of split alignment array
			for my $i (0..$#qseqs) {
				# Adjust sequence identity
				my $aln_mms = hamming_distance($qseqs[$i], $sseqs[$i]);
				my $aln_len = length($qseqs[$i]);
				my $aln_ide = round_float(($aln_len-$aln_mms)/$aln_len*100, 0);
				my $aln_qcv = round_float($aln_len/$hit->[5]*100, 0);
				my $gap_len = $qseqs[$i] =~ tr/-//;
				# Save split alignment as new hit to group
				push (@{$hit_groups->{$gid}}, [
					$hit->[0],	# Query id
					$hit->[1],	# Subject id
					$aln_ide,	# % identity
					$hit->[3],	# evalue
					$aln_qcv,	# % query coverage per subject (qcov per hit)
					$hit->[5],	# query length
					$qbegs[$i],	# q. start
					$qends[$i],	# q. end
					$hit->[8],	# subject length
					$sbegs[$i],	# s. start
					$sends[$i],	# s. end
					$hit->[11],	# subject strand
					$aln_len,	# alignment length
					$hit->[13],	# bit score
					$aln_mms,	# mismatches
					$gap_len,	# gap opens (total length of gaps)
					$qseqs[$i],	# q. alignment sequence
					$sseqs[$i]	# s. alignment sequence
				]);
			}
			# Add hits to array of hits to be deleted
			push(@hits_to_delete, $hit);
		}
		# Delete original hit that were previously split
		for my $hit (@hits_to_delete) {
			my $index = 0;
			$index++ until $hit_groups->{$gid}->[$index] eq $hit;
			splice(@{$hit_groups->{$gid}}, $index, 1);
		}
	}
	# Return modified hit groups hash
	return $hit_groups;
}

# Input: hit groups hash, hit alignments hash
# Output: modified hit groups hash, modified hit alignments
sub adjust_exon_bounds {
	# Arguments
	my($hit_groups) = @_;
	# Go through each group
	for my $gid (sort {$a <=> $b} keys %{$hit_groups}) {
		# Bool for corrected exon bounds
		my $bounds_corrected = 0;
		# Loop counter
		my $loop_i = 0;
		# While exon bounds are not corrected, try to adjust
		while (!$bounds_corrected) {
			# Adjust exon boundaries for current group
			$hit_groups->{$gid} = adjust_exon_bounds_group($hit_groups->{$gid});
			# Test for correct exon boundaries
			$bounds_corrected = test_exon_bounds_group($hit_groups->{$gid});
			# Count loop iteration
			$loop_i++;
			# Break if maxed out loop count
			last if $loop_i > 10;
		}
	}
	# Return modified hit groups and hit alignment hashes
	return $hit_groups;
}

sub test_exon_bounds_group {
	# Arguments
	my($hit_group) = @_;
	# Bool for corrected exon bounds
	my $bounds_corrected = 1;
	# Sort hits in group by query position
	@{$hit_group} = sort {$a->[6] <=> $b->[6]} @{$hit_group};
	# Loop through indeces of hits in group
	for my $i (0 .. scalar(@{$hit_group})-2) {
		# Get adjacent hits
		my $hit_x = $hit_group->[$i];
		my $hit_y = $hit_group->[$i+1];
		# Hit query coordinates
		my $qbeg_x = $hit_x->[6];
		my $qend_x = $hit_x->[7];
		my $qbeg_y = $hit_y->[6];
		my $qend_y = $hit_y->[7];
		# Hit subject coordinates
		my $schr_x = $hit_x->[1];
		my $schr_y = $hit_y->[1];
		# Adjust subject coordinates for direction
		my $sbeg_x = $hit_x->[9] < $hit_x->[10] ? $hit_x->[9] : $hit_x->[10];
		my $send_x = $hit_x->[9] < $hit_x->[10] ? $hit_x->[10] : $hit_x->[9];
		my $sbeg_y = $hit_y->[9] < $hit_y->[10] ? $hit_y->[9] : $hit_y->[10];
		my $send_y = $hit_y->[9] < $hit_y->[10] ? $hit_y->[10] : $hit_y->[9];
		# Continue if query coordinates of adjacent hits do not overlap
		if ($qend_x < $qbeg_y) { next; }
		# Get alignments for overlapping hits
		my $qseq_x = $hit_x->[16];
		my $sseq_x = $hit_x->[17];
		my $qseq_y = $hit_y->[16];
		my $sseq_y = $hit_y->[17];
		# Init counts for matches in overlapping alignment region
		my %mms_per_cut_site = ();
		# Overlap coordinates
		my $ol_beg = $qbeg_y-1;
		my $ol_end = $qend_x <= $qend_y ? $qend_x : $qend_y;
		# Overlap length
		my $ol_len = $ol_end-$ol_beg;
		# If overlap detected: Exon bounds not corrected
		$bounds_corrected = 0 if $ol_len;
	}
	# Return bool for corrected exon bounds
	return $bounds_corrected;
}

sub adjust_exon_bounds_group {
	# Arguments
	my($hit_group) = @_;
	# Array of hits to be deleted after trimming (if trimmed to 0)
	my @hits_to_delete = ();
	# Sort hits in group by query position
	@{$hit_group} = sort {$a->[6] <=> $b->[6]} @{$hit_group};
	# Loop through indeces of hits in group
	for my $i (0 .. scalar(@{$hit_group})-2) {
		# Get adjacent hits
		my $hit_x = $hit_group->[$i];
		my $hit_y = $hit_group->[$i+1];
		# Skip hit if marked for deletion
		if (my ($matched) = grep $_ eq $hit_x, @hits_to_delete) { next; }
		if (my ($matched) = grep $_ eq $hit_y, @hits_to_delete) { next; }
		# Hit query coordinates
		my $qbeg_x = $hit_x->[6];
		my $qend_x = $hit_x->[7];
		my $qbeg_y = $hit_y->[6];
		my $qend_y = $hit_y->[7];
		# Hit subject coordinates
		my $schr_x = $hit_x->[1];
		my $schr_y = $hit_y->[1];
		# Adjust subject coordinates for direction
		my $sbeg_x = $hit_x->[9] < $hit_x->[10] ? $hit_x->[9] : $hit_x->[10];
		my $send_x = $hit_x->[9] < $hit_x->[10] ? $hit_x->[10] : $hit_x->[9];
		my $sbeg_y = $hit_y->[9] < $hit_y->[10] ? $hit_y->[9] : $hit_y->[10];
		my $send_y = $hit_y->[9] < $hit_y->[10] ? $hit_y->[10] : $hit_y->[9];
		# Continue if query coordinates of adjacent hits do not overlap
		if ($qend_x < $qbeg_y) { next; }
		# Delete hit y if it is entirely within hit x
		if ($qend_x >= $qend_y) {
			# Add hit to array of hits to be deleted
			push(@hits_to_delete, $hit_y);
			# Skip further overlap processing
			next;
		}
		# Get alignments for overlapping hits
		my $qseq_x = $hit_x->[16];
		my $sseq_x = $hit_x->[17];
		my $qseq_y = $hit_y->[16];
		my $sseq_y = $hit_y->[17];
		# Init counts for matches in overlapping alignment region
		my %mms_per_cut_site = ();
		# Overlap coordinates
		my $ol_beg = $qbeg_y-1;
		my $ol_end = $qend_x;
		# Overlap length
		my $ol_len = $ol_end-$ol_beg;
		# print("$ol_len: $ol_beg - $ol_end\n");
		my $qgaps_x = 0;
		my $qgaps_y = 0;
		# Go though each position in the overlap window
		for my $pos ($ol_beg .. $ol_end) {
			# Get alignment position for hit x
			my $apos_x = $pos-$qbeg_x+1;
			my $qaas_x = substr($qseq_x, $apos_x, 1);
			# Get alignment position for hit y
			my $apos_y = $pos-$qbeg_y+1;
			my $qaas_y = substr($qseq_y, $apos_y, 1);
			# Count gaps in query
			$qgaps_x += 1 if $qaas_x eq '-';
			$qgaps_y += 1 if $qaas_y eq '-';
		}
		# Go though each position in the overlap window
		for my $pos ($ol_beg .. $ol_end) {
			# Get length of each side within the overlap window
			my $lt_len = $pos-$ol_beg;
			my $rt_len = $ol_end-$pos;
			# Get alignment position for hit x
			my $apos_x = $qbeg_y-$qbeg_x+1;
			my $qaas_x = substr($qseq_x, $apos_x-1, $lt_len);
			my $saas_x = substr($sseq_x, $apos_x-1, $lt_len);
			my $lt_mms = hamming_distance($qaas_x, $saas_x);
			# Get alignment position for hit y
			my $apos_y = $pos-$qbeg_y+2;
			my $qaas_y = substr($qseq_y, $apos_y-1, $rt_len);
			my $saas_y = substr($sseq_y, $apos_y-1, $rt_len);
			my $rt_mms = hamming_distance($qaas_y, $saas_y);
			# Count mismatches
			$mms_per_cut_site{$lt_len} = $lt_mms + $rt_mms;
			
		}
		# Get overlap split length with best match on both sides
		my $lt_len = (sort {$mms_per_cut_site{$a} <=> $mms_per_cut_site{$b}} keys %mms_per_cut_site)[0];
		# # Remove overlapping alignment pieces from hits
		# Length to be removed from hits x and y
		my $xd_len = $ol_len - $lt_len;
		my $yd_len = $lt_len;
		# Get number of subject gaps for each side
		my $xd_seq = substr($sseq_x, -$xd_len, $xd_len);
		my $yd_seq = substr($sseq_y, 0, $yd_len);
		my $xd_sgaps = $xd_seq =~ /-/ ? $xd_seq =~ tr/-// : 0;
		my $yd_sgaps = $yd_seq =~ /-/ ? $yd_seq =~ tr/-// : 0;
		if ($hit_x->[0] eq 'CG10071' && $hit_x->[1] eq 'utg000010l') {
			print("x=$xd_len y=$yd_len\n");
			print("$hit_x->[9] $hit_x->[10]\n");
			print("$hit_y->[9] $hit_y->[10]\n");
		}
		# Trim right side of hit x alignment if needed
		if ($xd_len) {
			# Adjust end coordinate of query
			$hit_x->[7] = $hit_x->[7] - $xd_len;
			# Adjust end coordinate of subject
			$hit_x->[10] = $hit_x->[10] - (($xd_len-$xd_sgaps) * 3) if $hit_x->[11] eq 'plus';
			$hit_x->[10] = $hit_x->[10] + (($xd_len-$xd_sgaps) * 3) if $hit_x->[11] eq 'minus';
			# Adjust alignment, remove overlap length from end
			$qseq_x = substr($qseq_x, 0, -$xd_len);
			$sseq_x = substr($sseq_x, 0, -$xd_len);
			$hit_x->[16] = $qseq_x;
			$hit_x->[17] = $sseq_x;
			# Adjust sequence identity
			my $hit_x_mms = hamming_distance($qseq_x, $sseq_x);
			my $hit_x_len = length($qseq_x);
			$hit_x->[2] = $hit_x_len ? int((($hit_x_len-$hit_x_mms)/$hit_x_len*100)*1000)/1000 : 0;
			# Prepare hits for deletion if trimmed to 0
			if ($hit_x_len < 1 || !$hit_x->[2]) {
				# Add hit to array of hits to be deleted
				push(@hits_to_delete, $hit_x);
			}
		}
		# Trim left side of hit y alignment if needed
		if ($yd_len) {
			# Adjust start coordinate of query
			$hit_y->[6] = $hit_y->[6] + $yd_len;
			# Adjust start coordinate of subject
			$hit_y->[9] = $hit_y->[9] + (($yd_len-$yd_sgaps) * 3) if $hit_y->[11] eq 'plus';
			$hit_y->[9] = $hit_y->[9] - (($yd_len-$yd_sgaps) * 3) if $hit_y->[11] eq 'minus';
			# Adjust alignment, remove overlap length from start
			$qseq_y = substr($qseq_y, $yd_len);
			$sseq_y = substr($sseq_y, $yd_len);
			$hit_y->[16] = $qseq_y;
			$hit_y->[17] = $sseq_y;
			# Adjust sequence identity
			my $hit_y_mms = hamming_distance($qseq_y, $sseq_y);
			my $hit_y_len = length($qseq_y);
			$hit_y->[2] = $hit_y_len ? int((($hit_y_len-$hit_y_mms)/$hit_y_len*100)*1000)/1000 : 0;
			# Prepare hits for deletion if trimmed to 0
			if ($hit_y_len < 1 || !$hit_y->[2]) {
				# Add hit to array of hits to be deleted
				push(@hits_to_delete, $hit_y);
			}
		}
	}
	# Delete original hit that were previously split
	for my $hit (@hits_to_delete) {
		my $index = 0;
		$index++ until $hit_group->[$index] eq $hit;
		splice(@{$hit_group}, $index, 1);
	}
	# Return modified hit group
	return $hit_group;
}

# Remove any hits that are <28 amino acids AND
# differ in %identity from other hits in the same gene by more than 50%
sub remove_odd_hits {
	# Arguments
	my($hit_groups) = @_;
	# Go through each hit group
	for my $gid (sort keys %{$hit_groups}) {
		# Total hit group length
		my $grp_len = 0;
		# Hit identities hash
		my %hit_ides = ();
		# Go through each hit
		for my $hit (@{$hit_groups->{$gid}}) {
			# Hit identity
			$hit_ides{$hit} = $hit->[2];
			# Add to hit group length
			my $hit_len = $hit->[7] - $hit->[6] + 1;
			$grp_len += $hit_len;
		}
		# Init index
		my $i = 0;
		# Index list for deletion
		my @del_indeces = ();
		# Go through each hit
		for my $hit (@{$hit_groups->{$gid}}) {
			# Hit identity and length
			my $hit_ide = $hit->[2];
			my $hit_len = $hit->[7] - $hit->[6] + 1;
			# Hit is smaller than 28aa
			if ($hit_len < $grp_len/2 && $hit_len < 28) {
				# Max hit identity in group
				my @hit_ides = sort {$b <=> $a} values %hit_ides;
				my $max_ide = $hit_ides[0];
				# Hit differs in %identity from other hits by more than 50%
				if ($hit_ide < $max_ide-50) {
					# Save hit index for removal
					push(@del_indeces,$i);
					print("$hit->[0] $hit->[1] # $hit_len $grp_len | $hit_ide $max_ide\n");
				}
			}
			# Increment index number
			$i++;
		}
		# Delete hits in index list
		foreach my $j (reverse(@del_indeces)) {
			splice(@{$hit_groups->{$gid}},$j,1);
		}
	}
	# Return modified hit groups hash
	return $hit_groups;
}

# Input: hit groups hash
# Output: hit groups identity hash
sub hit_group_identities {
	# Arguments
	my($hit_groups) = @_;
	# Hit groups sequence identity
	my %hit_groups_idents = ();
	# Go through each hit group
	for my $gid (sort keys %{$hit_groups}) {
		# Total hit group length
		my $grp_len = 0;
		# Product for total hit group identity calculation
		my $len_x_ident = 0;
		# Go through each hit
		for my $hit (@{$hit_groups->{$gid}}) {
			# Calculate hit length
			my $hit_len = $hit->[7] - $hit->[6] + 1;
			# Add to hit group length
			$grp_len += $hit_len;
			# Add to sum of length * identity products
			$len_x_ident += $hit_len * $hit->[2];
		}
		# Calculate total hit group identity
		$hit_groups_idents{$gid} = $len_x_ident/$grp_len;
	}
	# Return hash with total sequence identity for each group
	return \%hit_groups_idents;
}

sub dna2peptide {

	my($dna,$frame) = @_;

	# Initialize variables
	my $protein = '';
	my $codon;
	if (not $frame) { $frame = 1 }
	if ($frame && $frame>3) { $frame = 1 }
	# Translate each three-base codon into an amino acid, and append to a protein 
	for(my $i=($frame-1); $i < (length($dna) - 2) ; $i += 3) {
   		$codon = substr($dna,$i,3);
    	$protein .= codon2aa($codon);
	}
	return $protein;
}

sub codon2aa {

    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
        return '?';
    }
}

sub rev_com {
	my($seq) = @_;
	
	my $revcom = reverse $seq;
	$revcom =~ tr/ACGTUNacgtun/TGCAANtgcaan/;
	$revcom =~ s/\s//g;
	
	return $revcom;
}

# Round decimal number to desired decimal places
sub round_float {
	# Take value and number of decimal places
	my($val, $n_dec) = @_;
	# Truncate to desired decimal places
	my $trunc_val = int($val*(10**$n_dec)+0.5)/(10**$n_dec);
	# Return truncated value
	return $trunc_val;
}

# Save fasta data as hash
# Usage: my $sequences = get_fasta_seqs($infile);
sub get_fasta_seqs {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $name = '';
	my %sequences = ();
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($name) = ($line =~ />(.*)$/);
			if ($short) { $name =~ s/\s.*// }
		} else {
			$sequences{$name} .= $line;
		}
	}
	return \%sequences;
}

# Open input file
# Usage: my $in = open_infile($infile);
sub open_infile {
	# Take input file name
    my($file) = @_;
    # Open input file
    my $fh;
    if ($file =~ /.gz$/) {
		open($fh, "gunzip -c $file |") or die("Cannot open file '$file': $!\n");
	} else {
    	open($fh, '<', $file) or die("Cannot open file '$file': $!\n");
    }
    # Return filehandle
    return $fh;
}

# Open output file
# Usage: my $out = open_outfile($outfile);
sub open_outfile {
	# Take output file name
    my($file) = @_;
    # Open output file
    open(my $fh, '>', $file) or die("Cannot open file '$file': $!\n");
    # Return filehandle
    return $fh;
}

# Extract file data and save in array
# Usage: my @filedata = get_file_data_array($file);
sub get_file_data_array {
	# Take input file name
    my($file,$ref_opt) = @_;
    my @filedata = ();
    $ref_opt = 0 unless $ref_opt;
	# Open input file
    my $fh = open_infile($file);
	# Extract lines and save in array
    while (my $line = <$fh>) {
    	$line =~ s/\s+$//; #better chomp
    	push(@filedata, $line);
    }
	# Close file
    close($fh) or die("Unable to close: $!\n");
	# Return array containing file data
    if ($ref_opt) {
    	return \@filedata;
    } else {
    	return @filedata;
    }
}