#!/usr/bin/perl
use strict;

use Bio::SearchIO;
use Bio::SeqIO;
use File::Copy;

my $pidentity_threshold = 75;

($#ARGV == 5) || die "\nUsage: $0 <blast-output> <scaffold-file> <curr-protein-file> <curr-cluster-file> <iteration> <out-dir>\n\n";
my ($blast_file, $scaffold_file, $curr_protein_file, $curr_cluster_file, $iteration, $out_dir) = @ARGV;

my $new_genes_file = 	"$out_dir/new-genes.fna";
my $new_proteins_file =	"$out_dir/new-genes.faa";
my $new_cluster_file = 	"$out_dir/new-genes.faa.cluster.uc";

copy($curr_protein_file, $new_proteins_file);

my %new_genes = ();
my %curr_genes = ();
my %seq2size = ();

### Read sequence sizes ###
my $in = new Bio::SeqIO(-file => $scaffold_file);
while (my $seq_obj = $in->next_seq) {
	$seq2size{$seq_obj->display_id} = $seq_obj->length;
}

### Read new proteins ###
my $in = new Bio::SearchIO(-file => $blast_file);

print STDERR "Reading BLAST file ... ";
my $n = 0;
my ($ncomplete, $npartial, $ntotal) = (0, 0, 0);
while (my $result = $in->next_result) {
	$ntotal++;
	($result->query_name =~ /^(.+)\.(\d+)\.(\d+)$/) || die "\nUnexpected sequence name: " . $result->query_name . "\n\n";
	my ($scaf, $region_start, $region_end) = ($1, $2, $3);
	my $hit = $result->next_hit;
	defined($hit) || next;
	my $hsp = $hit->next_hsp;
	defined($hsp) || die "\nHit object defined but hsp obect not: " . $result->query_name . " vs. " . $hit->name . "\n\n";

	# If hit is not good enough - skip
	($hsp->percent_identity >= $pidentity_threshold) || next;

	my ($s, $e, $strand, $seq) = ($region_start+$hsp->start('query')-1, $region_start+$hsp->end('query')-1, $hsp->strand('query'), $hsp->query_string);
	
	# Too short
	next if(($e-$s+1)/3 < 60);
	# If this is a partial gene that is in the middle - skip
	next if(($s>15) && ($e<$seq2size{$scaf}-15) && (($e-$s+1)/3 < 0.9*$hit->length));

	# If the inside end of the scaffold matches a non-edge in the protein - skip it 
	next if(($s>15) && ((($hsp->strand == 1) && $hsp->hit->start > 10) || (($hsp->strand == -1) && ($hsp->hit->end < $hit->length-10))));
	next if(($e<$seq2size{$scaf}-15) && ((($hsp->strand == -1) && $hsp->hit->start > 10) || (($hsp->strand == 1) && ($hsp->hit->end < $hit->length-10))));

	if(($e-$s+1)/3 >= 0.9*$hit->length) {
		$ncomplete++;
	}
	else {
		$npartial++;
	}

	# We are fine, keep the protein
	$seq =~ s/\-//g;
	my $gene_index = 100000*$iteration + scalar(keys %{$new_genes{$scaf}}) + 1;
	$new_genes{$scaf}{$gene_index} = [$s, $e, $strand, $seq, $hit->name, $hsp->percent_identity];
	$n++;
}

print STDERR "ok, $n new proteins identified ($ncomplete complete, $npartial partial). ", ($ntotal-$n), " proteins rejected\nPreparing gene and protein files ... ";

my $in = new Bio::SeqIO(-file => $scaffold_file);
my $outp = new Bio::SeqIO(-file => ">>$new_proteins_file");
my $outg = new Bio::SeqIO(-file => ">$new_genes_file");
while(my $seq_obj = $in->next_seq) {
	foreach my $gene_index (keys %{$new_genes{$seq_obj->display_id}}) {
		my ($s, $e, $strand, $protein_seq, $hit_name, $pidentity) = @{$new_genes{$seq_obj->display_id}{$gene_index}};
		my $dna_seq = $seq_obj->subseq($s, $e);
		my $gc = int(100*($dna_seq =~ tr/GC/GC/)/length($dna_seq))/100;
		my $display_id = $seq_obj->display_id . "_$gene_index";
		my $desc = "# $s # $e # $strand # ;gc_cont=$gc";
		my $new_protein = new Bio::Seq(-display_id => $display_id, -seq => $protein_seq, -desc => $desc);
		$outp->write_seq($new_protein);
		my $new_gene = new Bio::Seq(-display_id => $display_id, -seq => $dna_seq, -desc => $desc);
		$outg->write_seq($new_gene);
		$curr_genes{$hit_name}{"$display_id $desc"} = [length($protein_seq), $pidentity, undef];
	}
}

print STDERR "ok\nRe-writing cluster file ... ";

open(IN, $curr_cluster_file) || die "\nCannot read $curr_cluster_file\n\n";
open(OUT, ">$new_cluster_file") || die "\nCannot write to $new_cluster_file\n\n";
my %clusters = ();
my %cluster2pidentity = ();

while(<IN>) {
	if($_ =~ /^#/) {
		print OUT $_;
		next;
	}

	chomp;
	# H	9	1370	100.0	.	0	0	1370M	mol-32-1605-077866_4 # 2445 # 6554 # 1 # ;gc_cont=0.631	mol-32-1605-015762_1 # 457 # 4566 # -1 # ;gc_cont=0.631
	# S	10	1342	*	*	*	*	*	mol-32-1605-029087_5 # 2115 # 6140 # -1 # ;gc_cont=0.654	*
	my @fs = split(/\t/);	
	my ($gene_id, $start, $end, $strand, $gc_cont_str) = split(/\s+#\s+/, $fs[8]);

	if($_ =~ /^[C]/) {
		$clusters{$fs[1]}{'line'} = $_;
		next;
	}
	elsif($fs[0] eq 'H') {
		$cluster2pidentity{$fs[1]}{'n'}++;
		$cluster2pidentity{$fs[1]}{'p'}+=$fs[3];
	}
	print OUT "$_\n";
	$clusters{$fs[1]}{'n'}++;
	if(exists($curr_genes{$gene_id})) {
		foreach my $new_gene (keys %{$curr_genes{$gene_id}}) {
			$clusters{$fs[1]}{'n'}++;
			$curr_genes{$gene_id}{$new_gene}[2] = $fs[1];
		}
	}
}

# Now write new proteins
foreach my $curr_gene (keys %curr_genes) {
	foreach my $new_gene (keys %{$curr_genes{$curr_gene}}) {
		defined($curr_genes{$curr_gene}{$new_gene}[2]) || die "\nFatal error: no cluster for $curr_gene/$new_gene\n\n";
		print OUT "H\t", $curr_genes{$curr_gene}{$new_gene}[2], "\t", $curr_genes{$curr_gene}{$new_gene}[0], "\t", int(10*$curr_genes{$curr_gene}{$new_gene}[1])/10, "\t.\tN/A\tN/A\tN/A\t$new_gene\t$curr_gene\n";
		$cluster2pidentity{$curr_genes{$curr_gene}{$new_gene}[2]}{'n'}++;
		$cluster2pidentity{$curr_genes{$curr_gene}{$new_gene}[2]}{'p'}+=int(10*$curr_genes{$curr_gene}{$new_gene}[1])/10;
	}
}

# And write cluster headers
foreach my $cluster (sort {$a <=> $b} keys %clusters) {
	my @fs = split(/\t/, $clusters{$cluster}{'line'});
	$fs[2] = $clusters{$cluster}{'n'};
	if($fs[2] > 1) {
		($cluster2pidentity{$cluster}{'n'} > 0) || die "\nError: $cluster seen 0 times ()\n", $cluster2pidentity{$cluster}{'p'}, "\n";
		$fs[3] = int(100*$cluster2pidentity{$cluster}{'p'}/$cluster2pidentity{$cluster}{'n'})/100;
	}
	print OUT join("\t", @fs), "\n";
}
close(OUT);

print STDERR "ok\n";
print STDERR "\nFinished successfully\n";
