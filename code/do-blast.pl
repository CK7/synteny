#!/usr/bin/perl
use strict;
use File::Copy;

use Bio::SeqIO;
use Bio::SearchIO;

my $min_region_size = 150;
my $num_processors = 10;

($#ARGV >= 3) || die "\nUsage: $0 <scaffold-file> <prodigal-protein-file> <clustering-file> <out-directory> [<glob-prev-regions-files>]\n\n";
my $scaffold_file = shift(@ARGV);
my $protein_file = shift(@ARGV);
my $clustering_file = shift(@ARGV);
my $out_dir = shift(@ARGV); 

my $regions_to_check_file = 	"$out_dir/scaffolds-to-check.fna";
my $blast_db = 			"$out_dir/blast/cluster-representatives.faa";
my $blast_output = 		"$out_dir/scaffolds-to-check.fna.vs.cluster-representatives.blast";

mkdir($out_dir);
(-d $out_dir) || die "\nError: could not create $out_dir or that is a file path\n\n";
mkdir("$out_dir/blast") || die "\nError: $out_dir/blast is not expected to be present\n\n";

my %cluster_seeds = ();
open(IN, $clustering_file) || die "\nError: could not read $clustering_file\n\n";
while(<IN>) {
	# C	3537	2	100.0	*	*	*	*	mol-32-1605-053740_8 # 7522 # 7611 # 1 # ;gc_cont=0.511	*
	($_ =~ /^C/) || next;
	my @fs = split(/\t/);
	($fs[8] =~ /^(\S+) #/) || die "\nUnexpected line in $clustering_file: $_\n\n";
	$cluster_seeds{$1} = $fs[2];
}
close(IN);

my $in = new Bio::SeqIO(-file => $protein_file);
my $out = new Bio::SeqIO(-file => ">$blast_db");
while(my $seq_obj = $in->next_seq) {
	if(exists($cluster_seeds{$seq_obj->display_id})) {
		my $new_seq_obj = new Bio::Seq(-display_id => $seq_obj->display_id, -seq => $seq_obj->seq, -desc => $cluster_seeds{$seq_obj->display_id});
		$out->write_seq($new_seq_obj);
		delete($cluster_seeds{$seq_obj->display_id});
	}
}
if(scalar(keys %cluster_seeds) > 0) {
	die "Error: could not find sequences for ", scalar(keys %cluster_seeds), " proteins\n\n";
}

system("formatdb -p T -o T -i $blast_db");

### Load previous checked regions if given so that we know not to check the same thing again ###
my %checked_regions = ();

if(@ARGV > 0) {
	foreach my $prev_regions_to_check_file (glob($ARGV[0])) {
		my $in = new Bio::SeqIO(-file => $prev_regions_to_check_file);
		while(my $seq_obj = $in->next_seq) {
			$checked_regions{$seq_obj->display_id} = 1;
		}
	}
}

### Find coordinates for genes already identified ###
my %scaf2protein = ();
my $in = new Bio::SeqIO(-file => $protein_file);
while(my $seq_obj = $in->next_seq) {
	($seq_obj->display_id =~ /(.+)_\d+$/) || die "\nUnexpected protein name: " . $seq_obj->display_id . "\n\n";
	my $scaf = $1;
	# 43 # 1074 # 1 # ;gc_cont=0.696
	# 20 # 436 # 1 # ;gc_cont=0.520	
	($seq_obj->description =~ /#\s+(\d+)\s+#\s+(\d+)\s+#\s+\-?\d+\s+#\s+/) || die "\nUnexpected description for " . $seq_obj->display_id . ": " . $seq_obj->description . "\n\n";
	my ($s, $e) = ($1, $2);
	$scaf2protein{$scaf}{$s} = $e;
}

$in = new Bio::SeqIO(-file => $scaffold_file);
my $out = new Bio::SeqIO(-file => ">$regions_to_check_file");

### Now generate the sequences ###
while(my $seq_obj = $in->next_seq) {
	# This will make the whole sequence analyzed
	$scaf2protein{$seq_obj->display_id}{$seq_obj->length+1} = $seq_obj->length+1;

	my $prev_e = 0;
	foreach my $s (sort {$a <=> $b} keys %{$scaf2protein{$seq_obj->display_id}}) {
		if(($s-$prev_e-1) >= $min_region_size) {
			my $display_id = $seq_obj->display_id . "." . ($prev_e+1) . "." . ($s-1);
			if(!exists($checked_regions{$display_id})) {
				my $new_seq_obj = new Bio::Seq(-display_id => $display_id, -seq => $seq_obj->subseq($prev_e+1, $s-1));
				$out->write_seq($new_seq_obj);
			}
		}
		$prev_e = $scaf2protein{$seq_obj->display_id}{$s};
	}
}

print STDERR "ok\nBlasting ... ";
system("blastall -p blastx -d $blast_db -i $regions_to_check_file -F F -e 1e-30 -b 1 -v 1 -a $num_processors -o $blast_output");
print STDERR "ok\n";
##system("../3.predict-genes/get-genes-from-blast.pl <blast-output> <scaffold-file> <iteration> <out-dir>");
