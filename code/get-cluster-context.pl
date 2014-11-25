#!/usr/bin/perl
use strict;

($#ARGV >= 1) || die "\nUsage: $0 <seq2claster_file> <cluster> [<cluster2> ... <clustern>]\n\n";
my $seq2cluster_file = shift(@ARGV);

my %target_clusters = ();
my %target_clusters2colors = ();

foreach my $cluster (@ARGV) {
	$target_clusters{$cluster} = 0;
	$target_clusters2colors{$cluster} = 30+scalar(keys %target_clusters2colors);
}

my %seqs = ();

open(IN, $seq2cluster_file) || die "\nCannot read $seq2cluster_file\n\n";
while(<IN>) {
	my @fs = split(/\t/);
	$fs[0] =~ s/_\d+$//;
	if(exists($target_clusters{$fs[1]})) {
		$seqs{$fs[0]}{$fs[1]} = 1;
		$target_clusters{$fs[1]}++;
	}
}
close(IN);

my @candidate_seqs = keys %seqs;
foreach my $seq (@candidate_seqs) {
	(scalar(keys %{$seqs{$seq}}) == scalar(keys %target_clusters)) || delete($seqs{$seq});
}

print "Clusters found:\n";
foreach my $cluster (keys %target_clusters) {
	my $color = $target_clusters2colors{$cluster};
	$_ = "$cluster\t" . $target_clusters{$cluster};
	system("echo -en '\\E[47;$color m'\"\\033[1m$_\\033[0m\"");
	print "\n";
}

$_ = "All together: " . scalar(keys %seqs);
system("echo -e \"\\033[1m$_\\033[0m\"");
print "------------------------------------------------------------------\n";

open(IN, $seq2cluster_file) || die "\nCannot read $seq2cluster_file\n\n";
my $prev = undef;
while(<IN>) {
	my @fs = split(/\t/);
	$fs[0] =~ s/_\d+$//;
	if(exists($seqs{$fs[0]})) {
		chomp;
		print "------------------------------------------------------------------\n" if(defined($prev) && ($prev ne $fs[0]));
		if(exists($target_clusters{$fs[1]})) {
			my $color = $target_clusters2colors{$fs[1]};
			system("echo -en '\\E[47;$color m'\"\\033[1m$_\\033[0m\"");
			print "\n";
#			system("echo -e \"\\033[1m$_\\033[0m\"");
		}
		else {
			print "$_\n";
		}
		$prev = $fs[0];
	}
}
close(IN);
