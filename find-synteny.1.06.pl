#!/usr/bin/perl
use strict;

###############################################################################################################################################################
# v1.01: added the component functionality, including the identification of components with adjacent paths (in addition to linear components) 
# v1.02: modified the definition of adjacent components to also include dead-ends
# v1.06: added strand
###############################################################################################################################################################

BEGIN {
	# Replace the following path to the program's directory 
	unshift(@INC, '<synteny-directory-path>');
}

use Node;
use Component;
use DotWriter;
use ExcelWriter;
use Component;
use Metacomponent;

my $max_adjacent_component_length = 10;	# maximal distance from start to end of components with adjacent paths

($#ARGV == 3) || die "\nUsage: $0 <uclust-file> <connection-threshold> <out-dir> <out-prefix>\n\n";
my ($uclust_file, $connection_threshold, $out_dir, $out_prefix) = @ARGV;

mkdir($out_dir);
mkdir("$out_dir/log-files");

my $graph_file = 	"$out_dir/$out_prefix.graph.dot";
my $excel_file = 	"$out_dir/$out_prefix.excel.txt";
my $nodes_file = 	"$out_prefix.nodes";
my $connections_file = 	"$out_prefix.connections";
my $components_file = 	"$out_prefix.components";
my $synteny_file = 	"$out_prefix.synteny.txt";
my $seq2clusters_file =	"$out_prefix.seq2clusters";

my $min_node_threshold = 2;

my %all_nodes = ();
my %all_sequences = ();
my @all_edges = ();

my %cluster2size = ();
my %cluster2seqs = ();
my %seq2clusters = ();

my %all_metacomponents = ();
my %all_components = ();

read_clustering_file($uclust_file, \%all_nodes, \%all_sequences);
write_seq2cluster(\%all_sequences, "$out_dir/log-files/01.$seq2clusters_file.after-clustering-file-reading.txt");

construct_graph(\%all_nodes, \@all_edges, \%all_sequences);

### Write nodes and connections files ###
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/02.$nodes_file.graph-constructions.txt", "$out_dir/log-files/02.$connections_file.graph-constructions.txt");

resolve_junctions(\%all_nodes);
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/03.$nodes_file.junctions-resolved.txt", "$out_dir/log-files/03.$connections_file.junctions-resolved.txt");
write_seq2cluster(\%all_sequences, "$out_dir/log-files/03.$seq2clusters_file.junctions-resolved.txt");

print STDERR "Removing nodes with less than $min_node_threshold instances ... ";
my ($nbefore_nodes, $nbefore_edges) = (scalar(keys %all_nodes), scalar(@all_edges));
remove_rare_nodes(\%all_nodes, \@all_edges, \%all_sequences, $min_node_threshold);
my ($nafter_nodes, $nafter_edges) = (scalar(keys %all_nodes), scalar(@all_edges));
print STDERR "ok, number of nodes is $nafter_nodes ($nbefore_nodes before), number of edges is $nafter_edges ($nbefore_edges before)\n";
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/04.$nodes_file.lt

$min_node_threshold instances" . "-removed.txt", "$out_dir/log-files/04.$connections_file.lt$min_node_threshold instances" . "-removed.txt");
write_seq2cluster(\%all_sequences, "$out_dir/log-files/04.$seq2clusters_file.lt$min_node_threshold instances" . "-removed.txt");

### Look for components with single in and out ###
print STDERR "Looking for bubbles ... ";

### First scan for components: try to collect everything regradless of whether it is a tip node or not ###
my ($before_ncomponents, $before_nnodes_assigned) = (scalar(keys %all_components), num_nodes_assigned(\%all_nodes));
search_for_bubble_components(\%all_nodes, \%all_components, undef);
my ($after_ncomponents, $after_nnodes_assigned) = (scalar(keys %all_components), num_nodes_assigned(\%all_nodes));
print STDERR "Added ", ($after_ncomponents-$before_ncomponents), " new components by clustering ", ($after_nnodes_assigned-$before_nnodes_assigned), " nodes\n";
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/05.$nodes_file.after-bubbles-all.txt", "$out_dir/log-files/05.$connections_file.after-bubbles-all.txt");
write_components(\%all_components, "$out_dir/log-files/05.$components_file.after-bubbles-all.txt");

### Look for linear components ###
print STDERR "Looking for linear components ... ";
($before_ncomponents, $before_nnodes_assigned) = (scalar(keys %all_components), num_nodes_assigned(\%all_nodes));
search_for_linear_components(\%all_nodes, \%all_components);
($after_ncomponents, $after_nnodes_assigned) = (scalar(keys %all_components), num_nodes_assigned(\%all_nodes));
print STDERR "Added ", ($after_ncomponents-$before_ncomponents), " new components by clustering ", ($after_nnodes_assigned-$before_nnodes_assigned), " nodes\n";
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/06.$nodes_file.after-linear-all.txt", "$out_dir/log-files/06.$connections_file.after-linear-all.txt");
write_components(\%all_components, "$out_dir/log-files/06.$components_file.after-linear-all.txt");

### Second scan for components: avoid tip nodes ###
identify_tip_nodes(\%all_nodes);
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/07.$nodes_file.after-tip-detection.txt", "$out_dir/log-files/07.$connections_file.after-tip-detection.txt");

print STDERR "Looking for bubble (internal nodes only) ... ";
($before_ncomponents, $before_nnodes_assigned) = (scalar(keys %all_components), num_nodes_assigned(\%all_nodes));
search_for_bubble_components(\%all_nodes, \%all_components, 'in');
($after_ncomponents, $after_nnodes_assigned) = (scalar(keys %all_components), num_nodes_assigned(\%all_nodes));
print STDERR "Added ", ($after_ncomponents-$before_ncomponents), " new components by clustering ", ($after_nnodes_assigned-$before_nnodes_assigned), " nodes\n";
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/08.$nodes_file.after-bubbles-all.txt", "$out_dir/log-files/08.$connections_file.after-bubbles-all.txt");
write_components(\%all_components, "$out_dir/log-files/08.$components_file.after-bubbles-all.txt");

### Look for linear components ###
print STDERR "Looking for linear components (internal nodes only) ... ";
($before_ncomponents, $before_nnodes_assigned) = (scalar(keys %all_components), num_nodes_assigned(\%all_nodes));
search_for_linear_components(\%all_nodes, \%all_components, 'in');
($after_ncomponents, $after_nnodes_assigned) = (scalar(keys %all_components), num_nodes_assigned(\%all_nodes));
print STDERR "Added ", ($after_ncomponents-$before_ncomponents), " new components by clustering ", ($after_nnodes_assigned-$before_nnodes_assigned), " nodes\n";
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/09.$nodes_file.after-linear-in.txt", "$out_dir/log-files/09.$connections_file.after-linear-in.txt");
write_components(\%all_components, "$out_dir/log-files/09.$components_file.after-linear-in.txt");

### Assign tip nodes to components ###
print STDERR "Assigning tip nodes to components ... ";
$before_nnodes_assigned = num_nodes_assigned(\%all_nodes);
print STDERR "\n*** Not assigning nodes to components ***\n";
##assign_nodes_to_components(\%all_nodes, \@all_edges, \%all_components, 'tip');
$after_nnodes_assigned = num_nodes_assigned(\%all_nodes);
print STDERR "Added ", ($after_nnodes_assigned-$before_nnodes_assigned), " nodes to existing components\n";
write_connections_and_nodes(\%all_nodes, \@all_edges, "$out_dir/log-files/10.$nodes_file.after-assigning-tip-nodes.txt", "$out_dir/log-files/10.$connections_file.after-assigning-tip-nodes.txt");
write_components(\%all_components, "$out_dir/log-files/10.$components_file.after-assigning-tip-nodes.txt");

### Create meta-components ###
print STDERR "Creating meta-components ... ";
create_meta_components(\%all_components, \%all_metacomponents);
print STDERR (scalar(keys %all_metacomponents)), " meta-components created\n";

### Assign strands to nodes with no current assignment ###
set_strands_to_nodes(\%all_nodes);

### Finally, write output ###
my $writer = new DotWriter;
$writer->write(\%all_metacomponents, \%all_components, \%all_nodes, \@all_edges, $graph_file);

$writer = new ExcelWriter;
$writer->write(\%all_metacomponents, \%all_components, \%all_nodes, \%all_sequences, $excel_file);

print STDERR scalar(keys %all_metacomponents), " meta-components\n";
print STDERR scalar(keys %all_components), " components (including components in meta-components)\n";
print STDERR num_nodes_assigned(\%all_nodes), " nodes were assigned to clusters\n";

###############################################################################################################################################################
sub num_nodes_assigned {
	my $all_nodes_ref = shift;
	my $n = 0;

	foreach my $node_obj (values %{$all_nodes_ref}) {
		$n++ if(defined($node_obj->component(5)) || defined($node_obj->component(3)));
	}
	return $n;
}

###############################################################################################################################################################
sub write_seq2cluster {
	my ($all_sequences_ref, $seq2clusters_file) = @_;
	open(OUT, ">$seq2clusters_file") || die "\nCannot write to $seq2clusters_file\n\n";
	foreach my $seq_obj (values %{$all_sequences_ref}) {
		foreach my $gene_obj (@{$seq_obj->genes}) {
			# mol-32-1605-085889_200001       307     (3, 476, -1)
			print OUT $gene_obj->id, "\t", $gene_obj->node->id, "\t(", $gene_obj->start, ", ", $gene_obj->end, ", ", $gene_obj->strand, ")\n"; 
		}
	}
	close(OUT);
}

###############################################################################################################################################################
sub write_components {
	my ($all_components_ref, $components_file) = @_;
	open(OUT, ">$components_file") || die "\nCannot write to $components_file\n\n";
	foreach my $component_obj (values %{$all_components_ref}) {
		print OUT $component_obj->id, "\t", $component_obj->type, "\t", $component_obj->size, "\t", $component_obj->node1->id, ":", $component_obj->side1, "\t",
				$component_obj->node2->id, ":", $component_obj->side2, "\t";
		my $i = 1;
		if(defined($component_obj->in_node($i))) {
			print OUT $component_obj->in_node($i++)->id;
			while (defined($component_obj->in_node($i))) {
				print OUT ", ", $component_obj->in_node($i++)->id;
			}
		}
		print OUT "\n";
	}
	close(OUT);
}

###############################################################################################################################################################
sub write_connections_and_nodes {
	my ($all_nodes_ref, $all_edges_ref, $nodes_file, $edges_file) = @_;

	open(OUT, ">$nodes_file") || die "\nCannot write to $nodes_file\n\n";
	foreach my $node_obj (values %{$all_nodes_ref}) {
		print OUT $node_obj->id, "\t", $node_obj->type, "\t", $node_obj->size, "\t";
		print OUT $node_obj->component(5)->id if(defined($node_obj->component(5)));
		print OUT "\t";
		print OUT $node_obj->component(3)->id if(defined($node_obj->component(3)));
		print OUT "\n";
	}
	close(OUT);

	open(OUT, ">$edges_file") || die "\nCannot write to $edges_file\n\n";
	foreach my $edge_obj (@{$all_edges_ref}) {
		print OUT $edge_obj->id, "\t", $edge_obj->nsequences, "\n";
	}
	close(OUT);
}

###############################################################################################################################################################
sub read_clustering_file {
	my ($uclust_file, $all_nodes_ref, $all_sequences_ref) = @_;

	# Read clustering file
	# First pass: create all node objects
	open(IN, $uclust_file) || die "\nCannot read $uclust_file\n\n";
	while(<IN>) {
		next if($_ !~ /^C\t/);
		chomp;
		# C	25	4	93.9	*	*	*	*	mol-32-1605-040921_1 # 55 # 7911 # -1 # ;gc_cont=0.580	*
		my @fs = split(/\t/);
		if($fs[2] >= $min_node_threshold) {
			my $node = new Node($fs[1], 'in');
			$all_nodes_ref->{$fs[1]} = $node;
		}
	}
	close(IN);

	# Second pass: get all edges
	my %genes_to_add = ();
	open(IN, $uclust_file) || die "\nCannot read $uclust_file\n\n";
	while(<IN>) {
		next if($_ =~ /^#/);
		chomp;
		# S	25	2619	*	*	*	*	*	mol-32-1605-040921_1 # 55 # 7911 # -1 # ;gc_cont=0.580	*
		# H	25	2560	93.1	.	0	0	55I2000MI19M6I371M3D167M	mol-32-1605-014159_1 # 18 # 7697 # 1 # ;gc_cont=0.584	mol-32-1605-040921_1 # 55 # 7911 # -1 # ;gc_cont=0.580
		# H	25	2112	91.3	.	0	0	391I1158M67I419M46I365M3I170M	mol-32-1605-002772_1 # 89 # 6424 # 1 # ;gc_cont=0.591	mol-32-1605-040921_1 # 55 # 7911 # -1 # ;gc_cont=0.580
		# H	25	37	97.3	.	0	0	2582I37M	mol-32-1605-066139_1 # 61 # 171 # 1 # ;gc_cont=0.577	mol-32-1605-040921_1 # 55 # 7911 # -1 # ;gc_cont=0.580
		# C	25	4	93.9	*	*	*	*	mol-32-1605-040921_1 # 55 # 7911 # -1 # ;gc_cont=0.580	*
		my @fs = split(/\t/);
		next if($fs[0] eq 'C');

		my ($gene_id, $start, $end, $strand, $gc_cont_str) = split(/\s+#\s+/, $fs[8]);
		($gene_id =~ /(.+)_\d+/) || die "\nUnexpected gene name $gene_id in the following line of $uclust_file:\n$_\n\n";
		my $seq = $1;

		$all_sequences_ref->{$seq} = new Sequence($seq) if(!exists($all_sequences_ref->{$seq}));
		my $gene = new Gene($gene_id, $all_sequences_ref->{$seq}, $start, $end, $strand, $all_nodes_ref->{$fs[1]});
		$genes_to_add{$seq}{$start} = $gene;
	}

	foreach my $seq_id (keys %genes_to_add) {
		my $true_index = 0;
		foreach my $start (sort {$a <=> $b} keys %{$genes_to_add{$seq_id}}) {
			$true_index++;
			my $gene_obj = $genes_to_add{$seq_id}{$start};
			next if(!defined($gene_obj->node));
			$all_sequences_ref->{$seq_id}->add_gene($gene_obj, $true_index);
			$gene_obj->node->add_gene($gene_obj);
		}
	}

	# We would like to remove all instances of sequences with only one gene. This will free us from a headache later.
	my @seq_ids = keys %{$all_sequences_ref};
	foreach my $seq_id (@seq_ids) {
		my $seq_obj = $all_sequences_ref->{$seq_id};
		if($seq_obj->ngenes() == 0) {
			delete($all_sequences_ref->{$seq_id});
			next;
		}
		next if($seq_obj->ngenes() != 1);
		my ($gene_obj) = @{$seq_obj->genes};
		defined($gene_obj) || die "\nError EMain0015: no gene objects were found for " . $seq_obj->id . " even though the number of genes is reported to be " . $seq_obj->ngenes() . "\n\n";
		my $node_obj = $gene_obj->node;

		delete($all_sequences_ref->{$seq_id});
		$node_obj->remove_gene($gene_obj);
		delete($all_nodes_ref->{$node_obj->id}) if($node_obj->size == 0);
	}

	close(IN);
}

###############################################################################################################################################################
sub construct_graph {
	my ($all_nodes_ref, $all_edges_ref, $all_sequences_ref) = @_;

	# Construct graph
	my %edges_seen = ();

	foreach my $seq_obj (values %{$all_sequences_ref}) {
		foreach my $i (1 .. $#{$seq_obj->genes}) {
			my ($prev_node_obj, $curr_node_obj) = ($seq_obj->genes->[$i-1]->node, $seq_obj->genes->[$i]->node);

			# We ignore self connections, these are most likely due to broken genes
			next if($prev_node_obj->id eq $curr_node_obj->id);
			my $prev_end = ($seq_obj->genes->[$i-1]->strand == 1)? 3 : 5; 
			my $curr_end = ($seq_obj->genes->[$i]->strand == 1)? 5 : 3; 

			my $edge_obj = new Edge($prev_node_obj, $prev_end, $curr_node_obj, $curr_end);

			if(!exists($edges_seen{$edge_obj->id})) {
				push(@{$all_edges_ref}, $edge_obj);
				$prev_node_obj->add_edge($prev_end, $edge_obj); 
				$curr_node_obj->add_edge($curr_end, $edge_obj); 
				$edges_seen{$edge_obj->id} = $edge_obj;
			}
			else {
				$edge_obj = $edges_seen{$edge_obj->id};
			}

			$edge_obj->add_sequence($seq_obj);
		}
	}
}

###############################################################################################################################################################
# Identify all nodes that lead to deadends. Don't touch nodes that are already assigned to components
###############################################################################################################################################################
sub identify_tip_nodes {
	my ($all_nodes_ref) = @_;

	my $round = 0;
	while($round++ < 10)  {
		my @node_objs = values %{$all_nodes_ref};
		my @tip_node_objs = ();

		foreach my $node_obj (@node_objs) {
			next if(defined($node_obj->component(5)) || defined($node_obj->component(3)));
			my ($ends5_ref, $ends3_ref) = ($node_obj->connected(5, 'in'), $node_obj->connected(3, 'in'));
			push(@tip_node_objs, $node_obj) if((@{$ends5_ref} == 0) || (@{$ends3_ref} == 0));
		}
		foreach my $node_obj (@tip_node_objs) {
			$node_obj->type('tip');
		}
		last if(@tip_node_objs == 0);
	}
}

###############################################################################################################################################################
sub get_node_name {
        my ($group_id, $node) = @_;
        return (10000000+$node*100+$group_id);
}

###############################################################################################################################################################
sub remove_rare_nodes {
	my ($all_nodes_ref, $all_edges_ref, $all_sequences_ref, $min_node_threshold) = @_;

	my $change = 1;
	while ($change) {
		$change = 0;
		# We will update $all_sequences_ref and the genes on each sequence and then rebuild $all_edges_ref and $all_nodes_ref
		my @node_objs = values %{$all_nodes_ref};
		foreach my $node_obj (@node_objs) {
			next if($node_obj->size >= $min_node_threshold);
			foreach my $gene_obj (@{$node_obj->genes}) {
				my $seq_obj = $gene_obj->seq;

				$seq_obj->remove_gene($gene_obj);
			}
			delete($all_nodes_ref->{$node_obj->id});
			$change = 1;
		}

		my @seq_objs = values %{$all_sequences_ref};
		foreach my $seq_obj (@seq_objs) {
			if($seq_obj->ngenes < 2) {
				foreach my $gene_obj (@{$seq_obj->genes}) {
					$gene_obj->node->remove_gene($gene_obj);
				}
				delete($all_sequences_ref->{$seq_obj->id});
				$change = 1;
			}
		}
	}

	# Remove all edges - we are about to create new ones
	splice(@{$all_edges_ref}, 0, scalar(@{$all_edges_ref}));

	# Remove all nodes and create new ones. We do this in order to get rid of the edges inside the nodes. Update reference in the gene objects
	my @node_objs = values %{$all_nodes_ref};
	foreach my $node_obj (@node_objs) {
		my $new_node_obj = new Node($node_obj->id, $node_obj->type);
		foreach my $gene_obj (@{$node_obj->genes}) {
			$gene_obj->node($new_node_obj);
			$new_node_obj->add_gene($gene_obj);
		}
		delete($all_nodes_ref->{$node_obj->id});
		$all_nodes_ref->{$new_node_obj->id} = $new_node_obj;
	}

	construct_graph($all_nodes_ref, $all_edges_ref, $all_sequences_ref);

	# Sanity check
	foreach my $node_id (keys %{$all_nodes_ref}) {
		defined($all_nodes_ref->{$node_id}) || die "\nError EMain0025: node $node_id is undef\n\n";
		($all_nodes_ref->{$node_id}->size >= $min_node_threshold) || die "\nError EMain0026: node $node_id has only " . $all_nodes_ref->{$node_id}->size . " instances\n\n"; 
	}

	foreach my $seq_id (keys %{$all_sequences_ref}) {
		my $seq_obj = $all_sequences_ref->{$seq_id};
		defined($seq_obj) || die "\nError EMain0029: sequence $seq_id has no object\n\n"; 
		($seq_obj->ngenes >= $min_node_threshold) || die "\nError EMain0027: sequence " . $seq_obj->id . " has fewer than $min_node_threshold genes\n\n";
		foreach my $gene_obj (@{$seq_obj->genes}) {
			defined($gene_obj) || die "\nError EMain0030: undefined gene object found\n\n";
			my $node_obj = $gene_obj->node;
			defined($gene_obj) || die "\nError EMain0031: undefined node object retrieved from gene " . $gene_obj->id . "\n\n";
			exists($all_nodes_ref->{$gene_obj->node()->id}) || die "\nError EMain0028: node " . $gene_obj->node()->id . " still exists in " . $gene_obj->id . "\n\n";
		}
	}
}

###############################################################################################################################################################
sub resolve_junctions {
	my ($all_nodes_ref) = @_;

	my $change = 1;
	my $nresolved = 0;

	while ($change) {
		$change = 0;
		foreach my $node_id (keys %{$all_nodes_ref}) {
			my $node_obj = $all_nodes_ref->{$node_id};
			# Is $node_obj a junction?
			# $node_obj has to have 2 connections or more on both the 3' and 5' ends
			((@{$node_obj->edges(5)} > 1) && (@{$node_obj->edges(3)} > 1)) || next;

			# These are all $node_obj's neighbors. We will be looking for them alone
			my %relevant_nodes = ($node_obj->id => 1);
			foreach my $s (5, 3) {
				foreach my $edge_obj (@{$node_obj->edges($s)}) {
					my ($node_obj2, $s2) = $edge_obj->other_end($node_obj->id, $s);
					$relevant_nodes{$node_obj2->id} = 1;
				}
			}

			# Now go over all sequences in which $node_obj appears and find connections between its neighbors
			my %connections = ();
			foreach my $seq_obj (@{$node_obj->sequences}) {
				my $node_objs_ref = $seq_obj->nodes;
				foreach my $j (0 .. $#{$node_objs_ref}) {
					my $n1_obj = $node_objs_ref->[$j];
					next if(!exists($relevant_nodes{$n1_obj->id}));
					foreach my $i (($j+1) .. $#{$node_objs_ref}) {
						my $n2_obj = $node_objs_ref->[$i];
						next if(!exists($relevant_nodes{$n2_obj->id}));
						$connections{$n1_obj->id}{$n2_obj->id}++; $connections{$n2_obj->id}{$n1_obj->id}++;
					}
				}
				
			}

			# The reason we keep connections to $node_obj is to give clusters that appear with $c alone a chance to be kept. We don't use
			# $node_obj itself for clustering, of course
			delete($connections{$node_obj->id});

			# Can we split the nodes into 2 or more clusters?
			my %node2group = ();
			my %group2nodes = ();
			my $group_id = 0;
			while(scalar(keys %connections) > 0) {
				my ($curr_n) = keys %connections;
				my %check = ($curr_n => 1);
				$group_id++;
				while(scalar(keys %check) > 0) {
					my ($curr_n) = keys %check;
					delete($check{$curr_n});
					!exists($node2group{$curr_n}) || die "\n Error EMain0001: clustering does not work right!!!\n\n";
					$node2group{$curr_n} = $group_id;
					$group2nodes{$group_id}{$curr_n} = 1;
					foreach my $other_n (keys %{$connections{$curr_n}}) {
						next if($other_n == $node_obj->id);
						if(exists($node2group{$other_n})) {
							($node2group{$other_n} == $group_id) || die "\n EMain0002: clustering does not work right!!!\n\n";
						}
						else {
							$check{$other_n} = 1;
						}
					}
					delete($connections{$curr_n});
				}
			}

			# No success: all clusters are grouped into a single set
			next if($group_id == 1);

			# Hurray, we can split the genes!
			$nresolved++;

			# We will remove $node_obj and create $group_id nodes instead
			delete($all_nodes_ref->{$node_obj->id});
			foreach my $group_id (keys %group2nodes) {
				my $new_node_obj = new Node(get_node_name($group_id, $node_obj->id), $node_obj->type);
				$all_nodes_ref->{$new_node_obj->id} = $new_node_obj;
			}

			# Get the new node for each sequence. No need to update sequence objects as they don't keep nodes
			my %seq_id2node_obj = ();
			foreach my $s (5, 3) {
				foreach my $edge_obj (@{$node_obj->edges($s)}) {
					my ($node_obj2, $s2) = $edge_obj->other_end($node_obj->id, $s);
					exists($node2group{$node_obj2->id}) || die "\nError EMain0003: did not find a group for " . $node_obj2->id . "\n\n";
					my $new_node_id = get_node_name($node2group{$node_obj2->id}, $node_obj->id);
					foreach my $seq_obj (@{$edge_obj->sequences}) {
						(!exists($seq_id2node_obj{$seq_obj->id}) || ($seq_id2node_obj{$seq_obj->id}->id eq $new_node_id)) || die "\nError EMain0004: same sequence assigned to two different new nodes\n\n";
						$seq_id2node_obj{$seq_obj->id} = $all_nodes_ref->{$new_node_id};
					}
				}
			}

			# Update all genes
			foreach my $gene_obj (@{$node_obj->genes}) {
				$gene_obj->node($seq_id2node_obj{$gene_obj->seq->id});
				$seq_id2node_obj{$gene_obj->seq->id}->add_gene($gene_obj);
			}

			# Update all edges
			foreach my $s (5, 3) {
				foreach my $edge_obj (@{$node_obj->edges($s)}) {
					my ($node_obj2, $s2) = $edge_obj->other_end($node_obj->id, $s);
					my $new_node_obj = $all_nodes_ref->{get_node_name($node2group{$node_obj2->id}, $node_obj->id)};

					$edge_obj->replace_end($node_obj, $s, $new_node_obj, $s);
					$new_node_obj->add_edge($s, $edge_obj);
				}
			}
			$change = 1;
			last;
		}
	}
	print STDERR "Resolved $nresolved junctions\n";
}

###############################################################################################################################################################
# This function makes sure that all paths from $a end at $b (both are node ends - either 5' or 3'). Starting from version 1.02 it is considered legal to
# have dead-ends in the component 
###############################################################################################################################################################
sub all_paths_from_a_to_b {
	my ($a_obj, $a_side, $b_obj, $b_side, $type, $component_end_ref) = @_;

	my %check = ($a_obj->id . ":$a_side" => 0);

	my %component_node_objs = ($a_obj->id => $a_obj, $b_obj->id => $b_obj);

	foreach my $e (keys %{$component_end_ref}) {
		delete($component_end_ref->{$e});
	}
	$component_end_ref->{$a_obj->id . ":$a_side"} = 1;
	$component_end_ref->{$b_obj->id . ":$b_side"} = 1;

	while(scalar(keys %check) > 0) {
		my ($e1) = keys %check;
		my ($node_id, $side) = split(/:/, $e1);
		my $dist = $check{$e1};
		my $node_obj = $component_node_objs{$node_id};

		# Sanity check - we must have all nodes
		defined($node_obj) || die "\nError EMain0005: could not get a node object for $node_id (edge $e1)\n\n";

		delete($check{$e1});
		$component_end_ref->{$e1} = $dist;

		my $ends_ref = $node_obj->connected($side, $type);

		# If this end is a tip then there's no reason to continue
		if(@{$ends_ref} == 0) {
			$component_end_ref = {};
			return 0;
		}

		foreach my $node_and_side (@{$ends_ref}) {
			my ($new_node_obj, $new_side) = @{$node_and_side};
			my $e2 = $new_node_obj->id . ":$new_side";

			# Have we reached this end already?
			next if(exists($component_node_objs{$new_node_obj->id}));

			# If this end is already assigned to a different cluster or we reached maximal number of steps - no good, $e_end cannot be an end  
			if(defined($new_node_obj->component($new_side)) || ($dist == $max_adjacent_component_length)) {
				$component_end_ref = {};
				return 0;
			}
			
			$component_end_ref->{$e2} = $dist+1;
			$component_node_objs{$new_node_obj->id} = $new_node_obj;

			# The following will prevent reaching $b 
			my $next_end = $new_node_obj->id . ":" . (8-$new_side);
			# Sanity check: we are not supposed to have this happen
			!exists($check{$next_end}) || die "Error EMain0006: seems like a circle of connections\n\n";
			$check{$next_end} = $dist+1;
		}
	}

	# If we are here then all is well!
	return 1;
}

###############################################################################################################################################################
# Bubble search: find start and end nodes that are the only in/out nodes to their component. Only nodes of type $type are considered - rest are ignored as 
# if they don't exist. Node type will not be considered if !defined($type).
###############################################################################################################################################################
sub search_for_bubble_components {
	my ($all_nodes_ref, $all_components_ref, $type) = @_;
	my $id = scalar(keys %{$all_components_ref}) + 1;

	foreach my $start_node_id (keys %{$all_nodes_ref}) {
		my $start_node_obj = $all_nodes_ref->{$start_node_id};

		# Skip if this node is not of the right type
		next if(defined($type) && ($start_node_obj->type ne $type));

		# Here we are looking for potential starting nodes for components
		foreach my $start_side (5, 3) {
			# No need to check if this is part of an existing component
			my $e_start = "$start_node_id:$start_side";
			next if(defined($start_node_obj->component($start_side)));

			my $ends_ref = $start_node_obj->connected($start_side, $type);

			# In order to be considered this end must have at least 2 connections
			next if(@{$ends_ref} <= 1);

			# We've got a candidate. Collect all nodes within $max_adjacent_component_length steps from this node as we allow only component that are this big.
			# There's no easy way I could think of for figuring out what are the start and end points for a component. Therefore we
			# will collect all nodes that are 10 steps or less from our candidate starting end and try each one of them as an end
			my %candidate_ends = ();
			my %check = ($e_start => 0);
			while(scalar(keys %check) > 0) {
				my ($curr_e) = keys %check;
				my $dist = $check{$curr_e};
				my ($curr_node_id, $curr_side) = split(/:/, $curr_e);

				delete($check{$curr_e});
				next if($dist == $max_adjacent_component_length);

				exists( $all_nodes_ref->{$curr_node_id}) || die "\nError EMain0024: could not find node $curr_node_id\n\n";
				my $curr_ends_ref = $all_nodes_ref->{$curr_node_id}->connected($curr_side, $type);

				# Nowhere to go from this end
				next if(@{$curr_ends_ref} == 0);

				foreach my $other_end (@{$curr_ends_ref}) {
					my ($other_node_obj, $other_side) = @{$other_end};
					my $other_e = $other_node_obj->id . ":$other_side";

					# Have we reached this end already?
					next if(exists($candidate_ends{$other_e}));

					# Maybe this end is already assigned to a different cluster? 
					next if(defined($other_node_obj->component($other_side)));

					# Keep it as a potential end
					$candidate_ends{$other_e} = $dist+1;

					# Now look to continue the search from the other side of this node
					my $next_e = $other_node_obj->id . ":" . (8-$other_side);

					# The following is a bit weird as it indicates that there are circles in our data. I am not sure as for how
					# to handle this at this point
					next if($other_node_obj->component(8-$other_side));

					# We are already checking this
					next if(exists($check{$next_e}));

					# Passed all hurdles - check it
					$check{$next_e} = $dist + 1;
				}
			}

			# Now check each of the candidates from the closest to the furthest. 
			foreach my $e_end (sort {$candidate_ends{$a} <=> $candidate_ends{$b}} keys %candidate_ends) {
				my ($end_node_id, $end_side) = split(/:/, $e_end);
				my $end_node_obj = $all_nodes_ref->{$end_node_id};
			
				defined($end_node_obj) || die "\nError EMain0007: node object for $end_node_id was not found ($e_end)\n";

				# Same node cannot be an end
				next if($start_node_id == $end_node_id);
				
				# an edge with a single connection cannot be an end
				my $end_node_connected_ref = $end_node_obj->connected($end_side, $type);
				next if(scalar(@{$end_node_connected_ref}) == 1); 

				# We have a pair of edges, now go from $e_start and make sure that all paths end in to $e_end (or in ends that lead to it)
				my %component_ends1 = ();
				my %component_ends2 = ();

				next if(!all_paths_from_a_to_b($start_node_obj, $start_side, $end_node_obj, $end_side, $type, \%component_ends1) || 
					!all_paths_from_a_to_b($end_node_obj, $end_side, $start_node_obj, $start_side, $type, \%component_ends2));

				# Make sure that we see only one side of each end node's end
				!exists($component_ends1{"$start_node_id:" . (8-$start_side)}) || next; 
				!exists($component_ends1{"$end_node_id:" . (8-$end_side)}) || next; 

				# Are the two components identical?
				# They have to be because nodes are picked based on $type (nodes that are not of type $type are ignored).
				# If the caller wants to add node later (e.g. tips?) they can call all_paths_from_a_to_b with type=undef. 
				my @component_in_node_objs = ();
				my $identical = 1;
				foreach my $e (sort {$component_ends1{$a} <=> $component_ends1{$b}} keys %component_ends1) {
					if(!exists($component_ends2{$e})) {
						$identical = 0;
						last;
					}
					($e =~ /(.+):[53]/) || die "\nError EMain0008: edge $e has an unexpected format\n\n";
					push(@component_in_node_objs, $all_nodes_ref->{$1}) if(($1 ne $end_node_id) && ($1 ne $start_node_id));
				}

				next if(!$identical);

				foreach my $e (sort {$component_ends2{$a} <=> $component_ends2{$b}} keys %component_ends2) {
					if(!exists($component_ends1{$e})) {
						$identical = 0;
						last;
					}
				}
				next if(!$identical);

				# Hurray, we've got a component! keep it
				my $component_id = scalar(keys %{$all_components_ref})+1;
				my $component_obj = new Component($component_id, $start_node_obj, $start_side, $end_node_obj, $end_side, 'bubble');

				# Start/end nodes are treated separately
				$start_node_obj->component($start_side, $component_obj);
				$end_node_obj->component($end_side, $component_obj);

				# We assign a component to the start/end nodes in $node2component_ref only if a node was not assigned already. A node can be
				# assigned to 2 different components if it has one of its ends in one component and the other in another node. This is not 
				# a healthy situation but right not it does not harm anything since $node2component_ref is used only for checking whether 
				# a node is assigned to a component or not.
				foreach my $node_obj (@component_in_node_objs) {
					$component_obj->add_node($node_obj);
					$node_obj->component(5, $component_obj);
					$node_obj->component(3, $component_obj);
				} 
				# Setting strands to nodes is done separately as the assignment of strands to one component can affect the others
				# $component_obj->set_strand_to_nodes();
				$all_components_ref->{$component_obj->id} = $component_obj;
				last;
			}
		}
	}
}

###############################################################################################################################################################
# Here we look for linear (non-bubble) components
###############################################################################################################################################################
sub search_for_linear_components {
	my ($all_nodes_ref, $all_components_ref, $type) = @_;

	# Go over all nodes, try each one as a starting point for a linear components (on both ends)
	foreach my $end1_node_id (keys %{$all_nodes_ref}) {
		my $end1_node_obj = $all_nodes_ref->{$end1_node_id};

		# Skip if this node is not of the right type
		next if(defined($type) && ($end1_node_obj->type ne $type));

		foreach my $side1 (5, 3) {
			### Step 1: check whether the current end can serve as a starting point for a component ###
			# Here we are looking for potential starting nodes for components
			my $e_end1 = "$end1_node_id:$side1";			# Starting point of the component

			# No need to check in this case...
			next if(defined($end1_node_obj->component($side1)));

			my $temp_ends_ref = $end1_node_obj->connected($side1, $type);

			# In order to be considered this end must have exactly one connection
			next if(@{$temp_ends_ref} != 1);

			# Other end of $end1_node_obj must have more than one end, or if it has one end than its connection must have more than one end,
			# or it belongs to a different component. If none of these true - this node cannot be an end node for a linear component
			if(!defined($end1_node_obj->component(8-$side1))) {
				my $temp_ends_ref = $end1_node_obj->connected(8-$side1, $type);
				if(@{$temp_ends_ref} == 1) {
					my ($temp_node_obj, $temp_side) = @{$temp_ends_ref->[0]};
		                        $temp_ends_ref = $temp_node_obj->connected($temp_side, $type);
					if(@{$temp_ends_ref} == 1) {
						# No good!
						next;
					}
				} 
			}

			# So far $end1_node_obj is a good candidate for a linear path starting node since it has exactly one connection on $side1 and more/less
			# than 1 connection on the (8-$side1) side, or the connection is not 1:1. Now we check the next node.
			my ($end2_node_obj, $side2) = @{$temp_ends_ref->[0]};

			my $e_end2 = $end2_node_obj->id . ":$side2";

			# If the following happens this means that we have a problem

			(!defined($end2_node_obj->component($side2)) || ($end2_node_obj->type ne $end1_node_obj->type)) || 
				die "\nError EMain0022: reached an end ($e_end2) that belongs to component " . $end2_node_obj->component($side2)->id . 
					" from another end ($e_end1) that does not belong to any cluster and has the same type as $e_end2\n\n"; 

			# This can happen if the nodes are not of the same type (confirmed this above) 
			!defined($end2_node_obj->component($side2)) || next;

			$temp_ends_ref = $end2_node_obj->connected($side2, $type);

			# Is this a 1:1 connection?			
			next if(@{$temp_ends_ref} != 1);
			# Hurray, we have a starting pait for the path.

			### Step 2: start elongating the path ###
			# We have a starting edge. Begin and end of the path are always ($e_end1, $e_end2).
			my @nodes_in_component = ();
			my %nodes_in_component_set = ();
			while(1) {
				my $e_next = $end2_node_obj->id . ":" . (8-$side2);

				# The following means that we have a circle
				last if($e_next eq $e_end1);

				# These are the same conditions as above
				last if(defined($end2_node_obj->component(8-$side2)));

				$temp_ends_ref = $end2_node_obj->connected((8-$side2), $type);
				last if(@{$temp_ends_ref} != 1);

				my ($temp_node_obj, $temp_side) = @{$temp_ends_ref->[0]};
				$temp_ends_ref = $temp_node_obj->connected($temp_side, $type);

				last if(@{$temp_ends_ref} != 1);

				!exists($nodes_in_component_set{$end2_node_obj->id}) || die "\nError EMain0016: a loop was detected, " . $end2_node_obj->id . " appears twice; " .
											"First end is " . $end1_node_obj->id . ":$side1, " . (scalar(@nodes_in_component)+1) . 
											" nodes in ccomponent so far\n\n";
				# ok, we can add this node to the component. This means that the node of $e_end2 belongs in its entirety to the component
				push(@nodes_in_component, $end2_node_obj->id);
				$nodes_in_component_set{$end2_node_obj->id} = 1;
				($end2_node_obj, $side2) = ($temp_node_obj, $temp_side);
			}

			# Hurray, we've got a component! keep it
			my $component_id = scalar(keys %{$all_components_ref})+1;
			my $component_obj = new Component($component_id, $end1_node_obj, $side1, $end2_node_obj, $side2, 'linear');

			# Start/end nodes are treated separately
			$end1_node_obj->component($side1, $component_obj);
			$end2_node_obj->component($side2, $component_obj);

			foreach my $node_id (@nodes_in_component) {
				$component_obj->add_node($all_nodes_ref->{$node_id}, "white");
				$all_nodes_ref->{$node_id}->component(5, $component_obj);
				$all_nodes_ref->{$node_id}->component(3, $component_obj);

			} 
			# Setting strands to nodes is done separately as the assignment of strands to one component can affect the others
			# $component_obj->set_strand_to_nodes();
			$all_components_ref->{$component_id} = $component_obj;		
		}
	}
}

###############################################################################################################################################################
sub assign_nodes_to_components {
	my ($all_nodes_ref, $all_edges_ref, $all_components_ref, $type) = @_;

	my %check = ();
	my %incoming_tip_ends = ();
	my %tip_assignments = ();

	# Start by assigning all $type nodes that are connected to non-tip nodes
	foreach my $edge_obj (@{$all_edges_ref}) {
		# We need to check if this edge changes that assignment of a tip node
		my ($node1_obj, $side1) = @{$edge_obj->end1};
		my ($node2_obj, $side2) = @{$edge_obj->end2};
		my ($nontip_node_obj, $nontip_side, $tip_node, $tip_side) = (undef, undef, undef, undef);

		if(($node1_obj->type ne $type) && ($node2_obj->type eq $type)) {
			($nontip_node_obj, $nontip_side, $tip_node, $tip_side) = ($node1_obj, $side1, $node2_obj, $side2);
		}
		elsif(($node2_obj->type ne $type) && ($node1_obj->type eq $type)) {
			($nontip_node_obj, $nontip_side, $tip_node, $tip_side) = ($node2_obj, $side2, $node1_obj, $side1);
		}
		else {
			next;
		}

		my $component_id = 'None';
		if(defined($nontip_node_obj->component(5))) {
			(!defined($nontip_node_obj->component(3)) || ($nontip_node_obj->component(5)->id == $nontip_node_obj->component(3)->id)) ||
				die "EMain0023: different component assignments to node " . $nontip_node_obj->id . " on both ends, expecting just one\n\n";
			$component_id = $nontip_node_obj->component(5)->id;
		}
		elsif(defined($nontip_node_obj->component(5))) {
			$component_id = $nontip_node_obj->component(3)->id;
		}

		if(exists($tip_assignments{$tip_node->id})) {
			!exists($incoming_tip_ends{$tip_node->id . ":" . (8-$tip_side)}) || die "\nError EMain0009: a node is assigned from two sides\n\n";
			if($tip_assignments{$tip_node->id} != $component_id) {
				$tip_assignments{$tip_node->id} = $incoming_tip_ends{$tip_node->id . ":$tip_side"} = "Multiple";
			}
		}
		else {
			$tip_assignments{$tip_node->id} = $incoming_tip_ends{$tip_node->id . ":$tip_side"} = $component_id;
			$check{$tip_node->id} = (8-$tip_side);
		}
	}

	# Now keep on passing the assignment until %check is empty
	while(scalar(keys %check) > 0) {
		my ($out_node_id) = keys %check;
		my $out_node_obj = $all_nodes_ref->{$out_node_id};
		my $out_side = $check{$out_node_id};
		my $component_assignment = $tip_assignments{$out_node_id};
		delete($check{$out_node_id});
		my $out_ends_ref = $out_node_obj->connected($out_side);

		foreach my $end_ref (@{$out_ends_ref}) {
			my ($in_node_obj, $in_side) = @{$end_ref};
			($in_node_obj->type eq $type) || die "\nError EMain0010: a non-$type accessed from a $type\n\n";
			if(exists($tip_assignments{$in_node_obj->id})) {
				exists($incoming_tip_ends{$in_node_obj->id . ":$in_side"}) || die "\nError EMain0011: a node is assigned from two sides\n\n";
				if($tip_assignments{$in_node_obj->id} ne $component_assignment) {
					$tip_assignments{$in_node_obj->id} = $incoming_tip_ends{$in_node_obj->id . ":$in_side"} = "Multiple";
					$check{$in_node_obj->id} = $in_side;
				}
			}
			else {
				$tip_assignments{$in_node_obj->id} = $incoming_tip_ends{$in_node_obj->id . ":$in_side"} = $component_assignment;
				$check{$in_node_obj->id} = $in_side;
			} 
		}
	}

	# Finally - add all nodes with new assignments
	foreach my $node_id (keys %tip_assignments) {
		next if(($tip_assignments{$node_id} eq 'Multiple') || ($tip_assignments{$node_id} eq 'None'));
		$all_components_ref-{$tip_assignments{$node_id}}->add_node($all_nodes_ref->{$node_id});
	}
}

###############################################################################################################################################################
sub get_other_end_of_node {
	my $e = shift;
	my $s = chop($e);
	(($s == 5) || ($s == 3)) || die "\nError EMain0012: unexpected side $s for $e$s\n\n";
	$e .= (8-$s);

	return $e;
}

###############################################################################################################################################################
# This subroutine will set strand to all nodes without strand assignment yet
sub set_strands_to_nodes {
	my $all_node_ref = shift(@_);
	
	foreach my $node_obj (values %{$all_node_ref}) {
		my $side = undef;
		next if(defined($node_obj->strand));
		my $node_strand = undef;
		foreach my $node_side (3, 5) {
			my $edges_ref = $node_obj->edges($node_side);
			foreach my $edge_obj (@{$edges_ref}) {
				my ($other_node_obj, $other_node_side) = $edge_obj->other_end($node_obj->id, $node_side);
				if(defined($other_node_obj->strand)) {
					if($other_node_obj->strand eq 'N/A') {
						$node_strand = 'N/A';
						last;
					}

					my $strand = '+';
					$strand = '-' if((($node_side == $other_node_side) && ($other_node_obj->strand eq '+')) ||
							 (($node_side != $other_node_side) && ($other_node_obj->strand eq '-')));
					if(!defined($node_side)) {
						$node_strand = $strand;
					}
					elsif($node_strand ne $strand) {
						$node_strand = 'N/A';
						last;
					}
				}
			}
			last if($node_strand eq 'N/A');
		}
		$node_strand = '+' if(!defined($node_strand));
		$node_obj->strand($node_strand);
	}
}

###############################################################################################################################################################
sub create_meta_components {
	my ($all_components_ref, $all_metacomponents_ref) = @_;

	# First - collect connections between components
	my %node_objects = ();
	my %end2component = ();
	while (my ($component_id, $component_obj) = each %{$all_components_ref}) {
		my ($node_obj1, $side1, $node_obj2, $side2) = ($component_obj->node1, $component_obj->side1, $component_obj->node2, $component_obj->side2);
		$node_objects{$node_obj1->id} = $node_obj1;
		$node_objects{$node_obj2->id} = $node_obj2;
		if(exists($end2component{$node_obj1->id})) {
			!exists($end2component{$node_obj1->id}{$side1}) || die "\nError EMain0018: " . $node_obj1->id . ":$side1 is assigned to two different components\n\n";
		}
		$end2component{$node_obj1->id}{$side1} = $component_obj;
		if(exists($end2component{$node_obj2->id})) {
			!exists($end2component{$node_obj2->id}{$side2}) || die "\nError EMain0019: " . $node_obj2->id . ":$side2 is assigned to two different components\n\n";
		}
		$end2component{$node_obj2->id}{$side2} = $component_obj;
	}

	# Second - look for meta-components
	my @node_ids = keys %end2component;
	foreach my $node1_id (@node_ids) {
		# Have we used this end already?
		next if(!exists($end2component{$node1_id}));

		my $side1 = undef;
		if(exists($end2component{$node1_id}{5}) && defined($end2component{$node1_id}{5}) && (!exists($end2component{$node1_id}{3}) || !defined($end2component{$node1_id}{3}))) {
			$side1 = 5;
		}
		elsif(exists($end2component{$node1_id}{3}) && defined($end2component{$node1_id}{3}) && (!exists($end2component{$node1_id}{5}) || !defined($end2component{$node1_id}{5}))) {
			$side1 = 3;
		}
		next if(!defined($side1));

		$end2component{$node1_id}{$side1}->set_strand_to_nodes();
		my $metacomponent_obj = extract_metacomponent_path($node1_id, $side1, \%end2component, \%node_objects, scalar(keys %{$all_metacomponents_ref})+1000001);

		if(defined($metacomponent_obj)) {
			$metacomponent_obj->type("meta-circular");
			$all_metacomponents_ref->{$metacomponent_obj->id} = $metacomponent_obj if(defined($metacomponent_obj));
		}
	}

	# Finally - check if we have circular meta
	@node_ids = keys %end2component;
	foreach my $node1_id (@node_ids) {
		next if(!exists($end2component{$node1_id}));
		if(exists($end2component{$node1_id}{5}) && defined($end2component{$node1_id}{5}) && exists($end2component{$node1_id}{3}) && defined($end2component{$node1_id}{3})) {
			my $metacomponent_obj = extract_metacomponent_path($node1_id, 3, \%end2component, \%node_objects, scalar(keys %{$all_metacomponents_ref})+1);
			defined($metacomponent_obj) || die "\nError EMain0013: could not reconstruct circular meta-component\n\n";
			$all_metacomponents_ref->{$metacomponent_obj->id} = $metacomponent_obj;
		}
	}
}

###############################################################################################################################################################
sub extract_metacomponent_path {
	my ($node1_id, $side1, $end2component_ref, $node_objects_ref, $metacomponent_id) = @_;
	my $c_obj = $end2component_ref->{$node1_id}{$side1};
	my ($node2_id, $side2) = ("$node1_id:$side1" eq $c_obj->end2)? ($c_obj->node1->id, $c_obj->side1) : ($c_obj->node2->id, $c_obj->side2);
	my @node_ids = ($node1_id, $node2_id);
	my @component_objs = ($c_obj);

	while(exists($end2component_ref->{$node2_id}{8-$side2}) && defined($end2component_ref->{$node2_id}{8-$side2})) {
		# The following can happen in cases of circular meta-components
		if($node2_id eq $node1_id) {
			($side1 == (8-$side2)) || die "\nError EMain0014: completed a loop but from the wrong direction\n\n";
			last;
		}

		$c_obj = $end2component_ref->{$node2_id}{8-$side2};
		($node2_id, $side2) = ("$node2_id:" . (8-$side2) eq $c_obj->end2)? ($c_obj->node1->id, $c_obj->side1) : ($c_obj->node2->id, $c_obj->side2);

		push(@component_objs, $c_obj);
		push(@node_ids, $node2_id);
		$c_obj->set_strand_to_nodes();
	}
	if(@component_objs > 1) {
		my $metacomponent_obj = new Metacomponent($metacomponent_id, $node_objects_ref->{$node1_id}, $side1, $node_objects_ref->{$node2_id}, $side2, 'meta');
		foreach my $component_obj (@component_objs) {
			$metacomponent_obj->add_node($component_obj);				
		}
		foreach my $node_id (@node_ids) {
			delete($end2component_ref->{$node_id});
		}
		return $metacomponent_obj;
	}
	return undef;
}
