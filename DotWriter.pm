#!/usr/bin/perl

use strict;
use Node;
use Edge;
use Component;
use Writer;

package DotWriter;

our @ISA = qw(Writer);

#######################################################################################################################
sub new {
        my $class = shift;

	my $self = $class->SUPER::new(@_);

        bless $self, $class;

        return $self;
}

#######################################################################################################################
sub write {
        my ( $self ) = shift;
	my ( $all_metacomponents_ref, $all_components_ref, $all_nodes_ref, $all_edges_ref, $outfile_name ) = @_;

	my %component_id2metacomponent_id = ();
	my %multinode_id2size = ();		# Will contain both components and metacomponents

	$self->SUPER::get_multinode_info($all_metacomponents_ref, $all_components_ref, \%component_id2metacomponent_id, \%multinode_id2size);

	open(OUTG, ">$outfile_name") || die "\nError EDotWriter0001: Cannot write to $outfile_name\n\n";

	print OUTG "digraph G {\n";
	# Write all components/metacomponents
	my %written_node_ids = ();
	my %written_edge_ids = ();

	foreach my $id (reverse sort {$multinode_id2size{$a} <=> $multinode_id2size{$b}} keys %multinode_id2size) {
		if(exists($all_metacomponents_ref->{$id})) {
			$self->write_metacomponent($all_metacomponents_ref->{$id}, \*OUTG, "\t", $multinode_id2size{$id}, \%written_node_ids, \%written_edge_ids);
		}
		else {
			$self->write_component($all_components_ref->{$id}, \*OUTG, "\t", \%written_node_ids, \%written_edge_ids);
		}
	}

	# Now write all component-less nodes...
	while (my ($node_id, $node_obj) = each %{$all_nodes_ref}) {
		$self->write_node($node_obj, \*OUTG, "\t", 'white', \%written_node_ids) if(!exists($written_node_ids{$node_id}));
	}

	# ... and edges
	foreach my $edge_obj (@{$all_edges_ref}) {
		$self->write_edge($edge_obj, \*OUTG, "\t", \%written_edge_ids) if(!exists($written_edge_ids{$edge_obj->id}));
	}
	print OUTG "}\n";
}

#######################################################################################################################
sub write_node {
	my ( $self, $node_obj, $OUT_ref, $tabs, $fillcolor, $written_node_ids_ref ) = @_;
	print $OUT_ref $tabs, $node_obj->id, " [label=\"", $node_obj->id, "\\nSize: ", $node_obj->size, "\",shape=box, fillcolor=$fillcolor, style=filled, color=black, penwidth=2]\n";
	$written_node_ids_ref->{$node_obj->id} = 1;
}

#######################################################################################################################
sub write_edge {
	my ( $self, $edge_obj, $OUT_ref, $tabs, $written_edge_ids_ref ) = @_;

	my ($node1_obj, $side1) = @{$edge_obj->end1};
	my ($node2_obj, $side2) = @{$edge_obj->end2};
	my $h1 = ($side1 == 5)? "normal" : "inv";
	my $h2 = ($side2 == 5)? "normal" : "inv";

	print $OUT_ref 	$tabs, $node1_obj->id, " -> ", $node2_obj->id, " [taillabel=\"$side1\", headlabel=\"$side2\", dir=both, label=\"",
			$edge_obj->nsequences, "\", arrowhead=$h2, arrowtail=$h1, color=black, penwidth=2]\n";
	$written_edge_ids_ref->{$edge_obj->id} = 1;
}

#######################################################################################################################
sub write_component {
        my ( $self, $component_obj, $OUT_ref, $tabs, $written_node_ids_ref, $written_edge_ids_ref ) = @_;

	print $OUT_ref $tabs, "subgraph cluster_", $component_obj->id, " {\n";
	print $OUT_ref $tabs, "\tfontcolor=\"black\" fontsize=16.0;\n";
	print $OUT_ref $tabs, "\tstyle=filled;\n";
	print $OUT_ref $tabs, "\tcolor=cyan;\n";
	print $OUT_ref $tabs, "\tlabel=\"", $component_obj->id, " (", $component_obj->size, " nodes, ", $component_obj->type, ")\"\n";

	### Write nodes ###
	# Write node1 if it was not written yet
	if(!exists($written_node_ids_ref->{$component_obj->node1->id})) {
		$self->write_node($component_obj->node1, $OUT_ref, "$tabs\t", 'red', $written_node_ids_ref);
	}

	for(my $i=1; $component_obj->in_node($i); $i++) {
		$self->write_node($component_obj->in_node($i), $OUT_ref, "$tabs\t", 'yellow', $written_node_ids_ref);
	}

	# Write node2 if it was not written yet
	if(!exists($written_node_ids_ref->{$component_obj->node2->id})) {
		$self->write_node($component_obj->node2, $OUT_ref, "$tabs\t", 'red', $written_node_ids_ref);
	}

	print $OUT_ref "\n";

	### Write edges ###
	for(my $i=1; $component_obj->in_node($i); $i++) {
		for my $s (5, 3) {
			foreach my $edge_obj (@{$component_obj->in_node($i)->edges($s)}) {
				$self->write_edge($edge_obj, $OUT_ref, "$tabs\t", $written_edge_ids_ref) if(!exists($written_edge_ids_ref->{$edge_obj->id}));
			}
		}
	}

	print $OUT_ref $tabs, "}\n\n";
}

#######################################################################################################################
sub write_metacomponent {
        my ( $self, $metacomponent_obj, $OUT_ref, $tabs, $metacomponent_size, $written_node_ids_ref, $written_edge_ids_ref ) = @_;

	print $OUT_ref "$tabs", "subgraph cluster_", $metacomponent_obj->id, " {\n";
	print $OUT_ref "$tabs\tfontcolor=\"black\" fontsize=16.0;\n";
	print $OUT_ref "$tabs\tstyle=filled;\n";
	print $OUT_ref "$tabs\tcolor=yellow;\n";
	print $OUT_ref "$tabs\tlabel=\"", $metacomponent_obj->id, " ($metacomponent_size nodes)\"\n";

	my $i = 1;
	while (my $component_obj = $metacomponent_obj->in_node($i++)) {
		if(!exists($written_node_ids_ref->{$component_obj->node1->id})) {
			$self->write_node($component_obj->node1, $OUT_ref, "$tabs\t", 'red', $written_node_ids_ref);
		}
		if(!exists($written_node_ids_ref->{$component_obj->node2->id})) {
			$self->write_node($component_obj->node2, $OUT_ref, "$tabs\t", 'red', $written_node_ids_ref);
		}		
	}

	$i = 1;
	while (my $component_obj = $metacomponent_obj->in_node($i++)) {
		$self->write_component($component_obj, $OUT_ref, "$tabs\t", $written_node_ids_ref, $written_edge_ids_ref);
	}
	print $OUT_ref "$tabs}\n";
}

1;
