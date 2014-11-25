#!/usr/bin/perl

use strict;
use Node;

package Component;

#######################################################################################################################
sub new {
	my $class = shift;
	my $self = {
		_id => shift,
		_node_obj1 => shift,
		_side1 => shift,
		_node_obj2 => shift,
		_side2 => shift,
		_type  => shift,
                _internal_nodes2index => {},
		_index2internal_nodes => {},
        };

        bless $self, $class;

        return $self;
}

#######################################################################################################################
sub id {
	my ( $self ) = shift;

	return $self->{_id};
}

#######################################################################################################################
sub reverse_strand_to_nodes {
	my ( $self ) = shift;

	($self->node1->strand eq '+')? $self->node1->strand('-') : $self->node1->strand('+');
	($self->node2->strand eq '+')? $self->node2->strand('-') : $self->node2->strand('+');

	foreach my $in_node_obj (values %{$self->{_index2internal_nodes}}) {
		($in_node_obj->strand eq '+')? $in_node_obj->strand('-') : $in_node_obj->strand('+');
	}
}

#######################################################################################################################
sub set_strand_to_nodes {
	my ( $self ) = shift;

	# First make assignment to node1 and node2. One (or both) of them may already have an assignment
	my $strand1 = ($self->side1 == 3)? "+" : "-";
	my $strand2 = ($self->side2 == 5)? "+" : "-";

	if((defined($self->node1->strand) && ($self->node1->strand ne $strand1)) || (defined($self->node2->strand) && ($self->node2->strand ne $strand2))) {
		$strand1 = ($strand1 eq '+')? '-' : '+';
		$strand2 = ($strand2 eq '+')? '-' : '+';
	}

	if(defined($self->node1->strand) && ($self->node1->strand ne $strand1)) {
		die "Error EComponent0001: node " . $self->node1->id . " already has a strand assigned to it but new assignment contradicts old assignment\n\n";
	}
	else {
		$self->node1->strand($strand1);
	}

	if(defined($self->node2->strand) && ($self->node2->strand ne $strand2)) {
		die "Error EComponent0002: node " . $self->node2->id . " already has a strand assigned to it but new assignment contradicts old assignment\n\n";
	}
	else {
		$self->node2->strand($strand2);
	}

	# Now make assignments for all internal nodes
	my $change = 1;
	while($change) {
		$change = 0;
		foreach my $in_node_obj (values %{$self->{_index2internal_nodes}}) {
			next if(defined($in_node_obj->strand));
			foreach my $in_node_side (3, 5) {
				my $edges_ref = $in_node_obj->edges($in_node_side);
				foreach my $edge_obj (@{$edges_ref}) {
					my ($other_node_obj, $other_node_side) = $edge_obj->other_end($in_node_obj->id, $in_node_side);
					if(defined($other_node_obj->strand)) {
						my $strand = '+';
						$strand = '-' if((($in_node_side == $other_node_side) && $other_node_obj->strand eq '+') || 
								 (($in_node_side != $other_node_side) && $other_node_obj->strand eq '-'));
						$in_node_obj->strand($strand);
						$change = 1;
						last;
					}
				}
			last if(defined($in_node_obj->strand));
			}
		}
	}

	# Now make sure that all strand assignments are consistent

	foreach my $in_node_obj (values %{$self->{_index2internal_nodes}}) {
		defined($in_node_obj->strand) || die "\nError EComponent0003: did not assign strand to " . $in_node_obj->id . "\n\n";
		foreach my $in_node_side (3, 5) {
			my $edges_ref = $in_node_obj->edges($in_node_side);
			foreach my $edge_obj (@{$edges_ref}) {
				my ($other_node_obj, $other_node_side) = $edge_obj->other_end($in_node_obj->id, $in_node_side);
				exists($self->{_index2internal_nodes}{$other_node_obj->id}) || next;

				if((($in_node_side != $other_node_side) && ($in_node_obj->strand ne $other_node_obj->strand)) ||
				   (($in_node_side == $other_node_side) && ($in_node_obj->strand eq $other_node_obj->strand))) {
					print STDERR "Error EComponent0004: nodes " . $in_node_obj->id . " (" . $in_node_obj->strand. ") and " . $other_node_obj->id . " (" . $other_node_obj->strand. ") have conflicting strands\n";
				}
			}
		}
	}
}

#######################################################################################################################
sub size {
	my ( $self ) = shift;
	my $num_nodes = scalar(keys %{$self->{_internal_nodes2index}}); 
	$num_nodes += ($self->{_node_obj1} eq $self->{_node_obj2})? 1 : 2;
	
	return $num_nodes;
}

#######################################################################################################################
sub add_node {
	my ( $self ) = shift;
	my $node_obj = shift;

	if(($node_obj ne $self->{_node_obj1}) && ($node_obj ne $self->{_node_obj2}) && (!exists($self->{_internal_nodes2index}{$node_obj}))) {
		$self->{_internal_nodes2index}{$node_obj} = scalar(keys %{$self->{_internal_nodes2index}})+1;
		$self->{_index2internal_nodes}{scalar(keys %{$self->{_index2internal_nodes}})+1} = $node_obj;
	}
}

#######################################################################################################################
sub type {
	my ( $self ) = shift;

	$self->{_type} = $_[0] if(@_ == 1);

	return $self->{_type}; 
}

#######################################################################################################################
# This subroutine will return any of the non-end (node1, node2) nodes
#######################################################################################################################
sub in_node {
	my ( $self ) = shift;
	my $index = shift;
	return (($index > 0) && ($index <= scalar(keys %{$self->{_index2internal_nodes}})))? $self->{_index2internal_nodes}{$index} : undef;
}

#######################################################################################################################
sub node1 {
	my ( $self ) = shift;

	return $self->{_node_obj1}; 
}

#######################################################################################################################
sub node2 {
	my ( $self ) = shift;

	return $self->{_node_obj2}; 
}

#######################################################################################################################
sub side1 {
	my ( $self ) = shift;

	return $self->{_side1}; 
}

#######################################################################################################################
sub side2 {
	my ( $self ) = shift;

	return $self->{_side2}; 
}

#######################################################################################################################
sub end1 {
	my ( $self ) = shift;

	return $self->{_node_obj1}->id . ":" . $self->{_side1}; 
}

#######################################################################################################################
sub end2 {
	my ( $self ) = shift;

	return $self->{_node_obj2}->id . ":" . $self->{_side2}; 
}

1;
