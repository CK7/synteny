#!/usr/bin/perl

use strict;
use Gene;
use Sequence;
use Node;

package Edge;

#######################################################################################################################
sub new {
        my $class = shift;

	(scalar(@_) == 4) || die "\nError EEdge0001: expected 4 arguments, received " . scalar(@_) . "\n\n";
        my $self = {
                _node1_obj => shift,
                _side1 => shift,
                _node2_obj  => shift,
                _side2 => shift,
                _sequences => (),
        };

	($self->{_node1_obj}, $self->{_side1}, $self->{_node2_obj}, $self->{_side2}) = ($self->{_node2_obj}, $self->{_side2}, $self->{_node1_obj}, $self->{_side1})
		if(($self->{_node1_obj}->id gt $self->{_node2_obj}->id) || (($self->{_node1_obj}->id eq $self->{_node2_obj}->id) && ($self->{_side1} < $self->{_side2})));
        bless $self, $class;

        return $self;
}

#######################################################################################################################
sub id {
        my ( $self ) = shift;

	return $self->{_node1_obj}->id . ":" . $self->{_side1} . " -> " . $self->{_node2_obj}->id . ":" . $self->{_side2};
}

#######################################################################################################################
sub add_sequence {
        my ( $self ) = shift;
	my $seq_obj = shift;

        push(@{$self->{_sequences}}, $seq_obj);
}

#######################################################################################################################
sub nsequences {
        my ( $self ) = shift;

        return scalar(@{$self->{_sequences}});
}

#######################################################################################################################
sub sequences {
        my ( $self ) = shift;

        return $self->{_sequences};
}

#######################################################################################################################
sub end1 {
        my ( $self ) = shift;

	return [$self->{_node1_obj}, $self->{_side1}];
}

#######################################################################################################################
sub end2 {
        my ( $self ) = shift;

	return [$self->{_node2_obj}, $self->{_side2}];
}

#######################################################################################################################
sub replace_end {
        my ( $self ) = shift;
	my ($old_node_obj, $old_s, $new_node_obj, $new_s) = @_;

	if(($self->{_node1_obj}->id eq $old_node_obj->id) && ($self->{_side1} == $old_s)) {
		($self->{_node1_obj}, $self->{_side1}) = ($new_node_obj, $new_s);
	}
	elsif(($self->{_node2_obj}->id eq $old_node_obj->id) && ($self->{_side2} == $old_s)) {
		($self->{_node2_obj}, $self->{_side2}) = ($new_node_obj, $new_s);
	}
	else {
		die "\nError EEdge0002: provided old edge does not match anything\n\n";
	}
}

#######################################################################################################################
sub other_end {
        my ( $self ) = shift;
	my $other_node_id = shift;
	my $other_side = shift;

	if(!defined($other_side)) {
		($other_node_id =~ /^(.+):([35])/) || die "\nError EEdge0003: unexpected input for Edge::other_end: $other_node_id\n\n";
		($other_node_id, $other_side) = ($1, $2);
	}

	if(($self->{_node1_obj}->id eq $other_node_id) && ($self->{_side1} == $other_side)) {
		return ($self->{_node2_obj}, $self->{_side2});
	}
	elsif(($self->{_node2_obj}->id eq $other_node_id) && ($self->{_side2} == $other_side)) {
		return ($self->{_node1_obj}, $self->{_side1});
	}

	die "\nError EEdge0004: $other_node_id:$other_side is not an end for the current edge\n\n";
}

1;
