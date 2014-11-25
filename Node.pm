#!/usr/bin/perl

use strict;
use Gene;
use Sequence;

package Node;

#######################################################################################################################
sub new {
	my $class = shift;

	($#_ == 1) || die "\nError ENode0001: expected 2 arguments in constructor but got " . scalar(@_) . " instead\n\n";

	my $self = {
		_id => shift,
		_type => shift,
		"_component:5" => undef,
		"_component:3" => undef,
		_5 => [],
		_3 => [],
		_genes => [],
		_strand => undef,
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
sub strand {
	my ( $self, $strand ) = @_;
	
	$self->{_strand} = $strand if(defined($strand));

	return $self->{_strand};
}

#######################################################################################################################
sub add_gene {
	my ( $self ) = shift;

	push(@{$self->{_genes}}, $_[0]);
}

#######################################################################################################################
sub remove_gene {
	my ( $self ) = shift;
	my $gene_obj = shift;

	foreach my $i (0 .. $#{$self->{_genes}}) {
		if($self->{_genes}->[$i]->id eq $gene_obj->id) {
			splice(@{$self->{_genes}}, $i, 1);
			last;
		}
	}
}

#######################################################################################################################
sub genes {
	my ( $self ) = shift;

	push(@{$self->{_genes}}, $_[0]) if(@_ == 1);

	return $self->{_genes};
}

#######################################################################################################################
sub sequences {
	my ( $self ) = shift;

	my @seqs = ();

	foreach my $gene_obj (@{$self->{_genes}}) {
		push(@seqs, $gene_obj->seq);
	}

	return \@seqs;
}

#######################################################################################################################
sub edges {
	my ( $self ) = shift;
	my $side = shift;

	return $self->{"_$side"};
}

#######################################################################################################################
sub connected {
	my ( $self ) = shift;
	my ($side, $type) = @_;
	
	my @ends = ();
	foreach my $edge_obj (@{$self->{"_$side"}}) {
		my ($other_node_obj, $other_side) = $edge_obj->other_end($self->id, $side);
		push(@ends, [$other_node_obj, $other_side]) if(!defined($type) || ($type eq $other_node_obj->type));
	}

	return \@ends;
}

#######################################################################################################################
sub add_edge {
	my ( $self ) = shift;
	my ( $side, $edge_obj) = @_;

	push(@{$self->{"_$side"}}, $edge_obj);
}

#######################################################################################################################
sub size {
	my ( $self ) = shift;

	return scalar(@{$self->{_genes}});
}

#######################################################################################################################
sub component {
	my ( $self ) = shift;
	my $side = shift;
	if(@_ == 1) {
		my $comp_obj = shift;
		$self->{"_component:$side"} = $comp_obj;
	}

	return $self->{"_component:$side"};
}

#######################################################################################################################
sub type {
	my ( $self ) = shift;

	$self->{_type} = shift if(@_ == 1);

	return $self->{_type};
}

1;
