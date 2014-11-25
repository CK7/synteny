#!/usr/bin/perl

use strict;
use Node;
use Gene;

package Sequence;

#######################################################################################################################
sub new {
        my $class = shift;
        my $self = {
                _id => shift,
		_genes => [],
		_gene_true_index => [],
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
sub genes {
        my ( $self ) = shift;

        return $self->{_genes};
}

#######################################################################################################################
sub ngenes {
        my ( $self ) = shift;

        return scalar(@{$self->{_genes}});
}

#######################################################################################################################
sub nodes {
        my ( $self ) = shift;

	my @nodes;

	foreach my $gene_obj (@{$self->{_genes}}) {
		push(@nodes, $gene_obj->node);
	}

        return \@nodes;
}

#######################################################################################################################
sub add_gene {
        my ( $self ) = shift;
	my $gene_obj = shift;
	my $gene_true_index = shift;

	defined($gene_true_index) || die "\nError ESequence0001: gene true index was not provided\n\n";

	if((@{$self->{_genes}} == 0) || ($gene_obj->start < $self->{_genes}->[0]->start)) {
		unshift(@{$self->{_genes}}, $gene_obj);
		unshift(@{$self->{_gene_true_index}}, $gene_true_index); 
		return;
	}
	elsif($gene_obj->start >= $self->{_genes}->[$#{$self->{_genes}}]->start) {
		push(@{$self->{_genes}}, $gene_obj);
		push(@{$self->{_gene_true_index}}, $gene_true_index);
		return;
	}
	else {
		foreach my $i (1 .. $#{$self->{_genes}}) {
			if(($gene_obj->start >= $self->{_genes}->[$i-1]->start) && ($gene_obj->start < $self->{_genes}->[$i]->start)) {
				splice(@{$self->{_genes}}, $i, 0, $gene_obj);
				splice(@{$self->{_gene_true_index}}, $i, 0, $gene_true_index);
				return;
			} 
		}
	}

	die "\nError ESequence:0002: was not able to insert a gene\n\n";
}

#######################################################################################################################
sub remove_gene {
        my ( $self ) = shift;
	my $gene_obj = shift;

	for(my $i=0; $i<=$#{$self->{_genes}}; $i++) {
		if($self->{_genes}->[$i]->id eq $gene_obj->id) {
			splice(@{$self->{_genes}}, $i, 1);
			splice(@{$self->{_gene_true_index}}, $i, 1);
			return;
		}
	}
	die "\nError ESequence:0003: tried to remove gene " . $gene_obj->id . " that does not seem to be part of sequence " . $self->id . "\n\n"; 
}

#######################################################################################################################
# Will return the number of genes between $gene_obj1 and $gene_obj2 INCLUDING genes that were not considered for this 
# sequence (e.g. because they appear in only a single copy), and also the distance in bps between the end of the first 
# gene and the beginning of the second one.
#######################################################################################################################
sub gene_distance {
        my ( $self ) = shift;
	my ($gene_obj1, $gene_obj2) = @_;
	my ($pos1, $pos2) = (undef, undef);

	for(my $i=0; $i<=$#{$self->{_genes}}; $i++) {
		$pos1 = $self->{_gene_true_index}->[$i] if($self->{_genes}->[$i]->id eq $gene_obj1->id);
		$pos2 = $self->{_gene_true_index}->[$i] if($self->{_genes}->[$i]->id eq $gene_obj2->id);
	}

	(defined($pos1) && defined($pos2)) || return (undef, undef);

	if($gene_obj1->start < $gene_obj2->start) {
		return ($pos2-$pos1, $gene_obj2->start-$gene_obj1->end);
	}
	return ($pos1-$pos2, $gene_obj1->start-$gene_obj2->end);
}
1;
