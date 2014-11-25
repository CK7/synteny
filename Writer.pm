#!/usr/bin/perl

use strict;
use Node;
use Edge;
use Component;

package Writer;

#######################################################################################################################
sub new {
	my $class = shift;
	my $self = {
	};

	bless $self, $class;

	return $self;
}

#######################################################################################################################
sub get_multinode_info {
        my ( $self ) = shift;
	my ( $all_metacomponents_ref, $all_components_ref, $component_id2metacomponent_id_ref, $multinode_id2size_ref ) = @_;

	# Figure out each component/meta-component size so that we can write them in decreasing size
	while(my ($metacomponent_id, $metacomponent_obj) = each %{$all_metacomponents_ref}) {
		$multinode_id2size_ref->{$metacomponent_id} = ($metacomponent_obj->node1 eq $metacomponent_obj->node2)? 0 : 1;
		my $i = 1;
		while (my $component_obj = $metacomponent_obj->in_node($i++)) {
			$component_id2metacomponent_id_ref->{$component_obj->id} = $metacomponent_id;
			$multinode_id2size_ref->{$metacomponent_id} += $component_obj->size() - 1;
		}
	}
	while(my ($component_id, $component_obj) = each %{$all_components_ref}) {
		next if(exists($component_id2metacomponent_id_ref->{$component_id}));
		$multinode_id2size_ref->{$component_id} = $component_obj->size;
	}
}

1;
