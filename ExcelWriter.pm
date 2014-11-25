#!/usr/bin/perl

use strict;
use Node;
use Edge;
use Component;
use Writer;

package ExcelWriter;

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
	my ( $all_metacomponents_ref, $all_components_ref, $all_nodes_ref, $all_equences_ref, $outfile_name ) = @_;

        my %component_id2metacomponent_id = ();
        my %multinode_id2size = ();             # Will contain both components and metacomponents

        $self->SUPER::get_multinode_info($all_metacomponents_ref, $all_components_ref, \%component_id2metacomponent_id, \%multinode_id2size);

	open(OUTG, ">$outfile_name") || die "\nError EExcelWriter0001: Cannot write to $outfile_name\n\n";

	# Write all components/metacomponents
	my %written_node_ids = ();
	my %written_edge_ids = ();

	foreach my $id (reverse sort {$multinode_id2size{$a} <=> $multinode_id2size{$b}} keys %multinode_id2size) {
		if(exists($all_metacomponents_ref->{$id})) {
			my ($out_line, $total_nbps, $total_ngenes, $total_nbackbone) = $self->write_metacomponent($all_metacomponents_ref->{$id}, \%written_node_ids);
			print OUTG "### METACOMPONENT_$id (", int($total_nbps), " bps, $total_ngenes genes, $total_nbackbone backbone genes) ###\n$out_line";
		}
		else {
			my ($out_line, $total_nbps, $total_ngenes, $total_nbackbone) = $self->write_component($all_components_ref->{$id}, 1, 1, \%written_node_ids);
			print OUTG "### COMPONENT_$id (", int($total_nbps), " bps, $total_ngenes genes, $total_nbackbone backbone genes) ###\n";
			print OUTG $out_line;
		}
	}

	# Now write all component-less nodes...
	print OUTG "### NO COMPONENT ###\n";
	while (my ($node_id, $node_obj) = each %{$all_nodes_ref}) {
		$self->write_node($node_obj, \*OUTG, \%written_node_ids) if(!exists($written_node_ids{$node_id}));
	}

	close(OUTG);
}

#######################################################################################################################
sub write_node {
	my $self = shift;
	my ($node_obj, $OUTG_ref, $written_node_ids_ref) = @_;

	# Compute average length
	my $avg_length = 0;
	foreach my $gene_obj (@{$node_obj->genes}) {
		$avg_length += $gene_obj->length;
	}
	$avg_length /= $node_obj->size;

	print $OUTG_ref $node_obj->id, "\tUNCLUSTERED\t+\t", int(10*$avg_length)/10, "\t", $node_obj->size, "\n";
}

#######################################################################################################################
sub write_metacomponent {
	my $self = shift;
	my ($metacomponent_obj, $written_node_ids_ref ) = @_;
	my ($out_line, $total_nbps, $total_ngenes, $total_nbackbone) = ('', 0, 0, 0);
	my $curr_node_obj = $metacomponent_obj->node1;
	for(my $i=1; defined($metacomponent_obj->in_node($i)); $i++) {
		my $component_obj = $metacomponent_obj->in_node($i);
		my ($direction, $next_node_obj) = ($component_obj->node1->id eq $curr_node_obj->id)? (1, $component_obj->node2) : (-1, $component_obj->node1);
		my ($line, $nbps, $ngenes, $nbackbone) = $self->write_component($component_obj, $direction, ($i==1), $written_node_ids_ref);
		$out_line .= $line;
		$total_nbps += $nbps;
		$total_ngenes += $ngenes;
		$total_nbackbone += $nbackbone;
		$curr_node_obj = $next_node_obj;
	}
	return ($out_line, $total_nbps, $total_ngenes, $total_nbackbone);
}

#######################################################################################################################
sub write_component {
	my $self = shift;
	my ($component_obj, $direction, $include_first, $written_node_ids_ref ) = @_;

	my @backbone_node_objs = ();
	my %node_info = ();

	if($direction == 1) {
		push(@backbone_node_objs, $component_obj->node1);
		push(@backbone_node_objs, $component_obj->node2);
	}
	else {
		push(@backbone_node_objs, $component_obj->node2);
		push(@backbone_node_objs, $component_obj->node1);
	}
	my $last_node_obj = pop(@backbone_node_objs);
	defined($component_obj->node1) || die "\nError EExcelWriter0004: $component_obj->node1 is undefined\n\n";
	defined($component_obj->node2) || die "\nError EExcelWriter0005: $component_obj->node2 is undefined\n\n";
	$node_info{$component_obj->node1->id}{'obj'} = $component_obj->node1;
	$node_info{$component_obj->node2->id}{'obj'} = $component_obj->node2;

	for(my $i=1; defined($component_obj->in_node($i)); $i++) {
		my $node_obj = $component_obj->in_node($i);
		defined($node_obj) || die "\nError EExcelWriter0002: failed to fetch node $i for component " . $component_obj->id . " with " . $component_obj->size . " nodes\n\n";
		$node_info{$node_obj->id}{'obj'} = $node_obj;
		if(($component_obj->type eq 'linear') && ($node_obj->type eq 'in')) {
			if($direction == 1) { 
				push(@backbone_node_objs, $node_obj);
			}
			else {
				splice(@backbone_node_objs, 1, 0, $node_obj);
			}
		}
	}
	push(@backbone_node_objs, $last_node_obj);

	# Compute average length
	foreach my $node_id (keys %node_info) {
		my $node_obj = $node_info{$node_id}{'obj'};
		$node_info{$node_id}{'ngenes'} = $node_obj->size;
		my $avg_length = 0;
		foreach my $gene_obj (@{$node_obj->genes}) {
			$avg_length += $gene_obj->length;
		}
		$node_info{$node_id}{'avg_length'} = $avg_length/$node_obj->size; 
	}

	# Compute distance between backbone genes
	foreach my $i (1 .. $#backbone_node_objs) {
		my ($this_node_obj, $prev_node_obj) = ($backbone_node_objs[$i], $backbone_node_objs[$i-1]);
		my %seq2gene = ();
		foreach my $gene_obj (@{$this_node_obj->genes}) {
			$seq2gene{$gene_obj->seq->id} = $gene_obj;
		}
		my ($npairs, $intergenebp, $intergenegenes) = (0, 0, 0);
		foreach my $gene_obj (@{$prev_node_obj->genes}) {
			if(exists($seq2gene{$gene_obj->seq->id})) {
				$npairs++;
				my ($ngenes, $nbps) = $gene_obj->seq->gene_distance($gene_obj, $seq2gene{$gene_obj->seq->id});
				$intergenebp += $nbps;
				$intergenegenes += $ngenes;
			}
		}
		$node_info{$this_node_obj->id}{'npairs'} = ($npairs != 0)? $npairs : 0;
		$node_info{$this_node_obj->id}{'nbps'} = ($npairs != 0)? $intergenebp/$npairs : 0;
		$node_info{$this_node_obj->id}{'ngenes'} = ($npairs != 0)? $intergenegenes/$npairs : 0;
	}


	my $first_node_obj = shift(@backbone_node_objs);
	my ($out_line, $total_nbps, $total_ngenes, $total_nbackbone) = ('', 0, 0, scalar(@backbone_node_objs));

	if($include_first) {
		$out_line .= $first_node_obj->id . "\tBACKBONE\t" . $first_node_obj->strand . "\t" . int(10*$node_info{$first_node_obj->id}{'avg_length'})/10 . "\t" . $first_node_obj->size . "\n";
		$written_node_ids_ref->{$first_node_obj->id} = 1;
		$total_nbps += $node_info{$first_node_obj->id}{'avg_length'};
		$total_ngenes += 1;
		$total_nbackbone++;
	}
	delete($node_info{$first_node_obj->id});
	my $last_node_obj = pop(@backbone_node_objs);
	
	if($component_obj->type eq 'linear') {
		while(@backbone_node_objs > 0) {
			my $node_obj = shift(@backbone_node_objs);
			$out_line .= $node_obj->id . "\tBACKBONE\t" . $node_obj->strand . "\t" . int(10*$node_info{$node_obj->id}{'avg_length'})/10 . "\t" . $node_obj->size . "\t" .
					$node_info{$node_obj->id}{'npairs'} . "\t" . int(10*$node_info{$node_obj->id}{'ngenes'})/10 . "\t" . 
					int(10*$node_info{$node_obj->id}{'nbps'})/10 . "\n";
			$written_node_ids_ref->{$node_obj->id} = 1;
			$total_nbps += ($node_info{$node_obj->id}{'avg_length'}+$node_info{$node_obj->id}{'nbps'});
			$total_ngenes += $node_info{$node_obj->id}{'ngenes'};
			delete($node_info{$node_obj->id});
		}
	}
	foreach my $node_id (keys %node_info) {
		next if($node_id eq $last_node_obj->id);
		exists($node_info{$node_id}) || die "\nError EExcelWriter0003: could not find $node_id\n\n";
		exists($node_info{$node_id}{'obj'}) || die "\nError EExcelWriter0004: could not find object $node_id\n\n";
		my $node_obj = $node_info{$node_id}{'obj'};
		my $node_type = ($node_obj->type eq 'in')? "BUBBLE" : "TIP";
		$out_line .= $node_obj->id . "\t$node_type\t" . $node_obj->strand . "\t" . int(10*$node_info{$node_obj->id}{'avg_length'})/10 . "\t" . $node_obj->size . "\n"; 
		$written_node_ids_ref->{$node_obj->id} = 1;
	}
	$out_line .= $last_node_obj->id . "\tBACKBONE\t" . $last_node_obj->strand . "\t" . int(10*$node_info{$last_node_obj->id}{'avg_length'})/10 . "\t" . $last_node_obj->size . "\t" .
			$node_info{$last_node_obj->id}{'npairs'} . "\t" . int(10*$node_info{$last_node_obj->id}{'ngenes'})/10 . "\t" . 
			int(10*$node_info{$last_node_obj->id}{'nbps'})/10 . "\n";
	$total_nbps += $node_info{$last_node_obj->id}{'avg_length'} + $node_info{$last_node_obj->id}{'nbps'};
	$total_ngenes += $node_info{$last_node_obj->id}{'ngenes'};
	$written_node_ids_ref->{$last_node_obj->id} = 1;

	return ($out_line, $total_nbps, $total_ngenes, $total_nbackbone);
}

1;
