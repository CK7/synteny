#!/usr/bin/perl

package Metacomponent;

use Component;
use strict;

our @ISA = qw(Component);

#######################################################################################################################
sub new {
	my $class = shift;

	my $self = $class->SUPER::new(@_);

	bless $self, $class;
	return $self;
}

1;
