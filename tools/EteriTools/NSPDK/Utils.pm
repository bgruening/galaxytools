package Array::Utils;

=head1 NAME

Array::Utils - small utils for array manipulation

=head1 SYNOPSIS

	use Array::Utils qw(:all);
	
	my @a = qw( a b c d );
	my @b = qw( c d e f );

	# symmetric difference
	my @diff = array_diff(@a, @b);

	# intersection
	my @isect = intersect(@a, @b);
	
	# unique union
	my @unique = unique(@a, @b);
	
	# check if arrays contain same members
	
	if ( !array_diff(@a, @b) ) {
		# do something
	}
	
=head1 DESCRIPTION

A small pure-perl module containing list manipulation routines. The module
emerged because I was tired to include same utility routines in numerous projects.

=head1 FUNCTIONS

=over 4

=item C<unique>

Returns an array of unique items in the arguments list.

=item C<intersect>

Returns an intersection of two arrays passed as arguments.

=item C<array_diff>

Return symmetric difference of two arrays passed as arguments.

=back

=head1 BUGS

None known yet

=head1 AUTHOR

Sergei A. Fedorov <zmij@cpan.org>

I will be happy to have your feedback about the module.

=head1 COPYRIGHT

This module is Copyright (c) 2007 Sergei A. Fedorov.
All rights reserved.

You may distribute under the terms of either the GNU General Public
License or the Artistic License, as specified in the Perl README file.

=head1 WARRANTY

This is free software. IT COMES WITHOUT WARRANTY OF ANY KIND.

=cut

use strict;

require Exporter;
our @ISA = qw(Exporter);

use Data::Dumper;

our %EXPORT_TAGS = (
	all	=> [ qw(
		&unique
		&intersect
		&array_diff
	) ],
);
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.3';

sub unique(@) {
	return keys %{{map {$_=>1} @_}};
}

sub intersect(\@\@) {
	my %e = map {$_=>1} @{$_[0]};
	return grep { $e{$_} } @{$_[1]};
}

sub array_diff(\@\@) {
	my %e = map {$_=>1} @{$_[1]};
	return @{[ ( grep { !delete $e{$_} } @{ $_[0] } ), keys %e ] };
}

1;
