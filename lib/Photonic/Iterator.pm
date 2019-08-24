package Photonic::Iterator;
$Photonic::Iterator::VERSION='0.011';
use base Exporter;
@EXPORT_OK=qw(iterator nextval);
%EXPORT_TAGS=(all=>\@EXPORT_OK);

sub nextval ($) { $_[0]->(); }
sub iterator(&) { return $_[0] }


=head1 NAME

Photonic::Iterator

=head1 VERSION

version 0.011

=head1 SYNOPSIS

    use Photonic::Iterator qw(iterator nextval);

=head1 DESCRIPTION

Iterator is an object interface to a list. Generates elements of a list as are solicited so that just a little part of the list will need to be in memory.
Uses nextval to replace the current element for the next one.

=back

=begin Pod::Coverage

=head2 BUILD
   
=end Pod::Coverage

=cut
