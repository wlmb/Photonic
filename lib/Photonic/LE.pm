#Longitudinal Epsilon
package Photonic::LE;
$Photonic::LE::VERSION='0.011';
use warnings;
use strict;
use Carp;

my @implementations=qw(NR2::OneH);
croak "Dont use Photonic::LE. Use a specific implementation. "
    . "Choose from: "
    . join ", ", map {"Photonic::LE::" . $_} @implementations;

0;
=head1 NAME

Photonic::LE

=head1 VERSION

version 0.011

=head1 SYNOPSIS

     use Photonic::LE::NR2::OneH;
     my $g1=Photonic::NR2::OneH->new(geometry=>$geometry);

=head1 DESCRIPTION

Implements calculation of a Haydock coefficients and Haydock states for
the longitunidal dielectric function.

=over 4

=item L<Photonic::LE::NR2::OneH>

Implementation for two arbitrary media in the non retarded approximation.

=back

Consult the documentation of each implementation.

=cut
