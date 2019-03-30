package Photonic::EL::NP;
$Photonic::EL::NP::VERSION='0.011';
use warnings;
use strict;
use Carp;

my @implementations=qw(OneH);
croak "Dont use Photonic::EL::NP. Use a specific implementation. "
    . "Choose from: "
    . join ", ", map {"Photonic::EL::NP::" . $_} @implementations;

0;
=head1 NAME

Photonic::EL::NP

=head1 VERSION

version 0.011

=head1 SYNOPSIS

     use Photonic::EL::NP::OneH;
     my $g1=Photonic::EL::NP::OneH->new(epsilon=>$e, geometry=>$g);

=head1 DESCRIPTION

Implements calculation of a Haydock coefficients and Haydock states for
the longitudinal dielectric function for N media and a macroscopic
initial state.

=over 4

=item L<Photonic::EL::NP::OneH>

Implementation for arbitrary numbers of arbitrary media in the non
retarded approximation assuming the initial state is the macroscopic
=state.

=back

Consult the documentation of each implementation.

=cut
