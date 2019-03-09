package Photonic::EL::S;
$Photonic::EL::S::VERSION='0.011';
use warnings;
use strict;
use Carp;

my @implementations=qw(OneH);
croak "Dont use Photonic::EL::S. Use a specific implementation. " 
    . "Choose from: "
    . join ", ", map {"Photonic::EL::S::" . $_} @implementations;
    
0;
=head1 NAME

Photonic::EL::S

=head1 VERSION

version 0.011

=head1 SYNOPSIS

     use Photonic::EL::S::OneH;
     my $g1=Photonic::S::OneH->new(epsilon=>$e, geometry=>$g);

=head1 DESCRIPTION

Implements calculation of a Haydock coefficients and Haydock states for
the longitudinal dielectric function.

=over 4

=item L<Photonic::EL::S::OneH>

Implementation for arbitrary numbers of arbitrary media in the non
retarded approximation.  

=back

Consult the documentation of each implementation.

=cut
