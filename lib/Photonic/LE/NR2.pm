package Photonic::EL::NR2;
$Photonic::EL::NR2::VERSION='0.010';
use warnings;
use strict;
use Carp;

my @implementations=qw(OneH);
croak "Dont use Photonic::EL::NR2. Use a specific implementation. " 
    . "Choose from: "
    . join ", ", map {"Photonic::EL::NR2::" . $_} @implementations;
    
0;
=head1 NAME

Photonic::EL::NR2

=head1 VERSION

version 0.010

=head1 SYNOPSIS

     use Photonic::EL::NR2::OneH;
     my $g1=Photonic::NR2::OneH->new(geometry=>$geometry);

=head1 DESCRIPTION

Implements calculation of a Haydock coefficients and Haydock states for
the longitunidal dielectric function.

=over 4

=item L<Photonic::EL::NR2::OneH>

Implementation for two arbitrary media in the non retarded approximation.

=back

Consult the documentation of each implementation.

=cut
