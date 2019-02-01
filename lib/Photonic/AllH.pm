package Photonic::AllH;
$Photonic::AllH::VERSION='0.010';
use warnings;
use strict;
use Carp;

my @implementations=qw(LE::NR2::AllH WE::R2::AllH);
croak "Dont use Photonic::AllH. Use a specific implementation. " 
    . "Choose from: "
    . join ", ", map {"Photonic::OneH::" . $_} @implementations;
    
0;

=head1 NAME

Photonic::AllH

=head1 VERSION

version 0.010

=head1 SYNOPSIS

     use Photonic::LE::NR2::AllH;
     my $g1=Photonic::LE::NR2::AllH->new(geometry=>$geometry, nh=>$nh);

=head1 DESCRIPTION

Iterates the calculation of Haydock coefficients and states. 

=over 4

=item L<Photonic::LE::NR2::AllH>

Implementation for the longitudinal epsilon in a binary media in the
non retarded approximation.  

=item L<Photonic::WE::R2::AllH>

Implementation for the wave equation in a binary media, particles
within a non-dissipative host with retardation, using a metric.

=back

Consult the documentation of each implementation.

=cut
