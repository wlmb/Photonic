package Photonic::OneH;
$Photonic::OneH::VERSION='0.010';
use warnings;
use strict;
use Carp;

my @implementations=qw(NR2 R2);
croak "Dont use Photonic::OneH. Use a specific implementation. " 
    . "Choose from: "
    . join ", ", map {"Photonic::OneH::" . $_} @implementations;
    
0;

=head1 NAME

Photonic::OneH

=head1 VERSION

version 0.010

=head1 SYNOPSIS

     use Photonic::OneH::NR2;
     my $g1=Photonic::OneH::NR2->new(geometry=>$geometry);
     use Photonic::OneH::R2;
     my $g1=Photonic::OneH::R2->new(metric=>$m, polarization=>$p);

=head1 DESCRIPTION

Implements calculation of a Haydock coefficients and Haydock states for
a non retarded Photonic calculation corresponding to a binary medium
in an arbitrary number of dimensions of the non retarded dielectric
function of arbitrary  periodic two component systems in arbitrary
=number of dimentions. One Haydock coefficient at a time.

=over 4

=item L<Photonic::OneH::NR2>

Implementation for two arbitrary media in the non retarded approximation.

=item L<Photonic::OneH::R2>

Implementation for one arbitrary medium within a non-dissipative host
with retardation, using a metric.

=back

Consult the documentation of each implementation.

=cut
