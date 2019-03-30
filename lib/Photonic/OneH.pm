package Photonic::OneH;
$Photonic::OneH::VERSION='0.011';
use warnings;
use strict;
use Carp;

my @implementations=qw(LE::NR2::OneH WE::R2::OneH);
croak "Dont use Photonic::OneH. Use a specific implementation. "
    . "Choose from: "
    . join ", ", map {"Photonic::OneH::" . $_} @implementations;

0;

=head1 NAME

Photonic::OneH

=head1 VERSION

version 0.011

=head1 SYNOPSIS

     use Photonic::LE::NR2::OneH;
     my $g1=Photonic::LE::NR2::OneH->new(geometry=>$geometry);
     use Photonic::WE::R2::OneH;
     my $g1=Photonic::WE::R2::OneH->new(metric=>$m, polarization=>$p);

=head1 DESCRIPTION

Implements calculation of a Haydock coefficients and Haydock states for
a non retarded Photonic calculation corresponding to a binary medium
in an arbitrary number of dimensions of the non retarded dielectric
function of arbitrary  periodic two component systems in arbitrary
=number of dimentions. One Haydock coefficient at a time.

=over 4

=item L<Photonic::LE::NR2::OneH>

Implementation for the longitudinal epsilon in a binary media in the
non retarded approximation.

=item L<Photonic::WE::R2::OneH>

Implementation for the wave equation in a binary media, particles
=within a non-dissipative host with retardation, using a metric.

=back

Consult the documentation of each implementation.

=cut
