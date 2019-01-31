package Photonic::WE;
$Photonic::WE::VERSION='0.010';
use warnings;
use strict;
use Carp;

my @implementations=qw(R2::OneH);
croak "Dont use Photonic::OneH. Use a specific implementation. " 
    . "Choose from: "
    . join ", ", map {"Photonic::WE::" . $_} @implementations;
    
0;

=head1 NAME

Photonic::WE

=head1 VERSION

version 0.010

=head1 SYNOPSIS

     use Photonic::WE::R2::OneH;
     my $g1=Photonic::WE::R2::OneH->new(metric=>$m, polarization=>$p);

=head1 DESCRIPTION

Implements calculation of a Haydock coefficients and Haydock states for
the wave operator in an arbitrary number of dimensions. One Haydock
coefficient at a time.  

=over 4

=item L<Photonic::WE::R2::OneH>

Implementation for two media, one arbitrary medium within a
non-dissipative host with retardation, using a metric.

=back

Consult the documentation of each implementation.

=cut
