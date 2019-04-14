package Photonic::WE;
$Photonic::WE::VERSION='0.011';

=head1 COPYRIGHT NOTICE 

Photonic - A perl package for calculations on photonics and
metamaterials. 

Copyright (C) 1916 by W. Luis Mochán

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 1, or (at your option)
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA  02110-1301 USA

    mochan@fis.unam.mx

    Instituto de Ciencias Físicas, UNAM
    Apartado Postal 48-3
    62251 Cuernavaca, Morelos
    México

=cut


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

version 0.011

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
