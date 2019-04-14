package Photonic::EL::S;
$Photonic::EL::S::VERSION='0.011';

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
