package Photonic::Geometry;
$Photonic::Geometry::VERSION='0.011';

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

my @implementations=qw(FromB FromImage2D);
croak "Dont use Photonic::Geometry. Use a specific implementation. "
    . "Choose from: "
    . join ", ", map {"Photonic::Geometry::" . $_} @implementations;

0;

=head1 NAME

Photonic::Geometry

=head1 VERSION

version 0.011

=head1 SYNOPSIS

use Photonic::Geometry::FromB;
my $g1=Photonic::Geometry::FromB->new(B=>$pdl);
use Photonic::Geometry::FromImage2D;
my $g2=Photonic::Geometry::FromImage2D->new(path=>$filename);

=head1 DESCRIPTION

Create a geometry object used in several modules of Photonic.
This is a stub. You should choose a specific implementation. Currently
you could use

=over 4

=item L<Photonic::Geometry::FromB>

Initialize the geometry from a characteriztic function

=item L<Photonic::Geometry::FromImage2D>

Initialize the geometry from a 2D black and white image.

=back

Consult the documentation of each implementation.

=cut
