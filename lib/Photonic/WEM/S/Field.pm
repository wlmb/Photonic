package Photonic::WEM::S::Field;
$Photonic::WEM::S::Field::VERSION = '0.024_01';

=encoding UTF-8

=head1 NAME

Photonic::WEM::S::Field

=head1 VERSION

version 0.024_01

=head1 COPYRIGHT NOTICE

Photonic - A perl package for calculations on photonics and
metamaterials.

Copyright (C) 2016 by W. Luis Mochán

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

=head1 SYNOPSIS

   use Photonic::WEM::S::Field;
   my $nrf=Photonic::WE::S::Field->new(...);
   my $field=$nrf->field;
   my $rawfield=$nrf->rawfield;

=head1 DESCRIPTION

Calculates the retarded electric field for a given fixed
Photonic::Geometry structure and given dielectric functions of
the components.

Consumes L<Photonic::Roles::Field>
- please see for attributes.

=cut

use namespace::autoclean;
use Moo;
#use MooX::StrictConstructor;
extends "Photonic::WE::S::Field";

__PACKAGE__->meta->make_immutable;

1;
