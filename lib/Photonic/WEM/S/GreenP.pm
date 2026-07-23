package Photonic::WEM::S::GreenP;
$Photonic::WEM::S::GreenP::VERSION = '0.02401';

=encoding UTF-8

=head1 NAME

Photonic::WE::S::GreenP

=head1 VERSION

version 0.02401

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

   use Photonic::WEM::S::GreenP;
   my $green=Photonic::WE::S::GreenP->new(haydock=>$h, nh=>$nh);
   my $greenProjection=$green->Gpp;
   my $WaveProjection=$green->waveOperator;
   my $EpsProjection=$green->epsilon;

=head1 DESCRIPTION

Calculates the dielectric function for a given fixed
L<Photonic::WE::ST::Haydock> structure as a function of the dielectric
functions of the components.

=head1 ATTRIBUTES

=over 4

=item * haydock

The L<Photonic::WE::ST::Haydock> structure (required).

=item * nh

The maximum number of Haydock coefficients to use.

=item * smallE

Criteria of convergence. 0 means don't check. (defaults to 1e-7)

=item * u

The spectral variable used in the calculation

=item * nhActual

The actual number of Haydock coefficients used in the calculation

=item * converged

Flags that the calculation converged before using up all coefficients

=item * waveOperator

The macroscopic wave operator calculated from the metric.

NOTE: Only works along principal directions, as it treats Green's
function as scalar.

=item * epsilon

The macroscopic dielectric projection

NOTE: Only works for polarizations along principal directions.

=back

=cut

use namespace::autoclean;
use Moo;
#use MooX::StrictConstructor;
extends "Photonic::WE::S::GreenP";

__PACKAGE__->meta->make_immutable;

1;
