package Photonic::WE::R2::WaveP;
$Photonic::WE::R2::WaveP::VERSION = '0.021';

=encoding UTF-8

=head1 NAME

Photonic::WE::R2::WaveP

=head1 VERSION

version 0.021

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

   use Photonic::WE::R2::WaveP;
   my $W=Photonic::WE::R2::WaveP->new(haydock=>$h, nh=>$nh, smallE=>$s, epsB=>$epsB);
   my $WaveProjection=$W->waveOperator;

=head1 DESCRIPTION

Calculates the macroscopic projected wave operator for a given fixed
L<Photonic::WE::R2::Haydock> structure as a function of the dielectric
functions of the components.

NOTE: Only works along principal directions, as it treats Green's
function as scalar.

Extends L<Photonic::WE::R2::GreenP>, please refer.

=head1 ATTRIBUTES

=over 4

=item * waveOperator

Returns the macroscopic wave operator for a given value of the
dielectric functions of the particle C<epsB>. The host's
response C<epsA> is taken from the metric.

=back

=cut

use namespace::autoclean;
use PDL::Lite;
use Photonic::Types;
use Moose;
use MooseX::StrictConstructor;

extends 'Photonic::WE::R2::GreenP';

has 'waveOperator' =>  (is=>'ro', isa=>'Photonic::Types::PDLComplex', init_arg=>undef,
             lazy=>1, builder=>'_build_waveOperator',
             documentation=>'Wave operator from last evaluation');

sub _build_waveOperator {
    1/shift->Gpp; #only works along principal directions!!
};

__PACKAGE__->meta->make_immutable;

1;
