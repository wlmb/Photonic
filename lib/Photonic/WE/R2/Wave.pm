package Photonic::WE::R2::Wave;
$Photonic::WE::R2::Wave::VERSION = '0.021';

=encoding UTF-8

=head1 NAME

Photonic::WE::R2::Wave

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

   use Photonic::WE::R2::Wave;
   my $W=Photonic::WE::R2::Wave->new(metric=>$m, nh=>$nh, epsB=>$epsB);
   my $WaveTensor=$W->waveOperator;

=head1 DESCRIPTION

Calculates the macroscopic wave operator for a given fixed
Photonic::WE::R2::Metric structure as a function of the dielectric
functions of the components.

Extends L<Photonic::WE::R2::Green>, please refer.

=head1 ATTRIBUTES

=over 4

=item * waveOperator

Returns the macroscopic wave operator for the dielectric functions of the
particle C<epsB>. The host's response C<epsA> is taken from the metric.

=back

=cut

use namespace::autoclean;
use Photonic::Utils qw(wave_operator);
use Photonic::Types;
use Moose;
use MooseX::StrictConstructor;

extends 'Photonic::WE::R2::Green';

has 'waveOperator' =>  (is=>'ro', isa=>'Photonic::Types::PDLComplex', init_arg=>undef,
             lazy=>1, builder=>'_build_waveOperator',
             documentation=>'Wave operator from last evaluation');

sub _build_waveOperator {
    my $self=shift;
    wave_operator($self->greenTensor, $self->geometry->ndims);
}

__PACKAGE__->meta->make_immutable;

1;
