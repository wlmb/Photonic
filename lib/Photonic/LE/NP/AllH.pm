package Photonic::LE::NP::AllH;
$Photonic::LE::NP::AllH::VERSION = '0.017';

=encoding UTF-8

=head1 NAME

Photonic::LE::NP::AllH

=head1 VERSION

version 0.017

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

   use Photonic::LE::NP::AllH;
   my $iter=Photonic::LE::NP::AllH->new(
       epsilon=>$e, geometry=>$geometry, nh=>$Nh, keepStates=>$save);
   $iter->run;
   my $haydock_as=$iter->as;
   my $haydock_bs=$iter->bs;
   my $haydock_b2s=$iter->b2s;

=head1 DESCRIPTION

Iterates the calculation of Haydock coefficients and (optionally) states and saves
them for later retrieval. This class corresponds to the nonretarded longitudinal
dielectric function calculation for possibly more than two phases.

=head2 WARNING

This class only works for the default macroscopic initial
state. firstState should not be modified.

=head1 METHODS

=over 4

=item * new(epsilon=>$e, geometry=>$g, nh=>$nh, keepStates=>$k)

Initializes an Ph::NR::NP::AllH object. $nh is the maximum number of
desired coefficients, $k is a flag, non zero to save the Haydock
states. All other arguments are as in Photonic::LE::NP::OneH.

=item * run

Runs the iteration to completion.

=item * All the Photonic::LE::NP::OneH methods

Implements calculation of Haydock coefficients and states. See implementation documentation.

=item * All the Photonic::Roles::AllH methods

Iterates the calculation of one Haydock coefficient and state for
non-retarded system at a time and saves them to later retrival.

=back

=head1 ACCESSORS (read only)

=over 4

=item * nh

Maximum number of desired Haydock 'a' coefficients and states. The
number of b coefficients is one less.

=item * keepStates

Flag to keep (1) or discard (0) Haydock states

=item * states

Array of Haydock states

=item * as

Array of Haydock a coefficients

=item * bs

Array of Haydock b coefficients

=item * b2s

Array of Haydock b coefficients squared

=back

=begin Pod::Coverage

=head2 BUILD

=end Pod::Coverage

=cut

use namespace::autoclean;
use Moose;
use MooseX::StrictConstructor;

extends 'Photonic::LE::NP::OneH';
with 'Photonic::Roles::AllH', 'Photonic::Roles::ReorthogonalizeC';

__PACKAGE__->meta->make_immutable;

1;
