package Photonic::WEM::ST::Metric;
$Photonic::WEM::S::Metric::VERSION = '0.024';

=encoding UTF-8

=head1 NAME

Photonic::WEM::S::Metric

=head1 VERSION

version 0.024

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

    use Photonic::WEM::S::Metric;
    my $gGG=Photonic::WEM::S::Metric->new(
            geometry=>$geometry, epsilon=>$eps,
            wavenumber => $q, $wavevector=>k);
    f($gGG->value);

=head1 DESCRIPTION

Calculates the retarded metric tensor g_{GG'}^{ij} for use in the
calculation of the retarded Haydock coefficients for the wave equation
in a binary medium where the host has no dissipation.

=head1 METHODS

=over 4

=item * new(geometry=>$g, epsilon=>$e, $wavenumber=>$q, $wavevector=>$k);

Create a new Ph::WEM::S::Metric object with Geometry $g, dielectric
function of the host $e, vacuum wavenumber $q=omega/c  and wavevector
$k. $q and $k are real.

= item * apply(psi)

Create and apply the magnetic metric to teh state psi provided in
Haydock when it calls this method.

=back


=cut

