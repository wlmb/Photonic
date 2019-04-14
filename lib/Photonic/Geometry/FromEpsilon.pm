package Photonic::Geometry::FromEpsilon;
$Photonic::Geometry::FromEpsilon::VERSION = '0.011';

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


use namespace::autoclean;
use Carp;
use Moose;
use MooseX::StrictConstructor;

has 'epsilon'=>(is=>'ro', isa=>'PDL::Complex', required=>1,
		documentation=>'Dielectric function as function of position');

has 'B' =>(is=>'ro', isa=>'PDL', init_arg=>undef, builder=>'_B', lazy=>1,
	   documentation=>'charateristic function');

with 'Photonic::Roles::Geometry';

#Filling fraction is meaningless in this case.
before 'f' => sub {
    croak "Filling fraction is meaningless for "
	. "Photonic::Geometry::FromEpsilon";
};

sub _B {
    my $self=shift;
    my $B=return $self->epsilon->re; #value is irrelevant. Only shape counts.
}

__PACKAGE__->meta->make_immutable; #faster execution

1;


=head1 NAME

Photonic::Geometry::FromB

=head1 VERSION

version 0.011

=head1 SYNOPSIS

     use Photonic::Geometry::FromB;
     $g=Photonic::Geometry::FromB->new(B=>$pdl);
     $B=$g->B;
     $G=$g->G;

=head1 DESCRIPTION

Create a geometry object to be used in a Photonic
calculation from a characteristic function.

=head1 METHODS

=over 4

=item * new(B=>$pdl[, L=>$L])

Creates a new Ph::G::FromB object

$pdl is a boolean array with 1's and 0's representing the characteriztic
function within the unit cell. Its dimensions must be odd. Its number
of dimensions is the dimension of space

$L is the size of the unit cell along the cartesian axes. By
default, it is the number of pixels.

=back

=head1 ACCESORS (read only)

=over 4

=item * B

The characteristic function as PDL

=item * L

Unit cell sizes as a B->ndims pdl.

=back

=head1 SEE ALSO

L<Photonic::Roles::Geometry>

=begin Pod::Coverage

=head2    PI

=end Pod::Coverage

=cut
