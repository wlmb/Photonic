package Photonic::Geometry::Geometry;
$Photonic::Geometry::Geometry::VERSION = '0.010';
use namespace::autoclean;
use Moose;

with 'Photonic::Roles::Geometry';

sub _forme {
    my $class=shift;
    my $args=shift;
    return defined $args->{B};
}

__PACKAGE__->meta->make_immutable; #faster execution

1;


=head1 NAME

Photonic::Geometry

=head1 VERSION

version 0.010

=head1 SYNOPSIS

     use Photonic::Geometry;
     $g=Photonic::Geometry->new(B=>$pdl);
     $B=$g->B;
     $G=$g->G;

=head1 DESCRIPTION

Create a geometry object to be used in a Photonic
calculation. 

You might want to use the related package
Photonic::Geometry::FromImage2D

=head1 METHODS

=over 4

=item * new(B=>$pdl[, L=>$L][, units=>$units])

Creates a new Ph::G object

$pdl is a boolean array with 1's and 0's representing the characteriztic
function within the unit cell. Its dimensions must be odd. Its number
of dimensions is the dimension of space

$L is the size of the unit cell along the cartesian axes. By
default, it is the number of pixels.

$units is an arrayref of pdl's, one for each basis vector of the
lattice. It defaults to cartesian unit vectors. B<Note:> Non-default
$units are not yet implemented.

=item * Vec2LC_G($v_G)

Returns the longitudinal component of a 'complex' vector field $v_G in
reciprocal space 

=item * LC2Vec_G($s_G)

longitudinal vector field from its longitudinal components in
reciprocal space. Scalar field to vector field.

=back

=head1 ACCESORS (read only)

=over 4

=item * B 

The characteristic function as PDL

=item * L 

Unit cell sizes as a B->ndims pdl.

=item * units

Basis C<e>_i of basis vectors for the lattice

=item * dims 

The dimensions [$X, $Y...] of the PDL B

=item * ndims 

The number of dimensions of the PDL B, i.e., the dimensionality of
space. 

=item * npoints 

Number of points within unit cell.

=item * scale

The distance between neighbor pixels along the axes.

=item * r

The position coordinates. In 2D, a 2,$X,$Y pdl. In 3D, a 3,$X,$Y,$Z pdl.

=item * G 

The reciprocal lattice. In 2D, a 2, $X, $Y pdl. B<G>.B<R>=multiple of 2\pi.

=item * GNorm

The reciprocal lattice, normalized to 1. In 2D, a 2, $X, $Y pdl. 

=item * f

Filling fraction of B region

=item * unitPairs

Normalized sum of pairs of unit vectors B<u>_{(ij)}=normalized
B<e>_ i+B<e>_j.

=item * unitDyads

Matrix of dyads of unit vector pairs
B<d>^{ij}_{kl}=B<u>^{i}_{kl}B<u>^{j}_{kl} as 2d matrix, adjusted for symmetry

=item * unitDyadsLU

LU decomposition of unitDyads. Used to obtain cartesian components of
dielectric tensor from longitudinal dielectric functions along the
directions given by unitPairs

=back

=head1 ACCESORS (read write)

=over 4

=item Direction0

Direction of the zero length wavevector

=back

=head1 PREDICATES

=over 4

=item has_Direction0

Test if Direction0 has been set.

=back

=for Pod::Coverage

=head2    PI


=cut
