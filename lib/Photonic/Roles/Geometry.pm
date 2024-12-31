package Photonic::Roles::Geometry;
$Photonic::Roles::Geometry::VERSION = '0.023';

=encoding UTF-8

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


use Moo::Role;

use PDL::Lite;
use PDL::NiceSlice;
use Photonic::Types -all;
use Photonic::Utils qw(any_complex lu_decomp make_dyads triangle_coords);
use PDL::MatrixOps qw();
use Carp;
use constant PI=>4*atan2(1,1);

requires 'B'; #characteristic function

has 'L' =>(is=>'lazy', isa => PDLObj,
	   documentation=>'array of unit cell size');
has 'units'=>(is=>'lazy', isa=>PDLObj,
     documentation=>'Basis of unit vectors');
has 'primitive'=>(is=>'lazy', isa=>PDLObj, required=>1, # required?
		  documentation=>'Primitive directions');
has 'primitiveNorm'=>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
		  documentation=>'Normalized primitive vectors');
has 'dualNorm'=>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
		 documentation=>'Normalized reciprocal primitive vectors');
has 'dims' =>(is=>'lazy', isa=>ArrayRef[Int],
           init_arg=>undef,
           documentation=>'list of dimensions of B');
has 'ndims' =>(is=>'lazy', isa=>Int,
           init_arg=>undef,
           documentation=>'number of dimensions of B');
has 'npoints' =>(is=>'lazy', isa=>Int,
           init_arg=>undef,
           documentation=>'number of points within B');
has 'scale'=>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
	      documentation=>'distances between pixels');
has 'r' =>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
	   documentation=>'array of positions x_or_y, nx, ny');
has 'G' =>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
	   documentation=>'array of reciprocal vectors xy,nx,ny...');
has 'GNorm' =>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
     documentation=>'array of unit norm reciprocal vectors xy,nx,ny...');
has 'mGNorm' =>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
     documentation=>
	  'array of negated unit norm reciprocal vectors xy,nx,ny...');
has 'pmGNorm' =>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
     documentation=>
       'array of spinors of +- unit norm reciprocal vectors xy,pm,nx,ny');
has 'f'=>(is=>'lazy', init_arg=>undef,
     documentation=>'filling fraction of B region');
has 'unitPairs'=>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
     documentation=>'Normalized sum of pairs of basis vectors');
has 'cUnitPairs'=>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
     documentation=>'Normalized complex sum of pairs of basis vectors');
has 'unitDyads'=>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
     documentation=>'Matrix of dyads of unit vector pairs');
has 'unitDyadsLU'=>(is=>'lazy', isa=>ArrayRef,
     documentation=>'LU decomposition of unitDyads');
has 'Direction0' =>(is => 'ro', isa => PDLObj,
     predicate=>'has_Direction0');

sub _build_L {
    my $self=shift;
    my $L=PDL->pdl($self->B->dims);
    return $L;
}

sub _build_units {
    my $self=shift;
    PDL::MatrixOps::identity($self->B->ndims); # unit vectors
}

sub _build_primitive {
    my $self=shift;
    $self->units->copy; #xy:n
}

sub _build_primitiveNorm {
    my $self=shift;
    return $self->primitive->norm;
}

sub _build_dualNorm {
    my $self=shift;
    my $pn=$self->primitiveNorm; #xy:n
    my $dual=$pn->inv->transpose->norm; #xy:n
}

sub _build_dims {
    my $self=shift;
    return [$self->B->dims];
}

sub _build_ndims {
    my $self=shift;
    return $self->B->ndims;
}

sub _build_npoints {
    my $self=shift;
    return $self->B->nelem;
}

sub _build_scale {
    my $self=shift;
    return $self->L/PDL->pdl($self->dims);
}

sub _build_r {
    my $self=shift;
    my $coords=$self->B->ndcoords; #i:nx:ny
    my $pr=$self->scale->(*1)*$self->primitiveNorm; #xy:i
    my $r=(
	$coords->(:,*1) #i:*:nx:ny
	*$pr->transpose #i:xy
	)->sumover; #xy:nx:ny
    return $r;
}

sub _build_G {
    my $self=shift;
    my @N=map {($_-($_%2))/2} @{$self->dims}; #N
    my $coords=($self->B->ndcoords)-(PDL->pdl(@N)); #vector from center.
    my $rec=(2*PI/$self->L)->(*1)*$self->dualNorm;
    my $G=(
	$coords->(:,*1) #i:*:nx:ny
	*$rec->transpose #i:xy
	)->sumover; #xy:nx:ny
    foreach(1..$G->ndims-1){
	$G=$G->mv($_,0)
	    ->rotate(-$N[$_-1])->mv(0,$_); #move origin to 0,0
    }
    return $G;
}


sub _build_GNorm { #origin set to zero here.
    my $self=shift;
    confess "Can't normalize reciprocal lattice unless Direction0 is set"
	unless $self->has_Direction0;
    return $self->G->norm;
}

sub _build_mGNorm { #normalized negated reciprocal lattice. Leave
    #direction 0 invariant.
    my $self=shift;
    return -$self->GNorm;
}
sub _build_pmGNorm { #normalized +- reciprocal lattice. Leave
    #direction 0 invariant.
    #xy,pm,nx,ny
    my $self=shift;
    return PDL->pdl($self->GNorm, $self->mGNorm)->mv(-1,1);
}

sub BUILD {
  my ($self) = @_;
  if ($self->has_Direction0) {
    my $value = $self->Direction0;
    confess "Direction0 must be of dimension ".$self->ndims." vector" unless
      $value->dim(0)==$self->ndims and $value->ndims==1;
    confess "Direction must be non-null" unless $value->inner($value)>0;
    my $arg=":". (",(0)" x $self->ndims); #:,(0),... dimension of space times
    $value=$value->norm; #normalize
    #Set element 0,0 for normalized arrays.
    $self->GNorm->slice($arg).=$value; #Normalized 0.
    $self->mGNorm->slice($arg).=$value; #Don't change sign for mGNorm!
    $self->pmGNorm->mv(1,-1)->slice($arg).=$value;
  }
  1; # cover -test causes a boolean test on the return value so make not pdl
}

sub _build_f { #calculate filling fraction
    my $self=shift;
    return $self->B->davg;
}

sub _build_unitPairs {
    my $self=shift;
    my $nd=$self->ndims;
    my $units=$self->units;
    my $pairs = PDL->zeroes($nd, $nd * ($nd + 1) / 2);
    my $pairs_slice = $pairs->mv(1, 0);
    my $indexes = triangle_coords($nd, 1);
    $pairs_slice .= (
      $units->indexND($indexes(0))+
      $units->indexND($indexes(1))
    )->transpose->norm->transpose;
    $pairs;
}

sub _build_cUnitPairs {
    my $self=shift;
    my $nd=$self->ndims;
    my $units=$self->units;
    my $cpairs = PDL->zeroes($nd, $nd * ($nd - 1) / 2)->r2C;
    my $cpairs_slice = $cpairs->mv(1, 0);
    my $indexes = triangle_coords($nd);
    $cpairs_slice .= PDL::czip(
      $units->indexND($indexes(0)),
      $units->indexND($indexes(1))
    );
    $cpairs_slice /= sqrt($cpairs_slice->abs2->transpose->sumover);
    $cpairs;
}

sub _build_unitDyads {
    my $self=shift;
    make_dyads($self->ndims, $self->unitPairs);
}

sub _build_unitDyadsLU {
    my $self=shift;
    [lu_decomp($self->unitDyads->r2C)];
}

sub Vec2LC_G { #longitudinal component of 'complex' vector field in
               #reciprocal space
    my $self=shift;
    my $field=shift; # vector field to project
    confess "Can't project unless Direction0 is set" unless
	$self->has_Direction0;
    my $gnorm=$self->GNorm;
    ($field * $gnorm->r2C)->sumover;
}

sub LC2Vec_G { #longitudinal vector field from its longitudinal
	       #components in reciprocal space
    my $self=shift;
    my $field=shift; # scalar field of longitudinal components
    confess "Can't project unless Direction0 is set" unless
	$self->has_Direction0;
    my $gnorm=$self->GNorm;
    #$gnorm is XorY nx, ny...
    #$field is nx ny nz
    #$result XoY nx ny nz
    $field->(*1)*$gnorm;
}

no Moo::Role;

1;


=head1 NAME

Photonic::Roles::Geometry

=head1 VERSION

version 0.023

=head1 SYNOPSIS

    use Photonic::Geometry::FromB;
    $g=Photonic::Geometry->new(B=>$pdl);
    $B=$g->B;
    $G=$g->G;

=over 4

=item (for developers)

    package Photonic::Geometry::FromB;
    $Photonic::Geometry::Geometry::VERSION = '0.023';
    use namespace::autoclean;
    use Moo;
    has 'B' =>(is=>'ro', isa=>PDLObj, required=>1,
	       documentation=>'characteristic function');
    with 'Photonic::Roles::Geometry';

=back

=head1 DESCRIPTION

Roles consumed by geometry objects to be used in a Photonic
calculation. See also the specific implementations under
L<Photonic::Geometry>.

=head1 ACCESSORS (defined by the implementation)

=over 4

=item * B

The characteristic function as PDL

=back

=head1 ACCESSORS (read write)

=over 4

=item * Direction0

Direction of the zero length wavevector

=back

=head1 ACCESSORS (read only)

=over 4

=item * L

Unit cell sizes as a ndims pdl.

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

=head1 METHODS

=over 4

=item * Vec2LC_G($v_G)

Returns the longitudinal component of a 'complex' vector field $v_G in
reciprocal space.

=item * LC2Vec_G($s_G)

Longitudinal vector field from its longitudinal components in
reciprocal space. Scalar field to vector field.

=back

=head1 PREDICATES

=over 4

=item has_Direction0

Test if Direction0 has been set.

=back

=begin Pod::Coverage

=head2 PI

=head2 BUILD

=end  Pod::Coverage

=cut
