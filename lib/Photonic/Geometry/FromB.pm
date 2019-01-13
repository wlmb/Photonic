package Photonic::Geometry::FromB;
$Photonic::Geometry::FromB::VERSION = '0.010';
use namespace::autoclean;
use Moose;

has 'B' =>(is=>'ro', isa=>'PDL', required=>1,
	   documentation=>'charateristic function');

with 'Photonic::Roles::Geometry';


__PACKAGE__->meta->make_immutable; #faster execution

1;


=head1 NAME

Photonic::Geometry::FromB

=head1 VERSION

version 0.010

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
