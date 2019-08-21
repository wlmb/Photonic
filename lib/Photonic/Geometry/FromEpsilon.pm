package Photonic::Geometry::FromEpsilon;
$Photonic::Geometry::FromEpsilon::VERSION = '0.011';
use namespace::autoclean;
use Carp;
use Moose;
use MooseX::StrictConstructor;

has 'epsilon'=>(is=>'ro', isa=>'PDL::Complex', required=>1,
		documentation=>'Dielectric function as function of position');

has 'B' =>(is=>'ro', isa=>'PDL', init_arg=>undef, builder=>'_B', lazy=>1,
	   documentation=>'Charateristic function');

with 'Photonic::Roles::Geometry';

#Filling fraction is meaningless in this case.
before 'f' => sub {
    croak "Filling fraction is meaningless for "
	. "Photonic::Geometry::FromEpsilon"; 
};

sub _Bdims {
    my $self=shift;
    my $Bdims=return $self->epsilon->re; #value is irrelevant. Only shape counts.
}

__PACKAGE__->meta->make_immutable; #faster execution

1;


=head1 NAME

Photonic::Geometry::FromEpsilon

=head1 VERSION

version 0.011

=head1 SYNOPSIS

     use Photonic::Geometry::FromEpsilon;
     $g=Photonic::Geometry::FromEpsilon->new(epsilon=>$pdl);
     $epsilon=$g->epsilon;
     

=head1 DESCRIPTION

Creates a geometry object to be used in a Photonic
calculation from a dielectric function which specifies
the dielectric function for each media in the unit cell. 

=head1 METHODS

=over 4

=item * new(epsilon=>$pdl[, L=>$L])

Creates a new Ph::G::FromEpsilon object

$pdl is a boolean array with epsilon's representing the dielectric
function within the unit cell. Its dimensions must be odd. Its number
of dimensions of the real part is the dimension of space

$L is the size of the unit cell along the cartesian axes. By
default, it is the number of pixels.

=back

=head1 ACCESORS (read only)

=over 4

=item * epsilon

The dielectric function as PDL

=item * L 

Unit cell sizes as a Bdims->ndims pdl.

=back

=head1 SEE ALSO

L<Photonic::Roles::Geometry>

=begin Pod::Coverage

=head2    PI

=end Pod::Coverage

=cut
