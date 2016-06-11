package Photonic;

use 5.006;
use strict;
use warnings;

=head1 NAME

Photonic - A perl package for calculations on photonics and metamaterials.

=head1 VERSION

Version 0.005

=cut
    
$Photonic::VERSION = '0.005';

=head1 SYNOPSIS

  use Photonic::Geometry;
  use Photonic::NonRetarded::EpsTensor;

  my $g=Photonic::Geometry->new(B=>$b);
  my $eps=Photonic::Nonretarded::EpsTensor->new(geometry=>$g);
  my $epsValue=$eps->evaluate($epsA, $epsB);

Calculates the dielectric tensor of a metamaterial made up of two
materials with dielectric functions $epsA and $epsB with a geometry $g
corresponding to a characteristic funcion $b.

=head1 DESCRIPTION

Set of packages for the calculation of optical properties of
metamaterials. The included packages are:

=over 4

=item L<Photonic::CharacteristicFunctions>

Some examples of characteriztic functions.

=item L<Photonic::ExtraUtils>

Some useful fortran routines from linpack.

=item L<Photonic::Geometry>

Object to hold the geometry of the metamaterial.

=item L<Photonic::Geometry::FromImage2D>

Obtain the geometry from a black and white 2D image.

=item L<Photonic::Types>

Definition of some types.

=item L<Photonic::Utils>

Some useful functions.

=item L<Photonic::NonRetarded::One>

Object to calculate one Haydock coefficient.

=item L<Photonic::Nonretarded::AllH>

Object to calculate all Haydock coeffficients.

=item L<Photonic::NonRetarded::EpsL>

Object to calculate the longitudinal dielectric response in the non
retarded regime.

=item L<Photonic::NonRetarded::EpsTensor>

Object to calculate all components of the dielectric tensor in the non
retarded regime.

=item L<Photonic::NonRetarded::FieldH>

Object to calculate the microscopic electromagnetic fields in the non
retarded regime.

=item L<Photonic::NonRetarded::SHP>

Object to prepare the data for the calculation of the SH polarization
in the non retarded regime using the dipolium model.

=item L<Photonic::NonRetarded::SH>

Object to calculate the second harmonic polarization in the
nonretarded regime.

=item L<Photonic::NonRetarded::SHChiTensor>

Object to calculate the second harmonic susceptibility tensor in the
non retarded regime using the dipolium model.


=item L<Photonic::Roles::EpsParams>

Fields that are factored as they are common in different packages.

=item L<Photonic::Roles::KeepStates>

Flag that is factored as it is common in different packages.

=back

=cut
    
=head1 ACKNOWLEDGMENTS

This work was partially supported by DGAPA-UNAM under grants IN108413
and IN113016.   

=cut
    
1;
