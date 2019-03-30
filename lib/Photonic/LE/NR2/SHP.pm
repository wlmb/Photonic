=head1 NAME

Photonic::LE::NR2::SHP

=head1 VERSION

version 0.011

=head1 SYNOPSIS

   use Photonic::LE::NR2::SHP;
   my $nrshp=Photonic::LE::NR2::SHP->
             new(nrf=>$nrf, densityA=>$dA, densityB=>$dB));



=head1 DESCRIPTION

Prepares the data for the calculation of the non retarded SH
polarization of an arbitrary periodic composite made up of
centrosymmetric isotropic component materials, using the continuous
dipolium model.

=head1 METHODS

=over 4

=item * new(nrf=>$nrf, densityA=>$dA, densityB=>$dB)

Initializes the structure

$nrf Photonic::LE::NR2::Field is a Haydock field calculator for the
structure.

$dA is the density of polarizable entities in medium A

$dB is the density of polarizable entities in medium B

=back

=head1 ACCESORS (read only)

=over 4

=item * nrf

Photonic::LE::NR2::Field Haydock field calculator

=item * densityA, densityB

Normalized (to what?) dipole entities density in media A and B

=item * density

Density field over unit cell

=item * ndims

Number of dimensions of the system

=back

=begin Pod::Coverage

=head2 BUILD

=end Pod::Coverage

=cut

package Photonic::LE::NR2::SHP;
$Photonic::LE::NR2::SHP::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use PDL::FFTW3;
use PDL::Constants qw(PI);
use Photonic::Utils qw(RtoG HProd linearCombine);
use Photonic::ExtraUtils qw(cgtsl);
use Moose;
use MooseX::StrictConstructor;


has 'nrf'=>(is=>'ro', isa=>'Photonic::LE::NR2::Field', required=>1,
         documentation=>'Haydock field calculator');
has 'densityA'=>(is=>'ro', isa=>'Num', required=>1,
         documentation=>'Normalized dipole entities density in medium A');
has 'densityB'=>(is=>'ro', isa=>'Num', required=>1,
         documentation=>'Normalized dipole entities density in medium B');
has 'density'=>(is=>'ro', isa=>'PDL', writer=>'_density', init_arg=>undef,
         documentation=>'Normalized dipole entities density over unit cell');
has 'ndims' =>(is=>'ro', isa=>'Int', init_arg=>undef, lazy=>1,
         builder=>'_ndims',
         documentation=>'Number of dimensions of system');

sub BUILD {
    my $self=shift;
    my $B=$self->nrf->nr->B;
    $self->_density($self->densityA*(1-$B)+$self->densityB*$B);
}

sub _ndims {
    my $self=shift;
    return $self->nrf->nr->B->ndims;
}

__PACKAGE__->meta->make_immutable;

1;
