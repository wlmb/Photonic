=head1 NAME

Photonic::WE::S::Metric

=head1 VERSION

version 0.010

=head1 SYNOPSIS

    use Photonic::WE::S::Metric;
    my $gGG=Photonic::WE::S::Metric->new(
            geometry=>$geometry, epsilon=>$eps, 
            wavenumber => $q, $wavevector=>k);
    f($gGG->value);

=head1 DESCRIPTION

Calculates the retarded metric tensor g_{GG'}^{ij} for use in the
calculation of the retarded Haydock coefficients for the wave equation in a binary medium where the host has no dissipation.

=head1 METHODS

=over 4

=item * new(geometry=>$g, epsilon=>$e, $wavenumber=>$q, $wavevector=>$k);

Create a new Ph::WE::R2::Metric object with Geometry $g, dielectric
function of the host $e, vacuum wavenumber $q=omega/c  and wavevector
$k. $q and $k are real.

=back

=head1 ACCESORS (read only)

=over 4

=item * value 

The actual metric tensor as a complex PDL (d,d,n1,n2..nd)
the first and second indices over cartesian indices for 0 to d-1 in d
dimensions, the next d indices n1,n2...nd identify the wavevector G.

=back

=cut

package Photonic::WE::S::Metric;
$Photonic::WE::R2::Metric::VERSION = '0.010';
use namespace::autoclean;
use PDL::Lite;
use PDL::MatrixOps;
use PDL::NiceSlice;
use PDL::Complex;
use List::Util;
use Carp;
use Moose;
use Photonic::Types;

# Later make it complex
has 'value'     => (is=>'ro', isa=>'PDL', init_arg=>undef, lazy=>1,
                   builder=>'_value', 
                   documentation=>'Metric tensor');

with 'Photonic::Roles::Metric';

sub _value {
    # Evaluate the metric tensor. Eq. MQ83.
    my $self = shift;
    my $G=$self->G; #reciprocal lattice;
    # check epsilon, wavenumber, wavevector real
    my $q=$self->wavenumber;
    my $eps=$self->epsilon;
    my $k=$self->wavevector;
    croak "Wave vector must be ".$self->ndims."-dimensional vector" unless
	[$k->dims]->[0]==$self->ndims and $k->ndims==1;
    #might generalize to complex k
    my ($kPG, $kMG) = ($k+$G, $k-$G); #xy:nx:ny
    # (k+G)(k+G) diad
    my ($kPGkPG, $kMGkMG) = map {$_->outer($_)} ($kPG, $kMG); #xy:xy:nx:ny
    # interior product
    my ($kPG2,$kMG2) = map {$_->inner($_)}  ($kPG, $kMG); #nx:ny;
    my $id=identity($self->ndims);
    my $k02=$eps*$q*$q; # squared wavenumber in 'host'
    #xyz xyz nx ny nz
    #cartesian matrix for each wavevector.
    my $gPGG=($k02*$id-$kPGkPG)/(($k02-$kPG2)->(*1,*1)); #xy:xy:nx:ny
    my $gMGG=($k02*$id-$kMGkMG)/(($k02-$kMG2)->(*1,*1)); #xy:xy:nx:ny
    my $gGG=PDL->pdl($gPGG, $gMGG)->mv(-1,0); #pm:xy:xy:nx:ny
    return $gGG;
}
    

1;
