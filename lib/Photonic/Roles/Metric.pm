=head1 NAME

Photonic::Roles::Metric

=head1 VERSION

version 0.011

=head1 SYNOPSIS

    use Photonic::WE::R2::Metric;
    my $gGG=Photonic::WE::R2::Metric->new(
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
function of the reference $e, vacuum wavenumber $q=omega/c  and wavevector
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

package Photonic::Roles::Metric;
$Photonic::Roles::Metric::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use Moose::Role;
use Photonic::Types;

has 'geometry'  => (is=>'ro', isa=>'Photonic::Types::Geometry', required=>1,
                    handles=>[qw(B dims ndims r G GNorm L scale f)],
                    required=>1,
                    documentation=>'Geometry');
has 'epsilon'   => (is=>'ro', isa=>'PDL', required=>1,
		    default=>sub{PDL->pdl(1)}, 
                   documentation=>'Real reference dielectric function');
has 'wavenumber'=> (is=>'ro', isa=>'PDL', required=>1,
                   documentation=>'Vacuum wavenumber w/c');
has 'wavevector'=> (is=>'ro', isa=>'PDL', required=>1,
                   documentation=>'Wave vector');
requires qw(value); #provided by metric instances

1;
