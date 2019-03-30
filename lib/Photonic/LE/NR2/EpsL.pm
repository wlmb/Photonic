=head1 NAME

Photonic::LE::NR2::EpsL

=head1 VERSION

version 0.011

=head1 SYNOPSIS

   use Photonic::LE::NR2::EpsL;
   my $eps=Photonic::LE::NR2::EpsL->new(nr=>$nr, nh=>$nh);
   my $epsilonLongitudinal=$eps->evaluate($epsA, $epsB);

=head1 DESCRIPTION

Calculates the longitudinal dielectric function for a given fixed
Photonic::LE::NR2::AllH structure as a function of the dielectric
functions of the components.

=head1 METHODS

=over 4

=item * new(nr=>$nr, nh=>$nh, smallE=>$smallE)

Initializes the structure.

$nr is a Photonic::LE::NR2::AllH structure (required).

$nh is the maximum number of Haydock coefficients to use (required).

$smallE is the criteria of convergence for the continued fraction
(defaults to 1e-7)

=item * evaluate($epsA, $epsB)

Returns the macroscopic dielectric function for a given value of the
dielectric functions of the host $epsA and the particle $epsB.

=back

=head1 ACCESORS (read only)

=over 4

=item * nr

The LE::NR2::AllH structure

=item * epsA epsB

The dielectric functions of component A and component B used in the
last calculation.

=item * u

The spectral variable used in the last calculation

=item * epsL

The longitudinal macroscopic function obtained in the last calculation.

=item * nh

The maximum number of Haydock coefficients to use.

=item * nhActual

The actual number of Haydock coefficients used in the last calculation

=item * converged

Flags that the last calculation converged before using up all coefficients

=item * smallE

Criteria of convergence for continued fraction. 0 means don't
check. From Photonic::Roles::EpsParams

=back

=begin Pod::Coverage

=head2 BUILD

=end Pod::Coverage

=cut

package Photonic::LE::NR2::EpsL;
$Photonic::LE::NR2::EpsL::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::LE::NR2::AllH;
use Photonic::Types;
use Photonic::Utils qw(lentzCF);

use List::Util qw(min);

use Moose;
use MooseX::StrictConstructor;

has 'epsA'=>(is=>'ro', isa=>'PDL::Complex', init_arg=>undef, writer=>'_epsA',
    documentation=>'Dielectric function of host');
has 'epsB'=>(is=>'ro', isa=>'PDL::Complex', init_arg=>undef, writer=>'_epsB',
        documentation=>'Dielectric function of inclusions');
has 'u'=>(is=>'ro', isa=>'PDL::Complex', init_arg=>undef, writer=>'_u',
    documentation=>'Spectral variable');
with 'Photonic::Roles::EpsL';

sub evaluate {
    my $self=shift;
    $self->_epsA(my $epsA=shift);
    $self->_epsB(my $epsB=shift);
    $self->_u(my $u=1/(1-$epsB/$epsA));
    my $as=$self->nr->as;
    my $b2s=$self->nr->b2s;
    my $min= min($self->nh, $self->nr->iteration);
    my ($fn, $n)=lentzCF([map {$u-$_} @$as], [map {-$_} @$b2s],
			 $min, $self->smallE);
    # Check this logic:
    my $converged=$n<$min || $self->nr->iteration<=$self->nh;
    $self->_converged($converged);
    $self->_nhActual($n);
    $self->_epsL($epsA*$fn/$u);
    return $self->epsL;
}

__PACKAGE__->meta->make_immutable;

1;
