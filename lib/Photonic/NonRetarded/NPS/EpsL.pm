=head1 NAME

Photonic::NonRetarded::NPS::EpsL

=head1 VERSION

version 0.010

=head1 SYNOPSIS

   use Photonic::NonRetarded::NPS::EpsL;
   my $eps=Photonic::NonRetarded::NPS::EpsL->new(nr=>$nr, nh=>$nh);
   my $epsilonLongitudinal=$eps->epsL;

=head1 DESCRIPTION

Calculates the longitudinal dielectric function for a given fixed
Photonic::NonRetarded::NPS::AllH structure as a function of the dielectric
functions of the components.

=head1 METHODS

=over 4

=item * new(nr=>$nr, nh=>$nh, smallE=>$smallE)

Initializes the structure.

$nr is a Photonic::NonRetarded::NPS::AllH structure (required).

$nh is the maximum number of Haydock coefficients to use (required).

$smallE is the criteria of convergence for the continued fraction 
(defaults to 1e-7)

=back

=head1 ACCESORS (read only)

=over 4

=item * epsL

The longitudinal macroscopic function.

=item * nr

The NonRetarded::NPS::AllH structure

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

package Photonic::NonRetarded::NPS::EpsL;
$Photonic::NonRetarded::NPS::EpsL::VERSION = '0.010';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::NonRetarded::NPS::AllH;
use Moose;
use Photonic::Types;

with 'Photonic::Roles::EpsParams';
has 'nr' =>(is=>'ro', isa=>'Photonic::NonRetarded::NPS::AllH', required=>1);
has 'epsL'=>(is=>'ro', isa=>'PDL::Complex', init_arg=>undef, writer=>'_epsL');
has 'nhActual'=>(is=>'ro', isa=>'Num', init_arg=>undef, 
                 writer=>'_nhActual');
has 'converged'=>(is=>'ro', isa=>'Num', init_arg=>undef, writer=>'_converged');

sub BUILD {
    my $self=shift;
    $self->nr->run unless $self->nr->iteration;
    my $as=$self->nr->as;
    my $b2s=$self->nr->b2s;
    # Continued fraction evaluation: Lentz method
    # Numerical Recipes p. 171
    my $tiny=1.e-30;
    my $converged=0;
#    b0+a1/b1+a2/...
#	lo debo convertir a
#	a0-b1^2/a1-b2^2/
#	entonces bn->an y an->-b_n^2
    my $fn=$as->[0];
    $fn=r2C($tiny) if $fn->re==0 and $fn->im==0;
    my $n=1;
    my ($fnm1,$Cnm1, $Dnm1)=($fn, $fn, r2C(0)); #previous coeffs.
    my ($fn, $Cn, $Dn); #current coeffs.
    my $Deltan;
    while($n<$self->nh && $n<$self->nr->iteration){
	$Dn=$as->[$n]-$b2s->[$n]*$Dnm1;
	$Dn=r2C($tiny) if $Dn->re==0 and $Dn->im==0;
	$Cn=$as->[$n]-$b2s->[$n]/$Cnm1;
	$Cn=r2C($tiny) if $Cn->re==0 and $Cn->im==0;
	$Dn=1/$Dn;
	$Deltan=$Cn*$Dn;
	$fn=$fnm1*$Deltan;
	last if $converged=$Deltan->approx(1, $self->smallE)->all;
	$fnm1=$fn;
	$Dnm1=$Dn;
	$Cnm1=$Cn;
	$n++;
    }
    #If there are less available coefficients than $self->nh and all
    #of them were used, there is no remaining work to do, so, converged 
    $converged=1 if $self->nr->iteration < $self->nh;
    $self->_converged($converged);
    $self->_nhActual($n);
    $self->_epsL($fn);
}

__PACKAGE__->meta->make_immutable;
    
1;
