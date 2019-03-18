=head1 NAME

Photonic::LE::S::EpsL

=head1 VERSION

version 0.011

=head1 SYNOPSIS

   use Photonic::LE::S::EpsL;
   my $eps=Photonic::LE::S::EpsL->new(nr=>$nr, nh=>$nh);
   my $epsilonLongitudinal=$eps->epsL;

=head1 DESCRIPTION

Calculates the longitudinal dielectric function for a given fixed
Photonic::LE::S::AllH structure as a function of the dielectric
functions of the components.

=head1 METHODS

=over 4

=item * new(nr=>$nr, nh=>$nh, smallE=>$smallE)

Initializes the structure.

$nr is a Photonic::LE::S::AllH structure (required).

$nh is the maximum number of Haydock coefficients to use (required).

$smallE is the criteria of convergence for the continued fraction 
(defaults to 1e-7)

=back

=head1 ACCESORS (read only)

=over 4

=item * epsL

The longitudinal macroscopic function.

=item * nr

The NonRetarded::LE::S::AllH structure

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

package Photonic::LE::S::EpsL;
$Photonic::LE::S::EpsL::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::LE::S::AllH;
use Photonic::Types;
use Photonic::Utils qw(lentzCF);
use List::Util qw(min);
use Moose;
use Moose;
use MooseX::StrictConstructor;

has 'nr' =>(is=>'ro', isa=>'Photonic::LE::S::AllH', required=>1);
has 'epsL'=>(is=>'ro', isa=>'PDL::Complex', init_arg=>undef, writer=>'_epsL');
has 'nhActual'=>(is=>'ro', isa=>'Num', init_arg=>undef, 
                 writer=>'_nhActual');
has 'converged'=>(is=>'ro', isa=>'Num', init_arg=>undef, writer=>'_converged');
with 'Photonic::Roles::EpsL';

after BUILD => sub {
    my $self=shift;
    my $as=$self->nr->as;
    my $b2s=$self->nr->b2s;
    my $min= min($self->nh, $self->nr->iteration);  
    my ($fn, $n)=lentzCF($as, [map {-$_} @$b2s], $min, $self->smallE);  
    # Check this logic:
    my $converged=$n<$min || $self->nr->iteration<=$self->nh;
    $self->_converged($converged);
    $self->_nhActual($n);
    $self->_epsL($fn);
    return $fn;
};

__PACKAGE__->meta->make_immutable;
    
1;
