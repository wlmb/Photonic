=head1 NAME

Photonic::Retarded::GreenF

=head1 VERSION

version 0.009

=head1 SYNOPSIS

   use Photonic::Retarded::GreenF;
   my $G=Photonic::Retarded::GreenF->new(metric=>$m, nh=>$nh);
   my $GreenTensor=$G->evaluate($epsB);

=head1 DESCRIPTION

Calculates the asymetric part of the retarded green's tensor for a given fixed
Photonic::Retarded::Metric structure as a function of the dielectric
functions of the components.

=head1 METHODS

=over 4

=item * new(metric=>$m, nh=>$nh, smallH=>$smallH, smallE=>$smallE,
keepStates=>$k)  

Initializes the structure.

$m Photonic::Retarded::Metric describing the structure and some parametres.

$nh is the maximum number of Haydock coefficients to use.

$smallH and $smallE are the criteria of convergence (default 1e-7) for
Haydock coefficients and continued fraction

$k is a flag to keep states in Haydock calculations (default 0)

=item * evaluate($epsB)

Returns the macroscopic Green's operator for a given value of the
dielectric functions of the particle $epsB. The host's
response $epsA is taken from the metric.  

=back

=head1 ACCESORS (read only)

=over 4

=item * keepStates

Value of flag to keep Haydock states

=item * epsA

Dielectric function of component A

=item * epsB

Dielectric function of componente B

=item * u 

Spectral variable

=item * haydock

Array of Photonic::Retarded::AllH structures, one for each polarization

=item * greenP

Array of Photonic::Retarded::GreenP structures, one for each direction.

=item * greenTensor

The Green's tensor of the last evaluation

=item * nh

The maximum number of Haydock coefficients to use.

=item * nhActual

The actual number of Haydock coefficients used in the last calculation

=item * converged

Flags that the last calculation converged before using up all coefficients

=item * smallH, smallE

Criteria of convergence of Haydock coefficients and continued
fraction. 0 means don't check. 

=back

=cut

package Photonic::Retarded::GreenF;
$Photonic::Retarded::GreenF::VERSION = '0.009';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use PDL::MatrixOps;
use Storable qw(dclone);
use PDL::IO::Storable;
use Photonic::Retarded::AllH;
use Photonic::Retarded::GreenP;
use Photonic::Retarded::Green;
use Moose;
use Photonic::Types;


extends 'Photonic::Retarded::Green';

has 'Chaydock' =>(is=>'ro', isa=>'ArrayRef[Photonic::Retarded::AllH]',
            init_arg=>undef, lazy=>1, builder=>'_build_Chaydock',
            documentation=>'Array of Haydock calculators');
has 'CgreenP'=>(is=>'ro', isa=>'ArrayRef[Photonic::Retarded::GreenP]',
             init_arg=>undef, lazy=>1, builder=>'_build_CgreenP',
             documentation=>'Array of projected G calculators');


after 'evaluate' => sub {
    my $self=shift;
    my @CCgreenP; #array of Green's projections along complex directions.
    my $converged=1;
    my $epsB=$self->epsB;
    foreach(@{$self->CgreenP}){
	push @CCgreenP, $_->evaluate($epsB);
	$converged &&=$_->converged;
    }
    $self->_converged($converged);
    my $nd=$self->geometry->B->ndims;
    my $greenTensor=$self->greenTensor;
    my $m=0;
    for my $i(0..$nd-2){
	for my $j($i+1..$nd-1){
	    $greenTensor->(:,($i),($j))+= -i*($CCgreenP[$m] 
	       -(1/2)*$greenTensor->(:,($i),($i))
		+(1/2)*$greenTensor->(:,($j),($j)) );
	    $greenTensor->(:,($j),($i))+= i*( $CCgreenP[$m]
	       -(1/2)*$greenTensor->(:,($i),($i))
		+(1/2)*$greenTensor->(:,($j),($j)) );
	    $m++
	}
     }
    $self->_greenTensor($greenTensor);
    return $greenTensor;
};


sub _build_Chaydock { # One Haydock coefficients calculator per direction0
    my $self=shift;
    my @Chaydock;
    # This must change if G is not symmetric
    foreach(@{$self->geometry->CunitPairs}){
	my $m=dclone($self->metric); #clone metric, to be safe
	my $e=$_; #polarization
	#Build a corresponding Photonic::Retarded::AllH structure
	my $cchaydock=Photonic::Retarded::AllH->new(
	    metric=>$m, polarization=>$e, nh=>$self->nh,
	    keepStates=>$self->keepStates, smallH=>$self->smallH);
	push @Chaydock, $cchaydock;
    }
    return [@Chaydock]
}

sub _build_CgreenP {
    my $self=shift;
    my @CgreenP;
    foreach(@{$self->Chaydock}){
	my $g=Photonic::Retarded::GreenP->new(
	    haydock=>$_, nh=>$self->nh, smallE=>$self->smallE);
	push @CgreenP, $g;
    }
    return [@CgreenP]
}



__PACKAGE__->meta->make_immutable;
    
1;

__END__
