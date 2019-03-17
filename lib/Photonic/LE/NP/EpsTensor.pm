=head1 NAME

Photonic::LE::NP::EpsTensor

=head1 VERSION

version 0.011

=head1 SYNOPSIS

   use Photonic::LE::NP::EpsTensor;
   my $eps=Photonic::LE::NP::EpsTensor->new(
                     epsilon=>$e, geometry=>$g);
   my $epsilonTensor=$epsTensor;

=head1 DESCRIPTION

Calculates the dielectric tensor for a given fixed
Photonic::Geometry structure as a function of the dielectric
functions of the components.

=head1 METHODS

=over 4

=item * new(epsilon=>$e, geometry=>$g, nh=>$nh, smallH=>$smallH, 
            smallE=>$smallE, keepStates=>$k) 

Initializes the structure.

$e PDL::Complex is the dielectric function as a complex scalar field
 
$g Photonic::Geometry describing the structure

$nh is the maximum number of Haydock coefficients to use.

$smallH and $smallE are the criteria of convergence (default 1e-7) for
the Haydock coefficients and the tensor calculations.

$k is a flag to keep states in Haydock calculations (default 0)

=back

=head1 ACCESORS (read only)

=over 4

=item * epsilon 

A PDL::Complex PDL giving the value of the dielectric function epsilon
for each pixel of the system

=item * keepStates

Value of flag to keep Haydock states

=item * nr

Array of Photonic::LE::NP::AllH structures, one for each direction

=item * epsL

Array of Photonic::LE::NP::EpsL structures, one for each direction.

=item * epsTensor

The valuated dielectric tensor 

=item * nh

The maximum number of Haydock coefficients to use.

=item * converged

Flags that the last calculation converged before using up all coefficients

=item * smallH smallE

Criteria of convergence for Haydock and epsilon calculations. 0 means
don't check. From Photonic::Roles::EpsParams.

    *Check last remark* 

=back

=cut

package Photonic::LE::NP::EpsTensor;
$Photonic::LE::NP::EpsTensor::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use PDL::MatrixOps;
use Storable qw(dclone);
use PDL::IO::Storable;
use Photonic::LE::NP::AllH;
use Photonic::LE::NP::EpsL;
use Moose;
use Photonic::Types;
with 'Photonic::Roles::EpsParams';

has 'epsilon'=>(is=>'ro', isa=>'PDL::Complex', required=>1);
has 'geometry'=>(is=>'ro', isa => 'Photonic::Geometry',
    handles=>[qw(B dims r G GNorm L scale f)],required=>1
);
with 'Photonic::Roles::KeepStates';
with 'Photonic::Roles::EpsParams';
has 'reorthogonalize'=>(is=>'ro', required=>1, default=>0,
         documentation=>'Reorthogonalize haydock flag');
has 'nr' =>(is=>'ro', isa=>'ArrayRef[Photonic::LE::NP::AllH]',
            init_arg=>undef, lazy=>1, builder=>'_build_nr',
            documentation=>'Array of Haydock calculators');
has 'epsL'=>(is=>'ro', isa=>'ArrayRef[Photonic::LE::NP::EpsL]',
             init_arg=>undef, lazy=>1, builder=>'_build_epsL',
             documentation=>'Array of epsilon calculators');
has 'epsTensor'=>(is=>'ro', isa=>'PDL', init_arg=>undef, lazy=>1,
		  builder=>'_build_epsTensor',
		  documentation=>'Dielectric Tensor');
has 'converged'=>(is=>'ro', init_arg=>undef, writer=>'_converged',
             documentation=>
                  'All EpsL evaluations converged in last evaluation'); 

sub _build_epsTensor {
    my $self=shift;
    my @eps; #array of @eps along different directions.
    my $converged=1;
    foreach(@{$self->epsL}){
	push @eps, $_->epsL;
	$converged &&=$_->converged;
    }
    $self->_converged($converged);
    my $reEpsL=PDL->pdl([map {$_->re} @eps]);
    my $imEpsL=PDL->pdl([map {$_->im} @eps]);
    my ($lu, $perm, $parity)=@{$self->geometry->unitDyadsLU};
    my $reEps=lu_backsub($lu, $perm, $parity, $reEpsL);
    my $imEps=lu_backsub($lu, $perm, $parity, $imEpsL);
    my $nd=$self->geometry->B->ndims;
    my $epsTensor=PDL->zeroes(2, $nd, $nd)->complex;
    my $n=0;
    for my $i(0..$nd-1){
	for my $j($i..$nd-1){
	    $epsTensor->(:,($i),($j)).=$reEps->($n)+i*$imEps->($n);
	    $epsTensor->(:,($j),($i)).=$reEps->($n)+i*$imEps->($n);
	    ++$n;
	}
    }
    return $epsTensor;
}

sub _build_nr { # One Haydock coefficients calculator per direction0
    my $self=shift;
    my @nr;
    foreach(@{$self->geometry->unitPairs}){
	my $g=dclone($self->geometry); #clone geometry
	$g->Direction0($_); #add G0 direction
	#Build a corresponding LE::NP::AllH structure
	my $nr=Photonic::LE::NP::AllH->new(
	    epsilon=>$self->epsilon, geometry=>$g, smallH=>$self->smallH, 
	    nh=>$self->nh, keepStates=>$self->keepStates,
	    reorthogonalize=>$self->reorthogonalize);
	push @nr, $nr;
    }
    return [@nr]
}

sub _build_epsL {
    my $self=shift;
    my @eps;
    foreach(@{$self->nr}){
	my $e=Photonic::LE::NP::EpsL->
	    new(nr=>$_, nh=>$self->nh, smallE=>$self->smallE);
	push @eps, $e;
    }
    return [@eps]
}


__PACKAGE__->meta->make_immutable;
    
1;

__END__
