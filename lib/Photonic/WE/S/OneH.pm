package Photonic::WE::S::OneH;
$Photonic::WE::S::OneH::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::FFTW3;
use PDL::Complex;
use List::Util;
use Carp;
use Moose;
use Photonic::Types;
use Photonic::Utils qw(VSProd);

has 'epsilon'=>(is=>'ro', isa=>'PDL::Complex', required=>1, lazy=>1,
		builder=>'_epsilon');  
has 'metric'=>(is=>'ro', isa => 'Photonic::WE::S::Metric',
	       handles=>{B=>'B', ndims=>'ndims', dims=>'dims',
			 geometry=>'geometry', epsilonR=>'epsilon'},
	       required=>1);
has 'polarization' =>(is=>'ro', required=>1, isa=>'PDL::Complex');
has 'normalizedPolarization' =>(is=>'ro', isa=>'PDL::Complex',
     init_arg=>undef, writer=>'_normalizedPolarization');
has 'complexCoeffs'=>(is=>'ro', init_arg=>undef, default=>1,
		      documentation=>'Haydock coefficients are complex');
with 'Photonic::Roles::OneH';

sub _epsilon {
    my $self=shift;
    die "Coudln't obtain dielectric function from geometry" unless
	$self->geometry->can('epsilon');
    return $self->geometry->epsilon;
}
    
#Required by Photonic::Roles::OneH

sub applyOperator {
    my $self=shift;
    my $psi=shift; #psi is ri:xy:pm:nx:ny
    my $gpsi=$self->applyMetric($psi);
    # gpsi is ri:xy:pm:nx:ny. Get cartesian and pm out of the way and
    # transform to real space. Note FFFTW3 wants real PDL's[2,...] 
    my $gpsi_r=ifftn($gpsi->real->mv(1,-1)->mv(1,-1), $self->ndims)->complex;
    #ri:nx:ny:xy:pm
    my $H=($self->epsilonR-$self->epsilon)/$self->epsilonR;
    my $Hgpsi_r=$H*$gpsi_r; #ri:nx:ny:xy:pm
    #Transform to reciprocal space, move xy and pm back and make complex, 
    my $psi_G=fftn($Hgpsi_r->real, $self->ndims)->mv(-1,1)->mv(-1,1)->complex;
    return $psi_G;
}

sub applyMetric {
    my $self=shift;
    my $psi=shift;
    #psi is ri:xy:pm:nx:ny
    my $g=$self->metric->value;
    #$g is xy:xy:pm:nx:ny  or ri:xy:xy:pm:nx:ny 
    #real or complex matrix times complex vector
    my $gpsi=($g*$psi(:,:,*1)) #ri:xy:xy:pm:nx:ny
	->sumover; #ri:xy:pm:nx:ny
    return $gpsi;
}

sub innerProduct {  #Return Hermitian product with metric
    my $self=shift;
    my $psi1=shift;
    my $psi2=shift;
    my $gpsi2=$self->applyMetric($psi2);
    return VSProd($psi1, $gpsi2);
}


sub magnitude {
    my $self=shift;
    my $psi=shift;
    return $self->innerProduct($psi, $psi)->abs->sqrt;
}
sub changesign { #don't change sign
    return 0;
}

sub _firstState { #\delta_{G0}
    my $self=shift;
    my $v=PDL->zeroes(2,2,@{$self->dims})->complex; #ri:pm:nx:ny...
    my $arg="(0),:" . ",(0)" x $self->ndims; #(0),(0),... ndims+1 times
    $v->slice($arg).=1/sqrt(2); 
    my $e=$self->polarization; #ri:xy
    my $d=[$e->dims]->[1];
    croak "Polarization has wrong dimensions. " .
	  " Should be $d-dimensional complex vector."
	unless $e->isa('PDL::Complex') && $e->ndims==2 &&
	[$e->dims]->[0]==2 && [$e->dims]->[1]==$d; 
    my $modulus2=$e->Cabs2->sumover;
    croak "Polarization should be non null" unless
	$modulus2 > 0;
    $e=$e/sqrt($modulus2);
    $self->_normalizedPolarization($e);
    #I'm using the same polarization for k and for -k. Could be
    #different (for chiral systems, for example
    my $phi=$e*$v(,*1); #initial state ordinarily normalized 
                       # ri:xy:pm:nx:ny
    return $phi;
}


__PACKAGE__->meta->make_immutable;
    
1;
