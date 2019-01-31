package Photonic::WE::R2::OneH;
$Photonic::WE::R2::OneH::VERSION = '0.010';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::FFTW3;
use PDL::Complex;
use List::Util;
use Carp;
use Moose;
#use Photonic::Types;
#with 'Photonic::Roles::EpsParams';
use Photonic::Utils qw(HProd);

with 'Photonic::Roles::OneHM';

has 'metric'=>(is=>'ro', isa => 'Photonic::WE::R2::Metric',
    handles=>[qw(B ndims dims epsilon)],required=>1);
has 'polarization' =>(is=>'ro', required=>1, isa=>'PDL::Complex');
has 'smallH'=>(is=>'ro', isa=>'Num', required=>1, default=>1e-7,
    	    documentation=>'Convergence criterium for Haydock coefficients');
has 'normalizedPolarization' =>(is=>'ro', isa=>'PDL::Complex',
     init_arg=>undef, writer=>'_normalizedPolarization');


sub applyOperator {
    my $self=shift;
    my $psi=shift;
    # psi is RorI xyz nx ny nz. Get cartesian out of the way and
    # transform to real space. Note FFFTW3 wants real PDL's[2,...] 
    my $psi_r=ifftn($psi->real->mv(1,-1), $self->ndims);
    #$psi_r is RorI nx ny nz  xyz, B is nx ny nz
    # Multiply by characteristic function
    my $Bpsi_r=Cscale($psi_r,$self->B);
    #Bpsi_r is RorI nx ny nz  xyz
    #Transform to reciprocal space, move xyz back and make complex, 
    my $psi_G=fftn($Bpsi_r, $self->ndims)->mv(-1,1)->complex;
    return $psi_G;
}

sub applyMetric {
    my $self=shift;
    my $psi=shift;
    #psi is RorI xy.. nx ny..
    my $g=$self->metric->value;
    #$g is xyz xyz nx ny nz 
    my $gpsi=($g*$psi(:,:,*1))->sumover; #matrix times vector
    #$gpsi is RorI xy.. nx ny..
    return $gpsi;
}

sub innerProduct {  #ignore $self. Return Hermitian product
    return HProd($_[1], $_[2]);
}

sub more { #check if I should continue
    my $self=shift;
    my $b2=shift;
    return $b2->re > $self->smallH;
}

sub changesign { #flag change sign required if b^2 negative
    return $_[1]->re < 0;
}

sub coerce { #Ignore $self. Take real part
    return $_[1]->re;
}

    



sub _initialState {
    my $self=shift;
    my $d=$self->ndims;
    my $v=PDL->zeroes(@{$self->dims}); #build a nx ny nz pdl
    my $arg="(0)" . ",(0)" x ($d-1); #(0),(0),... ndims times
    $v->slice($arg).=1; #delta_{G0}
    my $e=$self->polarization; #RorI xyz
    croak "Polarization has wrong dimensions. " .
	  " Should be $d-dimensional complex vector."
	unless $e->isa('PDL::Complex') && $e->ndims==2 &&
	[$e->dims]->[0]==2 && [$e->dims]->[1]==$d; 
    my $modulus2=$e->Cabs2->sumover;
    croak "Polarization should be non null" unless
	$modulus2 > 0;
    $e=$e/sqrt($modulus2);
    $self->_normalizedPolarization($e);
    my $phi=$e*$v(*1); #initial state ordinarily normalized 
                       # RorI xyz nx ny nz
    return $phi;
}

    
__PACKAGE__->meta->make_immutable;
    
1;

=head1 NAME

Photonic::OneH::R2

=head1 VERSION

version 0.010

=head1 SYNOPSIS

    use Photonic::OneH::R2;
    my $nr=Photonic::OneH::R2->new(metric=>$g, polarization=>$p);
    $nr->iterate;
    say $nr->iteration;
    say $nr->current_a;
    say $nr->next_b2;
    my $state=$nr->nextState;

=head1 DESCRIPTION

Implements calculation of Haydock coefficients and Haydock states for
the calculation of the retarded dielectric function of arbitrary
periodic two component systems in arbitrary number of dimentions. One
Haydock coefficient at a time.

=head1 METHODS

=over 4

=item * new(metric=>$m, polarization=>$e, [, smallH=>$s])

Create a new Ph::OneH::R2 object with PDL::Metric::R2 $m, with a
field along the complex direction $e and with smallness parameter  $s.

=back

=head1 ACCESORS (read only)

=over 4

=item * metric Photonic::Metric::R2 

A Photonic::Metric::R2 object defining the geometry of the
system, the charateristic function, the wavenumber, wavevector and
host dielectric function. Required in the initializer.

=item * polarization PDL::Complex

A non null vector defining the complex direction of the macroscopic
field. 

=item * smallH 

A small number used as tolerance to end the iteration. Small negative
b^2 coefficients are taken to be zero. From Photonic::Roles::EpsParams.

=item * B ndims dims epsilon

Accesors handled by metric (see Photonic::Metric::R2)

=item * previousState currentState nextState 

The n-1-th, n-th and n+1-th Haydock states; a complex vector for each
reciprocal wavevector

=item * current_a

The n-th Haydock coefficient a

=item * current_b2 next_b2 current_b next_b

The n-th and n+1-th b^2 and b Haydock coefficients

=item * next_c

The n+1-th c Haydock coefficient

=item * previous_g current_g next_g 

The n-1-th n-th and n+1-th g Haydock coefficients

=item * iteration

Number of completed iterations

=back

=head1 METHODS

=over 4

=item * iterate

Performs a single Haydock iteration and updates current_a, next_b,
next_b2, next_c, next_g, next_state, shifting the current values where
necessary. Returns 0 when unable to continue iterating. 
 
=begin Pod::Coverage

=head2 BUILD

=end Pod::Coverage

=back

=cut

