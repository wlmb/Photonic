=head1 NAME

Photonic::LE::NR2::OneH;

=head1 VERSION

version 0.011

=head1 SYNOPSIS

    use Photonic::LE::NR2::OneH;
    my $nr=Photonic::LE::NR2::OneH->new(geometry=>$geometry);
    $nr->iterate;
    say $nr->iteration;
    say $nr->current_a;
    say $nr->next_b2;
    my $state=$nr->nextState;

=head1 DESCRIPTION

Implements calculation of Haydock coefficients and Haydock states for
the calculation of the non retarded dielectric function of arbitrary
periodic two component systems in arbitrary number of dimentions. One
Haydock coefficient at a time.

=head1 METHODS

=over 4

=item * new(geometry=>$g[, smallH=>$s])

Create a new Ph::LE::NR2::OneH object with GeometryG0 $g and optional
smallness parameter  $s.

=back

=head1 ACCESORS (read only)

=over 4

=item * geometry Photonic::Types::GeometryG0 

A Photonic::Geometry object defining the geometry of the system,
the charateristic function and the direction of the G=0 vector. Should
be given in the initializer.

=item * B ndims dims r G GNorm L scale f

Accesors handled by geometry (see Photonic::Roles::Geometry)

=item * smallH

A small number used as tolerance to end the iteration. Small negative
b^2 coefficients are taken to be zero. Handled by Photonic::Roles::EpsParams

=item * previousState currentState nextState 

The n-1-th, n-th and n+1-th Haydock states; a complex number for each pixel

=item * current_a

The n-th Haydock coefficient a

=item * current_b2 next_b2 current_b next_b

The n-th and n+1-th b^2 and b Haydock coefficients

=item * iteration

Number of completed iterations

=back

=head1 METHODS

=over 4

=item * iterate

Performs a single Haydock iteration and updates current_a, next_b,
next_b2, next_state, shifting the current values where necessary. Returns 
0 when unable to continue iterating. 

=item * $s= _firstState($self)

Returns the fisrt state $v.

=item * $s=applyOperator($self, $psi_G)

Apply the Hamiltonian operator to state. State is ri:nx:ny... gnorm=i:nx:ny...

=item * $s=innerProduct($self, $left, $right)

Returns the inner product (Hamiltonian product) between states.

=item * $s=magnitude($self, $psi)

Returns the magnitude of a state gotten by taking the square root of the inner product of the state with itself, $self->innerProduct($psi, $psi)->abs->sqrt;.

=item * $c=changesign

Change sign to 

 
=back

=cut

package Photonic::LE::NR2::OneH;
$Photonic::LE::NR2::OneH::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::FFTW3;
use PDL::Complex;
use List::Util;
use Carp;
use Photonic::Types;
use Photonic::Utils qw(HProd);
use Moose;
use MooseX::StrictConstructor;

has 'geometry'=>(is=>'ro', isa => 'Photonic::Types::GeometryG0',
    handles=>[qw(B dims ndims r G GNorm L scale f)],required=>1);
has 'complexCoeffs'=>(is=>'ro', init_arg=>undef, default=>0,
		      documentation=>'Haydock coefficients are real');
with 'Photonic::Roles::OneH', 'Photonic::Roles::UseMask';

sub _firstState { #\delta_{G0}
    my $self=shift;
    my $v=PDL->zeroes(2,@{$self->dims})->complex; #RorI, nx, ny...
    my $arg="(0)" . ",(0)" x $self->ndims; #(0),(0),... ndims+1 times
    $v->slice($arg).=1; #delta_{G0}
    return $v;
}

sub applyOperator { 
    my $self=shift;
    my $psi_G=shift;
    my $mask=undef;
    $mask=$self->mask if $self->use_mask;
    # ri:nx:ny
    #state is ri:nx:ny:... gnorm=i:nx:ny...
    #Have to get cartesian out of the way, thread over it and iterate
    #over the rest 
    my $Gpsi_G=$psi_G*$self->GNorm->mv(0,-1); #^G |psi>
    #Gpsi_G is ri:nx:ny...:i
    #Take inverse Fourier transform over all space dimensions,
    #thread over cartesian indices
    my $Gpsi_R=ifftn($Gpsi_G->real, $self->ndims)->complex; #real space ^G|psi>
    #Gpsi_R is ri:nx:ny:...:i
    #Multiply by characteristic function. Thread cartesian
    my $BGpsi_R=$Gpsi_R*$self->B; #B^G|psi> in Real Space
    #BGpsi_R is ri:nx:ny:...:i
    #Transform to reciprocal space
    my $BGpsi_G=fftn($BGpsi_R->real, $self->ndims)->complex; #<G|B^G|psi>
    #BGpsi_G is ri:nx:ny:...:i
    #Scalar product with Gnorm
    my $GBGpsi_G=($BGpsi_G*$self->GNorm->mv(0,-1)) #^GB^G|psi>
	# ri:nx:ny:...i
	# Move cartesian to front and sum over
	->mv(-1,1)->sumover; #^G.B^G|psi>
    # Result is ri:nx:ny,...
    #Normalization should have been taken care of by fftw3
    $GBGpsi_G=$GBGpsi_G*$mask if defined $mask;
    return $GBGpsi_G;
}

sub innerProduct {
    return HProd($_[1], $_[2]); #skip self, Hermitian product
}

sub magnitude { #magnitude of a state
    return  sqrt($_[1]->Cabs2->sum);
    #could be innerProduct($_[1], $_[1]);
}

sub changesign { #don't change sign
    return 0;
}

__PACKAGE__->meta->make_immutable;
    
1;
