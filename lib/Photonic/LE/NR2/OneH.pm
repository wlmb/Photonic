=head1 NAME

Photonic::LE::NR2::OneH;

=head1 VERSION

version 0.010

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
 
=back

=cut

package Photonic::LE::NR2::OneH;
$Photonic::OneH::LE::NR2::OneH::VERSION = '0.010';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::FFTW3;
use PDL::Complex;
use List::Util;
use Carp;
use Moose;
use Photonic::Types;
use Photonic::Utils qw(HProd);
has 'geometry'=>(is=>'ro', isa => 'Photonic::Types::GeometryG0',
    handles=>[qw(B dims ndims r G GNorm L scale f)],required=>1);
has 'smallH'=>(is=>'ro', isa=>'Num', required=>1, default=>1e-7,
    	    documentation=>'Convergence criterium for Haydock coefficients');

with 'Photonic::Roles::OneH';

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
    # ri=real or imaginary, ij=cartesian
    #state is c:nx:ny:... gnorm=ij:nx:ny...
    #Have to get cartesian out of the way, thread over it and iterate
    #over the rest 
    my $Gpsi_G=$psi_G*$self->GNorm->mv(0,-1); #^G |psi>
    #the result is ri:nx:ny...:ij
    #Take inverse Fourier transform over all space dimensions,
    #thread over cartesian indices
    #Notice that (i)fftn wants a real 2,nx,ny... piddle, not a complex
    #one. Thus, I have to convert complex to real and back here and
    #downwards. 
    my $Gpsi_R=ifftn($Gpsi_G->real, $self->ndims); #real space ^G|psi>
    #the result is RorI, nx, ny,... cartesian
    #Multiply by characteristic function. Thread cartesian
    my $BGpsi_R=Cscale($Gpsi_R, $self->B); #B^G|psi> in Real Space
    #the result is RorI, nx, ny,... cartesian
    #Transform to reciprocal space
    my $BGpsi_G=fftn($BGpsi_R, $self->ndims); #reciprocal space B^G|psi>
    #the result is RorI, nx, ny,... cartesian
    #Scalar product with Gnorm
    my $GBGpsi_G=Cscale($BGpsi_G, $self->GNorm->mv(0,-1)) #^GB^G|psi>
	# RorI, nx, ny,... cartesian
	# Move cartesian to front and sum over
	->mv(-1,0)->sumover->complex; #^G.B^G|psi>
    # Result is RorI, nx, ny,...
    #Normalization should have been taken care of by fftw3
    return $GBGpsi_G;
}

sub innerProduct {
    return HProd($_[1], $_[2]); #skip self, Hermitian product
}

sub more {
    my $self=shift;
    my $next_b2=shift;
    return $next_b2->re > $self->smallH;
}

sub coerce {
    my $self=shift;
    my $x=shift;
    return $x->re;
}

__PACKAGE__->meta->make_immutable;
    
1;
