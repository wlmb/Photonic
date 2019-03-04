=head1 NAME

Photonic::LE::S::OneH

=head1 VERSION

version 0.010

=head1 SYNOPSIS

    use Photonic::LE::S::OneH;
    my $nr=Photonic::LE::S::OneH->new(epsilon=>$epsilon,
           geometry=>$geometry);  
    $nr->iterate;
    say $nr->iteration;
    say $nr->current_a;
    say $nr->next_b2;
    my $state=$nr->nextState;

=head1 DESCRIPTION

Implements calculation of Haydock coefficients and Haydock states for
the calculation of the non retarded dielectric function of arbitrary
periodic N component systems in arbitrary number of dimensions. One
Haydock coefficient at a time. Use k,-k spinors. MQ notes.

=head1 METHODS

=over 4

=item * new(epsilon=>$e, geometry=>$g[, smallH=>$s])

Create a new Photonic::LE::S::OneH object with GeometryG0 $g, dielectric
function $e and optional smallness parameter  $s.

=back

=head1 ACCESORS (read only)

=over 4

=item * epsilon 

A PDL::Complex PDL giving the value of the dielectric function epsilon
for each pixel of the system

=item * geometry Photonic::Types::GeometryG0 

A Photonic::Geometry object defining the geometry of the system,
the charateristic function and the direction of the G=0 vector. Should
be given in the initializer.

=item * B ndims dims r G GNorm L scale f

Accesors handled by geometry (see Photonic::Geometry)

=item * smallH

A small number used as tolerance to end the iteration. Small negative
b^2 coefficients are taken to be zero. Handled by Photonic::Roles::EpsParams

=item * previousState currentState nextState 

The n-1-th, n-th and n+1-th Haydock states; a complex 2-spinor for each
reciprocal vector.

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

Performs a single Haydock iteration and updates current_a, next_state,
next_b2, next_b, shifting the current values where necessary. Returns
0 when unable to continue iterating. 
 
=back

=cut

package Photonic::LE::S::OneH;
$Photonic::LE::S::OneH::VERSION = '0.010';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::FFTW3;
use PDL::Complex;
use List::Util;
use Carp;
use Moose;
use Photonic::Types;
use Photonic::Utils qw(SProd);

has 'epsilon'=>(is=>'ro', isa=>'PDL::Complex', required=>1, lazy=>1,
		builder=>'_epsilon');  
has 'geometry'=>(is=>'ro', isa => 'Photonic::Types::GeometryG0',
    handles=>[qw(B ndims dims r G GNorm L scale f pmGNorm)],required=>1
    );
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

sub _firstState { #\delta_{G0}
    my $self=shift;
    my $v=PDL->zeroes(2,2,@{$self->dims})->complex; #RorI,pmk, nx, ny...
    my $arg="(0),:" . ",(0)" x $self->ndims; #(0),(0),... ndims+1 times
    $v->slice($arg).=1/sqrt(2); 
    return $v;
}
sub applyOperator {
    my $self=shift;
    my $psi_G=shift;
    confess "State should be complex" unless $psi_G->isa('PDL::Complex');
    #Each state is a spinor with two wavefunctions \psi_{k,G} and
    #\psi_{-k,G}, thus the index plus or minus k, pmk.
    #Notation ri=real or imaginary, pmk=+ or - k, i,j=cartesian
    #state is ri:pmk:nx:ny... pmGnorm=i:pmk:nx:ny...
    #Multiply by vectors ^G and ^(-G).
    #Have to get cartesian out of the way, thread over it and iterate
    #over the rest 
    my $Gpsi_G=($psi_G*$self->pmGNorm->mv(0,-1))->mv(1,-1); #^G |psi>
    #the result is complex ri:nx:ny...i:pmk
    # Notice that I actually multiply by unit(k-G) instead of
    # unit(-k+G) when I use pmGNorm; as I do it two times, the result
    # is the same.
    #Take inverse Fourier transform over all space dimensions,
    #thread over cartesian and pmk indices
    #Notice that (i)fftn wants a real 2,nx,ny... piddle, not a complex
    #one. Thus, I have to convert complex to real and back here and
    #downwards. 
    #real space ^G|psi> 
    my $Gpsi_R=ifftn($Gpsi_G->real, $self->ndims)->complex; 
    # $Gpsi_R is ri:nx:ny...i:pmk
    # $self->epsilon is ri:nx:ny...
    #Multiply by the dielectric function in Real Space. Thread
    #cartesian and pm indices
    my $eGpsi_R=$self->epsilon*$Gpsi_R; #Epsilon could be tensorial!
    #$eGpsi_R is ri:nx:ny,...i:pmk
    #Transform to reciprocal space
    my $eGpsi_G=fftn($eGpsi_R->real, $self->ndims)->complex->mv(-1,1); 
    #$eGpsi_G is ri:pmk:nx:ny...:i
    #Scalar product with pmGnorm: i:pm:nx:ny...
    my $GeGpsi_G=($eGpsi_G*$self->pmGNorm->mv(0,-1)) #^Ge^G|psi>
	# ri:pmk:nx:ny...:i
	# Move cartesian to front and sum over
	->mv(-1,1)->sumover; #^G.epsilon^G|psi>
    return $GeGpsi_G;
}
sub innerProduct {
    #ignore self
    return SProd($_[1], $_[2]);
}
sub magnitude {
    my $self=shift;
    my $psi=shift;
    return $self->innerProduct($psi, $psi)->abs->sqrt;
}
sub changesign { #don't change sign
    return 0;
}

__PACKAGE__->meta->make_immutable;
    
1;
