=head1 NAME

Photonic::NonRetarded::NPS::OneH

=head1 VERSION

version 0.010

=head1 SYNOPSIS

    use Photonic::NonRetarded::NPS::OneH;
    my $nr=Photonic::NonRetarded::NPS::OneH->new(epsilon=>$epsilon,
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

Create a new Ph::NR::NPS::OneH object with GeometryG0 $g, dielectric
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

=item * B dims r G GNorm L scale f

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

package Photonic::NonRetarded::NPS::OneH;
$Photonic::NonRetarded::NPS::OneH::VERSION = '0.010';
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
with "Photonic::Roles::EpsParams";

has 'epsilon'=>(is=>'ro', isa=>'PDL::Complex', required=>1);
has 'geometry'=>(is=>'ro', isa => 'Photonic::Types::GeometryG0',
    handles=>[qw(B dims r G GNorm L scale f pmGNorm)],required=>1
);
has 'currentState' => (is=>'ro', isa=>'PDL::Complex', writer=>'_currentState',
      lazy=>1, init_arg=>undef,  default=>sub {0+i*0});
has 'previousState' =>(is=>'ro', isa=>'PDL::Complex', writer=>'_previousState',
    init_arg=>undef);
has 'nextState' =>(is=>'ro', isa=>'PDL::Complex|Undef', writer=>'_nextState',
    lazy=>1, default=>\&_firstState);
has 'current_a' => (is=>'ro', isa=>'PDL::Complex', writer=>'_current_a',
    init_arg=>undef);
has 'current_b2' => (is=>'ro', isa=>'PDL::Complex', writer=>'_current_b2',
    init_arg=>undef);
has 'next_b2' => (is=>'ro', isa=>'PDL::Complex', writer=>'_next_b2', 
		  init_arg=>undef, default=>sub{0+0*i});
has 'current_b' => (is=>'ro', isa=>'PDL::Complex', 
		    writer=>'_current_b', init_arg=>undef);
has 'next_b' => (is=>'ro', isa=>'PDL::Complex', writer=>'_next_b', 
		 init_arg=>undef, default=>sub {0+0*i});
has 'iteration' =>(is=>'ro', writer=>'_iteration', init_arg=>undef,
                   default=>0);
sub iterate { #single Haydock iteration in N=1,2,3 dimensions
    my $self=shift;
    #Note: calculate Current a, next b2, next b, next state
    #Done if there is no next state
    return 0 unless defined $self->nextState;
    $self->_iterate_indeed; 
}
sub _iterate_indeed {
    my $self=shift;
    #Shift and fetch results of previous calculations
    $self->_previousState($self->currentState);
    $self->_currentState( #state in reciprocal space 
			  my $psi_G = $self->nextState);
    $self->_current_b2($self->next_b2);
    $self->_current_b($self->next_b);
    #Each state is a spinor with two wavefunctions \psi_{k,G} and
    #\psi_{-k,G}, thus the index plus or minus k, pmk.
    #state is RorI, pmk, nx, ny... pmGnorm=cartesian,pmk,nx,ny...
    #Multiply by vectors ^G and ^(-G).
    #Have to get cartesian out of the way, thread over it and iterate
    #over the rest 
    my $Gpsi_G=Cscale($psi_G, $self->pmGNorm->mv(0,-1))->mv(1,-1); #^G |psi>
    #the result is complex RorI, nx, ny,... cartesian, pmk
    #Take inverse Fourier transform over all space dimensions,
    #thread over cartesian and pmk indices
    #Notice that (i)fftn wants a real 2,nx,ny... piddle, not a complex
    #one. Thus, I have to convert complex to real and back here and
    #downwards. 
    my $Gpsi_R=ifftn($Gpsi_G->real, $self->B->ndims)->complex; #real
				#space ^G|psi> 
    #the result is RorI, nx, ny,... cartesian
    # Multiply by the dielectric function in Real Space. Thread
    # cartesian and pm indices
    #the result is RorI, nx, ny,... cartesian, pmk
    my $eGpsi_R=$self->epsilon*$Gpsi_R; #Epsilon could be tensorial!
    #Transform to reciprocal space
    my $eGpsi_G=fftn($eGpsi_R->real, $self->B->ndims)
	->complex->mv(-1,1); #reciprocal space B^G|psi> 
    #the result is RorI, pmk, nx, ny,... cartesian
    #Scalar product with pmGnorm
    my $GeGpsi_G=Cscale($eGpsi_G, $self->pmGNorm->mv(0,-1)) #^Ge^G|psi>
	# RorI, pmk, nx, ny,... cartesian
	# Move cartesian to front and sum over
	->mv(-1,1)->sumover; #^G.epsilon^G|psi>
    #Result is ^G.epsilon^G|psi>, RorI, pmk, nx, ny,...
    #Normalization should have been taken care of by fftw3
    #Instead of conjugating the current \psi, change G to -G
    #First reverse all reciprocal dimensions
     # Calculate Haydock coefficients
    # current_a is the Euclidean product with G -> -G
    my $current_a=SProd($psi_G,$GeGpsi_G);
    # next b^2
    my $bpsi_G=$GeGpsi_G - $current_a*$psi_G 
	- $self->current_b*$self->previousState;
    my $next_b2= SProd($bpsi_G,$bpsi_G);
    my $next_b=sqrt($next_b2);
    my $next_state=undef;
    $next_state=$bpsi_G/$next_b if($next_b2->Cabs > $self->smallH);
    #save values
    $self->_current_a($current_a);
    $self->_next_b2($next_b2);
    $self->_next_b($next_b);
    $self->_nextState($next_state);
    $self->_iteration($self->iteration+1); #increment counter
    return 1;
}

sub _firstState { #\delta_{G0}
    my $self=shift;
    my $v=PDL->zeroes(2,2,@{$self->dims})->complex; #RorI,pmk, nx, ny...
    my $arg="(0),:" . ",(0)" x $self->B->ndims; #(0),(0),... ndims+1 times
    $v->slice($arg).=1/sqrt(2); 
    return $v;
}


__PACKAGE__->meta->make_immutable;
    
1;
