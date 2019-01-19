=head1 NAME

Photonic::Roles::OneH

=head1 VERSION

version 0.010

=head1 SYNOPSIS

    use Photonic::OneH::NR2;
    my $nr=Photonic::OneH::NR2->new(geometry=>$geometry);
    $nr->iterate;
    say $nr->iteration;
    say $nr->current_a;
    say $nr->next_b2;
    my $state=$nr->nextState;

=over 4

=item (for developers)

    package Photonic::OneH::NR2.pm;
    $Photonic::OneH::NR2::VERSION= '0.010';
    use namespace::autoclean;
    use Mosse;
    has...
    with 'Photonic::Roles::OneH';

=back

=head1 DESCRIPTION

Roles consumed by OneH objects to be used in a Photonic
calculation. See also the specific implementations under
=L<Photonic::OneH>. Calculation of Haydock coefficients and Haydock
states, one Haydock coefficient at a time.

=head1 METHODS

=over 4

=item * new(geometry=>$g[, smallH=>$s])

Create a new Ph::OneH::NR2 object with GeometryG0 $g and optional
smallness parameter  $s.

=back

=head1 ACCESORS (read only)

=over 4

=item * geometry Photonic::Types::GeometryG0 

A Photonic::Geometry object defining the geometry of the system,
the charateristic function and the direction of the G=0 vector. Should
be given in the initializer.

=item * B dims r G GNorm L scale f

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

package Photonic::Roles::OneH;
$Photonic::Roles::OneH::VERSION = '0.010';
use Moose::Role;

use PDL::Lite;
use PDL::NiceSlice;
use PDL::FFTW3;
use PDL::Complex;
use List::Util;
use Carp;
#use Photonic::Types;
use Moose::Util::TypeConstraints;

requires '_firstState'; #default first state
requires 'applyOperator'; #Apply Hamiltonian to current state
requires 'innerProduct'; #Apply Hamiltonian to current state
requires 'more'; #check if iterations are done early
requires 'coerce'; #if necessary coerce into desired type (ie real)

has 'previousState' =>(is=>'ro', isa=>'PDL::Complex', writer=>'_previousState',
    init_arg=>undef);
has 'currentState' => (is=>'ro', isa=>'PDL::Complex', writer=>'_currentState',
      lazy=>1, init_arg=>undef,  default=>sub {0+i*0});
has 'nextState' =>(is=>'ro', isa=>maybe_type('PDL::Complex'),
		   writer=>'_nextState',  lazy=>1, default=>\&_firstRState);
has 'current_a' => (is=>'ro', writer=>'_current_a',
    init_arg=>undef);
has 'current_b2' => (is=>'ro', writer=>'_current_b2',
    init_arg=>undef);
has 'next_b2' => (is=>'ro', writer=>'_next_b2', init_arg=>undef, default=>0);
has 'current_b' => (is=>'ro', writer=>'_current_b', init_arg=>undef);
has 'next_b' => (is=>'ro', writer=>'_next_b', init_arg=>undef, default=>0);
has 'iteration' =>(is=>'ro', writer=>'_iteration', init_arg=>undef,
                   default=>0);

sub iterate { #single Haydock iteration
    my $self=shift;
    #Note: calculate Current a, next b2, next b, next state
    #Done if there is no next state
    return 0 unless defined $self->nextState;
    $self->_iterate_indeed; 
}

sub _iterate_indeed {
    my $self=shift;
    #Shift and fetch results of previous calculations
    $self->_previousState(my $prevPsi=$self->currentState);
    $self->_currentState(my $psi=$self->nextState);
    $self->_current_b2($self->next_b2);
    $self->_current_b(my $b=$self->next_b);
    my $opPsi=$self->applyOperator($psi); 
    my $current_a=$self->innerProduct($psi, $opPsi);
    my $bpsi=$opPsi-$current_a*$psi-$b*$prevPsi;
    my $next_b2=$self->innerProduct($bpsi, $bpsi);
    my $next_b=sqrt($next_b2);
    my $next_state=undef;
    $next_state=$bpsi/$next_b if $self->more($next_b2);
    #save values
    $self->_current_a($self->coerce($current_a));
    $self->_next_b2($self->coerce($next_b2));
    $self->_next_b($self->coerce($next_b));
    $self->_nextState($next_state);
    $self->_iteration($self->iteration+1); #increment counter
    return 1;
}

sub _firstRState {
    my $self=shift;
    return $self->_firstState;
}

1;
