=head1 NAME

Photonic::Roles::OneHM

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

=over 4

=item (for developers)

    package Photonic::OneH::R2.pm;
    $Photonic::OneH::R2::VERSION= '0.010';
    use namespace::autoclean;
    use Mosse;
    has...
    with 'Photonic::Roles::OneHM';

=back

=head1 DESCRIPTION

Roles consumed by OneH objects to be used in a Photonic
calculation that require a Metric, as in Eqs. 4.29-4.33 of the thesis
of José Samuel Pérez. See also the specific implementations under 
=L<Photonic::OneH>. Calculation of Haydock coefficients and Haydock
states, one Haydock coefficient at a time.

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

package Photonic::Roles::OneHM;
$Photonic::Roles::OneHM::VERSION = '0.010';
use Moose::Role;

use PDL::Lite;
use PDL::NiceSlice;
use PDL::FFTW3;
use PDL::Complex;
use List::Util;
use Carp;
#use Photonic::Types;
use Moose::Util::TypeConstraints;

requires '_initialState'; #default first state
requires 'applyOperator'; #Apply Hamiltonian to state
requires 'applyMetric'; #Apply metric to state
requires 'innerProduct'; #Inner product between states
requires 'more'; #check if iterations are done early
requires 'changesign'; #if necessary changesign of b2
requires 'coerce'; #if necessary coerce into desired type (ie real)

has 'initialState' =>(is=>'ro', isa=>'PDL::Complex', 
    writer=>'_initial', lazy=>1, default=>\&_initialRState);
has 'previousState' =>(is=>'ro', isa=>'PDL::Complex',
    writer=>'_previousState', lazy=>1, init_arg=>undef, 
    default=>sub {0+i*0});
has 'currentState' => (is=>'ro', isa=>'PDL::Complex',
      writer=>'_currentState', 
      lazy=>1, init_arg=>undef,  default=>sub {0+i*0});
has 'nextState' =>(is=>'ro', isa=>maybe_type('PDL::Complex'), 
    writer=>'_nextState', init_arg=>undef);
has 'current_a' => (is=>'ro', writer=>'_current_a',
    init_arg=>undef);
has 'current_b2' => (is=>'ro', writer=>'_current_b2',
    init_arg=>undef);
has 'next_b2' => (is=>'ro', writer=>'_next_b2', init_arg=>undef, default=>0);
has 'current_b' => (is=>'ro', writer=>'_current_b', init_arg=>undef);
has 'next_b' => (is=>'ro', writer=>'_next_b', init_arg=>undef); 
has 'current_c' => (is=>'ro', writer=>'_current_c', init_arg=>undef); 
has 'next_c' => (is=>'ro', writer=>'_next_c', init_arg=>undef, default=>0);
has 'next_bc' => (is=>'ro', writer=>'_next_bc', init_arg=>undef, default=>0);
has 'previous_g' => (is=>'ro', writer=>'_previous_g', init_arg=>undef);
has 'current_g' => (is=>'ro', writer=>'_current_g', init_arg=>undef, 
     default=>0);
has 'next_g' => (is=>'ro', writer=>'_next_g', init_arg=>undef);
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
    #Notation: nm1 is n-1, np1 is n+1
    $self->_previousState(my $psi_nm1=$self->currentState);
    $self->_currentState(my $psi_n #state in reciprocal space 
			 =$self->nextState);
    $self->_current_b(my $b_n=$self->next_b);
    $self->_current_b2(my $b2_n=$self->next_b2);
    $self->_current_c(my $c_n=$self->next_c);
    $self->_previous_g(my $g_nm1=$self->current_g);
    $self->_current_g(my $g_n=$self->next_g);
    my $gPsi=$self->applyMetric($psi_n);
    my $psi_np1=$self->applyOperator($gPsi);
    # Eq. 4.41
    my $a_n=$g_n*$self->innerProduct($psi_n, $psi_np1);
    # Eq. 4.30
    #Modified algorithm:
    my $next_state=($psi_np1-$a_n*$psi_n-$c_n*$psi_nm1);
    my $b2_np1=$self->innerProduct($next_state,$next_state);
    my $g_np1=1;
    $g_np1=-1, $b2_np1=-$b2_np1 if $self->changesign($b2_np1);
    my $b_np1=sqrt($b2_np1);
    # Eq. 4.31
    my $c_np1=$g_np1*$g_n*$b_np1;
    my $bc_np1=$g_np1*$g_n*$b2_np1;
    # Eq. 4.33
    $next_state=undef unless $self->more($b2_np1);
    $next_state=$next_state/$b_np1 if defined $next_state;
    #save values
    $self->_nextState($next_state);
    $self->_current_a($self->coerce($a_n));
    $self->_next_b2($self->coerce($b2_np1));
    $self->_next_b($self->coerce($b_np1));
    $self->_next_g($g_np1);
    $self->_next_c($self->coerce($c_np1));
    $self->_next_bc($self->coerce($bc_np1));
    $self->_iteration($self->iteration+1); #increment counter
    return 1;
}

sub _initialRState {
    my $self=shift;
    return $self->_initialState;
}

sub BUILD {} #noop
after BUILD => sub {
    my $self=shift;
    my $phi=$self->initialState;
    my $g=1;
    my $b2=$self->innerProduct($phi,$phi);
    $g=-1, $b2=-$b2 if $self->changesign($b2);
    $b=sqrt($b2);
    $phi=$phi/$b; #first state;
    $self->_nextState($phi);
    #skip $self->current_a;
    $self->_next_b2($b2);
    $self->_next_b($b);
    $self->_next_c(0); #no c0
    $self->_next_bc(0); #no bc0
    $self->_next_g($g);
};

1;
