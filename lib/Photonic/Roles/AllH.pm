=head1 NAME

Photonic::Roles::AllH

=head1 VERSION

version 0.010

=head1 SYNOPSIS

   use Photonic::LE::NR2::AllH;
   my $iter=Photonic::LE;;NR2::AllH->new(geometry=>$geometry,nh=>$Nh,
            keepStates=>$save); 
   $iter->run;
   my $haydock_as=$iter->as;
   my $haydock_bs=$iter->bs;
   my $haydock_b2s=$iter->b2s;
   my $haydock_states=$iter->states;

=over 4

=item (for developers)

    package Photonic::LE::NR2::AllH.pm;
    $Photonic::LE::NR2::AllH::VERSION= '0.010';
    use namespace::autoclean;
    use Moose;
    has...
    with 'Photonic::Roles::OneH';

=back

=head1 DESCRIPTION

Roles consumed by AllH objects to be used in a Photonic
calculation. See also specific implementations. Iterates the
calculation of Haydock coefficients and states. 

=head1 METHODS

=over 4

=item * new(nh=>$nh, geometry=>$g, keepStates=>$k) 

Initializes an Ph::...::AllH object. $nh is the maximum number of desired
coefficients, $k is a flag, non zero to save the Haydock states. All
other arguments are as in Photonic::...::OneH.

=item * run

Runs the iteration to completion

=item * All the Photonic::NonRetarded::OneH methods

=back

=head1 ACCESORS (read only)

=over 4

=item * nh

Maximum number of desired Haydock 'a' coefficients and states. The
number of b coefficients is one less (but it starts with a dummy zero).

=item * keepStates

Flag to keep (1) or discard (0) Haydock states

=item * states

Array of Haydock states

=item * as

Array of Haydock a coefficients

=item * bs

Array of Haydock b coefficients.

=item * b2s

Array of Haydock b coefficients squared

=item * All the Photonic::...::OneH methods

=back

=begin Pod::Coverage

=head2 BUILD

=end Pod::Coverage

=cut

package Photonic::Roles::AllH;
$Photonic::Roles::AllH::VERSION = '0.010';
use Moose::Role;
use namespace::autoclean;
use Machine::Epsilon;
use PDL::Lite;
use PDL::NiceSlice;
#use Photonic::Utils qw(HProd);
#with 'Photonic::Roles::OneH';

has nh=>(is=>'ro', required=>1, 
         documentation=>'Maximum number of desired Haydock coefficients');
has 'keepStates'=>(is=>'ro', required=>1, default=>0, writer=> '_keepstates',
         documentation=>'flag to save Haydock states');
has states=>(is=>'ro', isa=>'ArrayRef[PDL::Complex]', 
         default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved states');
# why not pdl?
has as=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved a coefficients');
has bs=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved b coefficients');
has b2s=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved b^2 coefficients');
has reorthogonalize=>(is=>'ro', required=>1, default=>0,
         documentation=>'Reorthogonalize flag'); 
has 'previous_W' =>(is=>'ro', isa=>'PDL',
     writer=>'_previous_W', lazy=>1, init_arg=>undef,
     default=>sub{PDL->pdl([0])},
     documentation=>"Row of error matrix"
);
has 'current_W' =>(is=>'ro', isa=>'PDL',
     writer=>'_current_W', lazy=>1, init_arg=>undef,
     default=>sub{PDL->pdl([0])},
     documentation=>"Row of error matrix"
);
has 'next_W' =>(is=>'ro', isa=>'PDL',
     writer=>'_next_W', lazy=>1, init_arg=>undef,
     default=>sub {PDL->pdl([1])},
     documentation=>"Next row of error matrix"
);
has 'accuracy'=>(is=>'ro', default=>sub{machine_epsilon()},
                documentation=>'Desired or machine precision');
has 'noise'=>(is=>'ro', default=>sub{machine_epsilon()},
                documentation=>'Noise introduced each iteration');

sub BUILD {
    my $self=shift;
    # Can't reorthogonalize without previous states
    $self->_keepstates(1) if  $self->reorthogonalize;
}


#I use before and after trick (below), as a[n] is calculated together
#with b[n+1] in each iteration

#provided by OneH instance  
requires 'iterate', '_iterate_indeed', 'magnitude', 'innerProduct';

before 'states' => sub {
    my $self=shift;
    die "Can't return states unless keepStates!=0" unless $self->keepStates;
};

before '_iterate_indeed' => sub {
    my $self=shift;
    $self->_save_state;
    $self->_save_b2;
    $self->_save_b;
};
after '_iterate_indeed' => sub {
    my $self=shift;
    $self->_save_a;
    $self->_checkorthogonalize if $self->reorthogonalize;
};

sub run { #run the iteration
    my $self=shift;
    while($self->iteration < $self->nh && $self->iterate){
    }
}

sub _save_state {
    my $self=shift;
    return unless $self->keepStates;
    push @{$self->states}, $self->nextState;
}

sub _save_b {
    my $self=shift;
    push @{$self->bs}, $self->next_b;
}

sub _save_b2 {
    my $self=shift;
    push @{$self->b2s}, $self->next_b2;
}

sub _save_a {
    my $self=shift;
    push @{$self->as}, $self->current_a;
}

sub _checkorthogonalize {
    my $self=shift;
    return unless defined $self->nextState;
    my $a=PDL->pdl($self->as);
    my $b=PDL->pdl($self->bs);
    my $n=$self->iteration;
    $self->_previous_W(my $previous_W=$self->current_W);
    $self->_current_W(my $current_W=$self->next_W);
    my $next_W=PDL->pdl([]);
    if($n>=2){
	$next_W= $b->(1:-1)*$current_W->(1:-1) 
	    + ($a->(0:-2)-$a->(($n-1)))*$current_W->(0:-2)
	    - $b->(($n-1))*$previous_W;
	$next_W->(1:-1)+=$b->(1:-2)*$current_W->(0:-3) if ($n>=3);
	$next_W+=_sign($next_W)*2*$self->noise;
	$next_W/=$self->next_b;
    }
    $next_W=$next_W->append($self->noise) if $n>=1;
    $next_W=$next_W->append(1);
    $self->_next_W($next_W);
    return unless $n>=2;
    my $max=$next_W->(0:-2)->maximum;
    if($max > sqrt($self->accuracy)){
	$self->_orthogonalize;
	$current_W(0:-2).=$self->noise;
	$next_W(0:-2).=$self->noise;
    }
    $self->_current_W($current_W);
    $self->_next_W($next_W);
}

sub _orthogonalize {
    warn 'Hi there';
    my $self=shift;
    my $states=$self->states;
    my $currentState=$self->currentState;
    my $nextState=$self->nextState;
    pop(@{$states}); #get rid of currentState;
    for my $s (@{$states}){
	($nextState,$currentState)=
	    map {$_-$s*$self->innerProduct($s, $_)}($nextState, $currentState);
    }
    my $sc=$self->magnitude($currentState);
    # Is there a danger that orthogonalization finishes the iteration
    # early at current or next state?
    $currentState/=$sc;
    $self->_currentState($currentState);
    push @{$self->states}, $currentState;
    #necessary?: By the way, -= may not work yet
    $nextState=$nextState
	-$currentState*$self->innerProduct($currentState, $nextState); 
    my $sn=$self->magnitude($nextState);
    #What if nextstate disappears after reorthogonalization?
    $nextState/=$sn;
    $self->_nextState($nextState);
}

sub _sign {
    my $s=shift;
    return 2*($s>=0)-1;
}

1;
