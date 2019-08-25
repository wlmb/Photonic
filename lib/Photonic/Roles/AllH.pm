=head1 NAME

Photonic::Roles::AllH

=head1 VERSION

version 0.011

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

    package Photonic::LE::NR2::AllH;
    $Photonic::LE::NR2::AllH::VERSION= '0.011';
    use namespace::autoclean;
    use Moose;
    has...
    with 'Photonic::Roles::OneH';

=back

=head1 DESCRIPTION

Roles consumed by AllH objects to be used in a Photonic
calculation. See also specific implementations. Uses Iterator to iterate the
calculation of Haydock coefficients and states, and saves them for later retrieval.

=head1 METHODS

=over 4

=item * new(nh=>$nh, geometry=>$g, keepStates=>$k) 

Initializes an Ph::...::AllH object. $nh is the maximum number of desired
coefficients, $k is keepstates non zero flag to save the Haydock states. All other arguments are as in Photonic::...::OneH.

=item * run

Runs the iteration to completion, tells $iter to start the calculation and get values of variables you search.

=item * All the Photonic::...::OneH methods

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
$Photonic::Roles::AllH::VERSION = '0.011';
use Machine::Epsilon;
use Photonic::Iterator qw(:all);
use PDL::Lite;
use PDL::Complex;
use PDL::NiceSlice;
use IO::File;
use Storable qw(store_fd);
use PDL::IO::Storable;
use Carp;
use Moose::Role;

has nh=>(is=>'ro', required=>1, 
         documentation=>'Maximum number of desired Haydock coefficients');
has 'keepStates'=>(is=>'ro', required=>1, default=>0, writer=> '_keepstates',
         documentation=>'flag to keep Haydock states');
has _states=>(is=>'ro', isa=>'ArrayRef[PDL::Complex]', 
         default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved states');

# why not pdl?
has as=>(is=>'ro', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved a coefficients');
has bs=>(is=>'ro', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved b coefficients');
has b2s=>(is=>'ro', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved b^2 coefficients');
has cs=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved c coefficients');
has bcs=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved b*c coefficients');
has gs=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved g coefficients');
has reorthogonalize=>(is=>'ro', required=>1, default=>0,
         documentation=>'Reorthogonalize flag'); 
has 'stateFN'=>(is=>'ro', required=>1, default=>undef, 
		documentation=>'Filename to save Haydock states');
has 'stateFD'=>(is=>'ro', init_arg=>undef, builder=>'_build_stateFD',
		lazy=>1, documentation=>'Filedescriptor to save
		Haydock states');  

#provided by OneH instance  
requires qw(iterate _iterate_indeed magnitude innerProduct
    _checkorthogonalize);


sub BUILD {
    my $self=shift;
    # Can't reorthogonalize without previous states
    $self->_keepstates(1) if  $self->reorthogonalize;
}

sub state_iterator {
    my $self=shift;
    confess "Can't return states unless keepStates!=0"
	unless $self->keepStates;
    my $n=0;
    my $s=$self->_states;
    my $f=iterator { #closure
	return $s->[$n++];
    };
    return $f;
}

sub states {
    croak "Don't use 'states'. Use state_iterator instead";
}

sub _build_stateFD {
    my $self=shift;
    my $fn=$self->stateFN;
    croak "You didn't provide a stateFN" unless defined $fn;
    my $fh=IO::File->new($fn, "w+")
	or croak "Couldn't open ".$self->stateFN.": $!";
    return $fh;
}

    
#I use before and after trick (below), as a[n] is calculated together
#with b[n+1] in each iteration

#before 'states' => sub {
#    my $self=shift;
#    confess "Can't return states unless keepStates!=0" unless $self->keepStates;
#};

before '_iterate_indeed' => sub {
    my $self=shift;
    $self->_save_state;
    $self->_save_b2;
    $self->_save_b;
    $self->_save_c;
    $self->_save_bc;
    $self->_save_g;
};
after '_iterate_indeed' => sub {
    my $self=shift;
    $self->_save_a;
    $self->_checkorthogonalize;
};

sub run { #run the iteration
    my $self=shift;
    while($self->iteration < $self->nh && $self->iterate){
    }
}

sub _save_state {
    my $self=shift;
    push @{$self->_states}, $self->nextState if $self->keepStates;
    return unless defined $self->stateFN;
    my $fh=$self->stateFD;
    store_fd \$self->nextState, $fh or croak "Couldn't store state: $!";
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

sub _save_bc {
    my $self=shift;
    push @{$self->bcs}, $self->next_bc;
}

sub _save_c {
    my $self=shift;
    push @{$self->cs}, $self->next_c;
}

sub _save_g {
    my $self=shift;
    push @{$self->gs}, $self->next_g;
}

sub _pop { # undo the changes done after, in and before iteration, for
	   # reorthogonalization, in reverse order
    my $self=shift;
    $self->_nextState(pop @{$self->_states});
    $self->_currentState($self->_states->[-1]);
    $self->_previousState($self->_states->[-2]//r2C(0));
    $self->_current_a(pop @{$self->as});
    $self->_next_b2(pop @{$self->b2s});
    $self->_current_b2($self->b2s->[-1]);
    $self->_next_b(pop @{$self->bs});
    $self->_current_b($self->bs->[-1]);
    $self->_next_c(pop @{$self->cs});
    $self->_current_c($self->cs->[-1]);
    $self->_next_bc(pop @{$self->bcs});
    $self->_next_g(pop @{$self->gs});
    $self->_current_g($self->gs->[-1]);
    $self->_previous_g($self->gs->[-2]);
    $self->_iteration($self->iteration-1); #decrement counter
}    

no Moose::Role;

1;
