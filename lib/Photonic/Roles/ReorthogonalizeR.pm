package Photonic::Roles::ReorthogonalizeR;
$Photonic::Roles::ReorthogonalizeR::VERSION = '0.010';
use Moose::Role;
use Machine::Epsilon;
use PDL::Lite;
use PDL::NiceSlice;
use List::MoreUtils qw(pairwise);

has 'previous_W' =>(is=>'ro', 
     writer=>'_previous_W', lazy=>1, init_arg=>undef,
     default=>sub{PDL->pdl([0])},
     documentation=>"Row of error matrix"
);
has 'current_W' =>(is=>'ro',
     writer=>'_current_W', lazy=>1, init_arg=>undef,
     default=>sub{PDL->pdl([0])},
     documentation=>"Row of error matrix"
);
has 'next_W' =>(is=>'ro',
     writer=>'_next_W', lazy=>1, init_arg=>undef,
     default=>sub {PDL->pdl([1])},
     documentation=>"Next row of error matrix"
);
has 'accuracy'=>(is=>'ro', default=>sub{machine_epsilon()},
                documentation=>'Desired or machine precision');
has 'noise'=>(is=>'ro', default=>sub{machine_epsilon()},
    documentation=>'Noise introduced each iteration to overlap matrix');
has 'normOp'=>(is=>'ro', required=>1, default=>1,
	       documentation=>'Estimate of operator norm'); 
has fullorthogonalize_N=>(is=>'ro', init_arg=>undef, default=>0,
			  writer=>'_fullorthogonalize_N',
			  documentation=>'# desired reorthogonalizations'); 
has 'orthogonalizations'=>(is=>'ro', init_arg=>undef, default=>0,
			   writer=>'_orthogonalizations');

around '_fullorthogonalize_indeed' => sub {
    my $orig=shift; #won't use
    my $self=shift;
    my $psi=shift; #state to orthogonalize
    return $psi unless $self->fullorthogonalize_N;
    $self->_fullorthogonalize_N($self->fullorthogonalize_N-1);
    $self->_orthogonalizations($self->orthogonalizations+1);
    foreach(pairwise {[$a, $b]} @{$self->states}, @{$self->gs}){
	#for every saved state
	my ($s, $g)=($_->[0], $_->[1]); #state, metric
	$psi=$psi-$g*$self->innerProduct($s, $psi)*$s;
    }
    return $psi;
};

sub _checkorthogonalize {
    my $self=shift;
    return unless defined $self->nextState;
    return unless $self->reorthogonalize;
    return if $self->fullorthogonalize_N; #already orthogonalizing
    my $n=$self->iteration;
    my $a=PDL->pdl($self->as);
    my $b=PDL->pdl($self->bs);
    my $c=PDL->pdl($self->cs);
    my $g=PDL->pdl($self->gs);
    $self->_previous_W(my $previous_W=$self->current_W);
    $self->_current_W(my $current_W=$self->next_W);
    my $next_W=PDL->pdl([]);
    if($n>=2){
	$next_W= $b->(1:-1)*$current_W->(1:-1) 
	    + ($a->(0:-2)-$a->(($n-1)))*$current_W->(0:-2)
	    - $c->(($n-1))*$previous_W;
	$next_W->(1:-1)+=$c->(1:-2)*$current_W->(0:-3) if ($n>=3);
	$next_W=$next_W+_sign($next_W)*2*$self->normOp*$self->noise;
	$next_W=$next_W/$self->next_b;
    }
    $next_W=$next_W->append($self->noise) if $n>=1;
    $next_W=$next_W->append(1);
    $self->_next_W($next_W);
    return unless $n>=2;
    my $max=$next_W->(0:-2)->maximum;
    if($max > sqrt($self->accuracy)){
	#recalculate the las two states with full reorthogonalization
	$self->_fullorthogonalize_N(3); #2 states, but check until 3d state
	$self->_pop; #undoes stack 2 steps
	$self->_pop;
	$current_W(0:-2).=$self->noise;
	$next_W(0:-2).=$self->noise;
    }
    $self->_current_W($current_W);
    $self->_next_W($next_W);
}

sub _sign {
    my $s=shift;
    return 2*($s>=0)-1;
}

1;
