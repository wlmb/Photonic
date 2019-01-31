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

=item * new(geometry=>$g, nh=>$nh) 

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
number of b coefficients is one less.

=item * as

Array of Haydock a coefficients

=item * bs

Array of Haydock b coefficients

=item * b2s

Array of Haydock b coefficients squared

=item * All the Photonic::...::OneH methods

=back

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
# why not pdl?
has as=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved a coefficients');
has bs=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved b coefficients');
has b2s=>(is=>'ro', isa=>'ArrayRef[Num]', default=>sub{[]}, init_arg=>undef,
         documentation=>'Saved b^2 coefficients');

#I use before and after trick (below), as a[n] is calculated together
#with b[n+1] in each iteration

requires 'iterate', '_iterate_indeed';

before '_iterate_indeed' => sub {
    my $self=shift;
    $self->_save_b2;
    $self->_save_b;
};
after '_iterate_indeed' => sub {
    my $self=shift;
    $self->_save_a;
};

sub run { #run the iteration
    my $self=shift;
    while($self->iteration < $self->nh && $self->iterate){
    }
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

1;
