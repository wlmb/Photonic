=head1 NAME

Photonic::WE::S::WaveP

=head1 VERSION

version 0.011

=head1 SYNOPSIS

   use Photonic::WE::S::WaveP;
   my $W=Photonic::WE::S::WaveP->new(haydock=>$h, nh=>$nh, smallE=>$s);
   my $WaveProjection=$W->evaluate($epsB);

=head1 DESCRIPTION

Calculates the macroscopic projected wave operator for a given fixed
Photonic::WE::S::AllH structure as a function of the dielectric
functions of the components. 

NOTE: Only works along principal directions, as it treats Green's
function as scalar.  

=head1 METHODS

=over 4

=item * new(haydock=>$h, nh=>$nh, smallE=>$smallE)  

Initializes the structure.

$h is Photonic::WE::S::AllH structure (required).

$nh is the maximum number of Haydock coefficients to use (required).

$smallE is the criteria of convergence (default 1e-7).

=item * evaluate($epsB)

Returns the macroscopic wave operator for a given value of the
dielectric functions of the particle $epsB. The host's
response $epsA is taken from the metric.  

=back

=head1 ACCESORS (read only)

=over 4

=item * waveOperator

The macroscopic wave operator of the last operation

=item * All accesors of Photonic::WE::S::GreenP


=back

=cut

package Photonic::WE::S::WaveP;
$Photonic::WE::S::WaveP::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Types;
use Moose;
use MooseX::StrictConstructor;

extends 'Photonic::WE::S::GreenP'; 

has 'waveOperator' =>  (is=>'ro', isa=>'PDL::Complex', init_arg=>undef,
             lazy=>1, builder=>'_build_waveOperator',   
             documentation=>'Wave operator');

sub _build_waveOperator {
    my $self=shift;
    my $greenP=$self->Gpp;
    my $wave=1/$greenP; #only works along principal directions!!
    return $wave;
};

__PACKAGE__->meta->make_immutable;
    
1;

__END__
