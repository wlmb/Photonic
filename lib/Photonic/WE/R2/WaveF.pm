=head1 NAME

Photonic::WE::R2::WaveF

=head1 VERSION

version 0.011

=head1 SYNOPSIS

   use Photonic::WE::R2::WaveF;
   my $W=Photonic::WE::R2::WaveF->new(metric=>$m, nh=>$nh);
   my $WaveTensor=$W->evaluate($epsB);

=head1 DESCRIPTION

Calculates the asymmetric macroscopic wave operator for a given fixed
Photonic::WE::R2::Metric structure as a function of the dielectric
functions of the components.

=head1 METHODS

=over 4

=item * new(metric=>$m, nh=>$nh, smallH=>$smallH, smallE=>$smallE, keepStates=>$k)

Initializes the structure.

$m Photonic::WE::R2::Metric describing the structure and some parametres.

$nh is the maximum number of Haydock coefficients to use.

$smallH and $smallE are the criteria of convergence (default 1e-7) for
Haydock coefficients and for continued fraction.

$k is a flag to keep states in Haydock calculations (default 0)

=item * evaluate($epsB)

Returns the macroscopic wave operator for a given value of the
dielectric functions of the particle $epsB. The host's
response $epsA is taken from the metric.

=back

=head1 ACCESORS (read only)

=over 4

=item * waveOperator

The macroscopic wave operator of the last operation

=item * All accesors of Photonic::WE::R2::Green


=back

=cut

package Photonic::WE::R2::WaveF;
$Photonic::WE::R2::WaveF::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use PDL::MatrixOps;
use Storable qw(dclone);
use PDL::IO::Storable;
#use Photonic::WE::R2::AllH;
use Photonic::Types;
use Moose;
use MooseX::StrictConstructor;

extends 'Photonic::WE::R2::GreenF';

has 'waveOperator' =>  (is=>'ro', isa=>'PDL::Complex', init_arg=>undef,
             writer=>'_waveOperator',
             documentation=>'Wave operator from last evaluation');

around 'evaluate' => sub {
    my $orig=shift;
    my $self=shift;
    my $green=$self->$orig(@_);
    #make a real matrix from [[R -I][I R]] to solve complex eq.
    my $greenreim=$green->re->append(-$green->im)
       ->glue(1,$green->im->append($green->re))->sever; #copy vs sever?
    my($lu, $perm, $par)=$greenreim->lu_decomp;
    my $d=$self->geometry->ndims;
    my $idreim=identity($d)->glue(1,PDL->zeroes($d,$d))->mv(0,-1);
    my $wavereim=lu_backsub($lu,$perm,$par,$idreim);
    my $wave=$wavereim->reshape($d,2,$d)->mv(1,0)->complex;
    $self->_waveOperator($wave);
    return $wave;
};

__PACKAGE__->meta->make_immutable;

1;

__END__
