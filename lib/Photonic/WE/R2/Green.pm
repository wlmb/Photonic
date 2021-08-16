package Photonic::WE::R2::Green;
$Photonic::WE::R2::Green::VERSION = '0.021';

=encoding UTF-8

=head1 NAME

Photonic::WE::R2::Green

=head1 VERSION

version 0.021

=head1 COPYRIGHT NOTICE

Photonic - A perl package for calculations on photonics and
metamaterials.

Copyright (C) 2016 by W. Luis Mochán

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 1, or (at your option)
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA  02110-1301 USA

    mochan@fis.unam.mx

    Instituto de Ciencias Físicas, UNAM
    Apartado Postal 48-3
    62251 Cuernavaca, Morelos
    México

=cut

=head1 SYNOPSIS

   use Photonic::WE::R2::Green;
   my $G=Photonic::WE::R2::Green->new(metric=>$m, nh=>$nh, epsB=>$epsB);
   my $GreenTensor=$G->greenTensor;
   my $WaveTensor=$G->waveOperator;

=head1 DESCRIPTION

Calculates the retarded green's tensor for a given fixed
L<Photonic::WE::R2::Metric> structure as a function of the dielectric
functions of the components. Includes the antysimmetric part, unless
it is not desired.

Can also calculate the macroscopic wave operator.

=head1 ATTRIBUTES

=over 4

=item * metric

L<Photonic::WE::R2::Metric> describing the structure and some parameters. (required)

=item * epsB

Dielectric function of component B (particle) (required)

=item * nh

The maximum number of Haydock coefficients to use.

=item * smallH, smallE

Criteria of convergence (default 1e-7) for Haydock coefficients and
continued fraction

=item * keepStates

Value of flag to keep Haydock states (default 0)

=item * symmetry

Flags that you only want the symmetric part of the Green's tensor, default 0.

=item * greenTensor

Macroscopic Green's operator for the dielectric functions of the particle
C<epsB>. The host's response C<epsA> is taken from the metric.

=item * epsA

Dielectric function of component A (host) taken from the metric.

=item * u

Spectral variable

=item * haydock

Array of L<Photonic::WE::R2::Haydock> structures, one for each
polarization, from the geometry.

=item * greenP

Array of L<Photonic::WE::R2::GreenP> structures, one for each direction

=item * nhActual

The actual number of Haydock coefficients used in the calculation

=item * converged

Flags that the calculation converged before using up all coefficients

=item * waveOperator

Returns the macroscopic wave operator for the Green tensor.

=back

=cut

use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use Photonic::WE::R2::Haydock;
use Photonic::WE::R2::GreenP;
use Photonic::Types;
use Photonic::Utils qw(tensor make_haydock make_greenp wave_operator);

use List::Util qw(all any);
use Moose;
use MooseX::StrictConstructor;

has 'metric'=>(is=>'ro', isa => 'Photonic::WE::R2::Metric',
       handles=>[qw(geometry B dims ndims r G GNorm L scale f)],required=>1);
has 'haydock' =>(is=>'ro', isa=>'ArrayRef[Photonic::WE::R2::Haydock]',
            init_arg=>undef, lazy=>1, builder=>'_build_haydock',
	    documentation=>'Array of Haydock calculators');
has 'greenP'=>(is=>'ro', isa=>'ArrayRef[Photonic::WE::R2::GreenP]',
             init_arg=>undef, lazy=>1, builder=>'_build_greenP',
             documentation=>'Array of projected G calculators');
has 'greenTensor'=>(is=>'ro', isa=>'PDL', init_arg=>undef,
             lazy=>1, builder=>'_build_greenTensor',
             documentation=>"Green's Tensor");
has 'converged'=>(is=>'ro', init_arg=>undef, writer=>'_converged',
             documentation=>
                  'All greenP evaluations converged in evaluation');
has 'reorthogonalize'=>(is=>'ro', required=>1, default=>0,
         documentation=>'Reorthogonalize haydock flag');
has 'nh' =>(is=>'ro', isa=>'Num', required=>1,
	    documentation=>'Desired no. of Haydock coefficients');
has 'smallH'=>(is=>'ro', isa=>'Num', required=>1, default=>1e-7,
	    documentation=>'Convergence criterium for Haydock coefficients');
has 'smallE'=>(is=>'ro', isa=>'Num', required=>1, default=>1e-7,
	    documentation=>'Convergence criterium for use of Haydock coeff.');
has 'epsA'=>(is=>'ro', isa=>'Photonic::Types::PDLComplex', init_arg=>undef, writer=>'_epsA',
    documentation=>'Dielectric function of host');
has 'epsB'=>(is=>'ro', isa=>'Photonic::Types::PDLComplex', required=>1,
        documentation=>'Dielectric function of inclusions');
has 'u'=>(is=>'ro', isa=>'Photonic::Types::PDLComplex', init_arg=>undef, writer=>'_u',
    documentation=>'Spectral variable');
has 'waveOperator' =>  (is=>'ro', isa=>'Photonic::Types::PDLComplex', init_arg=>undef,
             lazy=>1, builder=>'_build_waveOperator',
             documentation=>'Wave operator from evaluation');

with 'Photonic::Roles::KeepStates',  'Photonic::Roles::UseMask';

has 'cHaydock' =>(
    is=>'ro', isa=>'ArrayRef[Photonic::WE::R2::Haydock]',
    init_arg=>undef, lazy=>1, builder=>'_build_cHaydock',
    documentation=>'Array of Haydock calculators for complex projection');

has 'cGreenP'=>(
    is=>'ro', isa=>'ArrayRef[Photonic::WE::R2::GreenP]',
    init_arg=>undef, lazy=>1, builder=>'_build_cGreenP',
    documentation=>'Array of projected G calculators for complex projection');

has 'symmetric' => (
    is=>'ro', required=>1, default=>0,
    documentation=>'Flags only symmetric part required');

sub _build_greenTensor {
    my $self=shift;
    $self->_epsA(my $epsA=$self->metric->epsilon->r2C);
    my $epsB=$self->epsB;
    $self->_u(my $u=1/(1-$epsB/$epsA));
    $self->_converged(all { $_->converged } @{$self->greenP});
    my $greenTensor = tensor(pdl([map $_->Gpp, @{$self->greenP}]), $self->geometry->unitDyadsLU, $self->geometry->ndims, 2);
    #That's all unless you want the antisymmetric part
    return $greenTensor if $self->symmetric;
    my @greenPc = map $_->Gpp, @{$self->cGreenP}; ; #array of Green's projections along complex directions.
    $self->_converged(any { $_->converged } $self, @{$self->cGreenP});
    my $nd=$self->geometry->B->ndims;
    my $asy=$greenTensor->zeroes; #xy,xy, $ndx$nd
    my $cpairs=$self->geometry->cUnitPairs;
    my $m=0;
    for my $i(0..$nd-2){
	for my $j($i+1..$nd-1){
	    my $pair=$cpairs->(:,($m));
	    #$asy is xy,xy. First index is column
	    $asy(($i), ($j)).=PDL->i()*(
		$greenPc[$m]-
		($pair->conj->(*1) # column, row
		 *$pair->(:,*1)
		 *$greenTensor)->sumover->sumover
		);
	    $asy(($j), ($i)).=-$asy(($i),($j));
	    $m++
	}
     }
    #print $asy, "\n";
    $greenTensor+$asy;
}

sub _build_haydock { # One Haydock coefficients calculator per direction0
    my ($self) = @_;
    make_haydock($self, 'Photonic::WE::R2::Haydock', $self->geometry->unitPairs, 0, qw(reorthogonalize use_mask mask));
}

sub _build_greenP {
    make_greenp(shift, 'Photonic::WE::R2::GreenP', undef, qw(epsB));
}

sub _build_cHaydock {
    # One Haydock coefficients calculator per complex polarization
    my $self=shift;
    make_haydock($self, 'Photonic::WE::R2::Haydock', $self->geometry->cUnitPairs, 0);
}

sub _build_cGreenP {
    make_greenp(shift, 'Photonic::WE::R2::GreenP', 'cHaydock', qw(epsB));
}

sub _build_waveOperator {
    my $self=shift;
    wave_operator($self->greenTensor, $self->geometry->ndims);
}

__PACKAGE__->meta->make_immutable;

1;
