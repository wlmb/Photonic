package Photonic::WEM::S::Metric;
$Photonic::WEM::S::Metric::VERSION = '0.024';

=encoding UTF-8

=head1 NAME

Photonic::WEM::S::Metric

=head1 VERSION

version 0.024

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

    use Photonic::WEM::S::Metric;
    my $gGG=Photonic::WEM::S::Metric->new(
            geometry=>$geometry, epsilon=>$eps,
            wavenumber => $q, $wavevector=>k);
    f($gGG->value);

=head1 DESCRIPTION

Calculates the retarded metric tensor g_{GG'}^{ij} for use in the
calculation of the retarded Haydock coefficients for the wave equation
in a binary medium where the host has no dissipation.

=head1 METHODS

=over 4

=item * new(mu=>$mu, geometry=>$g, epsilon=>$e, $wavenumber=>$q, $wavevector=>$k);

Create a new Ph::WEM::S::Metric object with permeability $mu, Geometry $g, dielectric
function of the host $e, vacuum wavenumber $q=omega/c  and wavevector
$k. $q and $k are real.

= item * apply(psi)

Create and apply the magnetic metric to the state psi provided in
Haydock when it calls this method.

=back


=cut


use namespace::autoclean;
use PDL::Lite;
use PDL::Primitive;
use PDL::MatrixOps;
use PDL::NiceSlice;
use Carp;
use Photonic::Types -all;
use Photonic::Utils qw(any_complex mvN GtoR RtoG);
use Moo;
use MooX::StrictConstructor;

has 'mu' =>(is=>'ro', isa=>PDLObj, required=>1,
	    documentation=> 'Magnetic permeability, scalar function of position');

has 'LP'
    =>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
       documentation=>'Longitudinal Projector in reciprocal spinorial space');
has '_G_pm_k_div_sqrd'
    =>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
       documentation=>'Internal. Reciprocal vector in spinor space divided by its square');

has  'value'=>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
	       documentation=>'Metric Tensor');

with 'Photonic::Roles::Metric';
# Roles::Metric Pulls the following attributes: geometry, epsilonRef,
# wavenumber, and wavevector. All metrics need these parameters.

sub BUILD {
    my $self=shift;
    croak "Currenlty, only 3D is supported" unless $self->ndims==3;
    my $eps=$self->epsilon;
    croak "For the time being the reference epsilon should be 1 or not initialized"
	unless $eps==1;
    my $k=$self->wavevector; # bloch wavevector, xyz
    die "Wavevector should be 3D" unless $k->dims==1 && $k->dim(0)==3;
    my $mu=$self->mu;
    die "Mu should be a $self->dims scalar function" unless
	(pdl($mu->dims)==pdl($self->dims))->all; #exactly the same dimensions
}

sub apply{
    # Evaluate the metric tensor applied to the state
    my $self = shift;
    my $psi = shift; #xyz:pm:nx:ny:nz
    my $q=$self->wavenumber; # wavenumber
    my $mu=$self->mu; # magnetic permeability nx:ny:nz
    my $ndims=$self->ndims;
    # FIRST TERM of the metric
    # apply longitudinal projector to state psi
    my $long_part = $self->LP #xyz:xyz:pm:nx:ny:nz
	->inner($psi(:,*1)); #xyz:1:pm:nx:ny:nz;
	#xyz:pm:nx:ny:nz a this matrix vector mutiplication
    # SECOND TERM
    my $G_pm_k_div_sqr = $self->_G_pm_k_div_sqrd; #(G pm k)/|G pm k|^2
    # cross prod with psi
    my $G_X_psi_G = crossp($G_pm_k_div_sqr, $psi); #xyz:pm:nx:ny:nz
    # FT to real space and move dims
    my $G_X_psi_r = GtoR(mvN($G_X_psi_G, 0,1,-1), $ndims,0); #nx:ny:nz:xyz:pm
    # left-multiply by mu and move dims back
    my $mu_G_X_psi_r = mvN($mu*($G_X_psi_r),-2,-1,0); #xyz:pm:nx:ny:nz
    # mu G X psi in G space
    my $mu_G_X_psi_G = RtoG($mu_G_X_psi_r, $ndims,2); #xyz:pm:nx:ny:nz
    # cross produt with wavvectors
    my $G_X_mu_G_X_psi = crossp($G_pm_k_div_sqr, $mu_G_X_psi_G); #xyz::pm:nx:ny:nz
    # multiply by q**2
    my $trans_part = ($q*$q)*$G_X_mu_G_X_psi; #xyz::pm:nx:ny:nz
    # complete applied metric in G space
    return $long_part + $trans_part; #xyz:nx:ny:nz
}

#this is built here because it uses the k vector and it's an atributte
#because it is called in Haydock to build the Hamiltonian
sub _build_LP{
    #returns the value of the projector
    my $self=shift;
    my $G=$self->G; #reciprocal lattice from Geometry
    my $k=$self->wavevector; # bloch wavevector, solution
    return pdl(map{$_->outer($_)/$_->inner($_)->(*1,*1)}($G+$k, $G-$k))
	->mv(-1,2); # (G+-k)(G+-k)/(G+-k)^2 xy:xy:pm:nx:...
}

sub _build__G_pm_k_div_sqrd { # (G+-k)/(G+-k)^2
    my $self=shift;
    my $k=$self->wavevector;
    my $G=$self->G;
    return pdl(map{$_/$_->inner($_)->(*1)}($G+$k, $G-$k))->mv(-1,1);
}

sub _build_value {
    croak "Sorry, there is no value for the metric of magnetizable systems";
}


1;
