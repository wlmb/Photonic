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

=item * new(geometry=>$g, epsilon=>$e, $wavenumber=>$q, $wavevector=>$k);

Create a new Ph::WEM::S::Metric object with Geometry $g, dielectric
function of the host $e, vacuum wavenumber $q=omega/c  and wavevector
$k. $q and $k are real.

= item * applyMetric(psi)

Create and apply the magnetic metric to teh state psi provided in
Haydock when it calls this method.

=back


=cut


use namespace::autoclean;
use PDL::Lite;
use PDL::MatrixOps;
use PDL::NiceSlice;
use Carp;
use Photonic::Types -all;
use Photonic::Utils qw(any_complex);
use Moo;
use MooX::StrictConstructor;

has 'value' =>(is=>'lazy', isa=>PDLObj, init_arg=>undef,
	       documentation=>'Build a value so Roles does not complain');
has 'mu'    =>(is=>'ro', isa=>PDLObj, required=>1,
	       documentation=> 'Magnetic permeability scalar function');


with 'Photonic::Roles::Metric';
# Roles::Metric Pulls the following attributes: geometry, epsilonRef,
# wavenumber, and wavevector. All metrics need these parameters.


sub apply{
    # Evaluate the metric tensor applied to the state
    my $self = shift;
    my $psi = shift; #xyz:pm:nx:ny:nz
    my $G=$self->G; #reciprocal lattice, xyz:nx:ny:nz
    my $q=$self->wavenumber; # wavenumber
    my $eps=$self->epsilon; # this is the reference epsilon (from Roles)
    croak "For the time being the reference epsilon should be 1 or not initialized"
	unless $eps==1;
    my $k=$self->wavevector; # bloch wavevector, xyz
    croak "Wave vector must be ".$self->ndims."-dimensional vector" unless
	$k->dim(0)==$self->ndims;
    my $mu=$self->mu; # magnetic permeability nx:ny:nz
    croak "Mu must be a $self->dims array" unless
	(pdl($mu->dims)==pdl($self->dims))->all; #exactly the same dimensions
    # FIRST TERM of the metric
    # apply longitudinal projector to state psi
    my $ProjL_psi = ($self->longitudinal_projector #xyz:xyz:pm:nx:ny:nz
		     ->inner($psi(:,*1) #xyz:1:pm:nx:ny:nz
	); #xyz:pm:nx:ny:nz a this matrix vector mutiplication
    # SECOND TERM
    # wavevectors with +-k
    my ($kPG,$kMG) = ($G+$k, $G-$k); # xyz:nx:ny:nz
    # magnitude of Gpmk vectors squared
    my ($kPG2,$kMG2) = map $_->inner($_), $kPG, $kMG; #nx:ny:nz;
    # divide each vector by its norm squared and concatenate them
    my $kPMG_norm2 = cat($kPG/$kPG2->(*1),$kMG/$kMG2->(*1))->mv(-1,1); #xyz:pm:nx:ny:nz
    # cross prod with psi
    my $kPMGXpsi = crossp($kPMG_norm2,$psi); #xyz:pm:nx:ny:nz
    # FT to real space and move dims 
    my $kPMGXpsi_r = GtoR($kPMGXpsi,$self->ndims,2)->mv(0,-1)->mv(0,-1); #nx:ny:nz:xyz:pm
    # left-multiply by mu and move dims back
    my $mu_kPMGXpsi_r = (($mu)*($kPMGXpsi_r))->mv(-1,0)->mv(-1,0); #xyz:pm:nx:ny:nz
    # mu kPMGtimespsi in G space 
    my $mu_kPMGXgpsi_G = RtoG($mu_kPMGXpsi_r,$self->ndims,2); #xyz:pm:nx:ny:nz
    # cross produt with wavvectors 
    my $kPMGXmukPMGXpsi = crossp($kPMG_norm2, $mu_kPMGXpsi_G); #xyz::pm:nx:ny:nz
    # multiply by q**2
    my $mgn_part_G_psi = ($q*$q)*$kPMGXmukPMGXpsi; #xyz::pm:nx:ny:nz
    # complete applied metric in G space
    my $gpsi_G = $ProjL_psi + $mgn_part_G_psi; #xyz:nx:ny:nz
    return $gpsi_G;
}

#this is built here because it uses the k vector and it's an atributte
#because it is called in Haydock to build the Hamiltonian
sub longitudinal_projector{
    #returns the value of the projector
    my $self=shift;
    my $G=$self->G; #reciprocal lattice from Geometry
    my $k=$self->wavevector; # bloch wavevector, solution
    my ($kPG,$kMG) = ($k+$G, $k-$G); # xyz:pm:nx:ny:nz
    # spinorlike wavevectors,
    my ($kPGkPG, $kMGkMG) = map $_->outer($_), $kPG, $kMG; #xyz:xyz:nx:ny:nz
    #outer product GG to build
    #the longitudinal projector 
    my ($kPG2,$kMG2) = map $_->inner($_), $kPG, $kMG; #nx:ny:nz;
    #magnitude squared of G+-k vectors
    my $ProjLPMk = cat($kPGkPG/$kPG2->(*1),$kMGkMG/$kMG2->(*1)) # xyz:xyz:nx:ny:nz:pm
	->mv(-1,2); # xyz:xyz:pm:nx:ny:nz
    #normalized longitudinal projectors for each G+-k wave vectors
    return $ProjLPMk;
}


1;
