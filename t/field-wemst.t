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

use strict;
use warnings;
use PDL;
use PDL::Constants qw(PI);
use Photonic::WEM::ST::Haydock;
use Photonic::WEM::ST::Metric;
use Photonic::WEM::ST::Field;
use Photonic::Geometry::FromEpsilonTensor;
use Photonic::Geometry::FromB;

use Test::More tests => 4;
use lib 't/lib';
use TestUtils;

my $ea=1+2*i;
my $eb=3+4*i;
my $ident5 = identity(3)->slice(":,:,*1,*1,*1"); #xyz:xyz:1:1:1
my $mua = $ident5->r2C; #xyz:xyz:1:1:1
my $mub = $ident5->r2C; #xyz:xyz:1:1:1

my $da = 0.4;
my $db = (1-$da);

# Check for simple 1D system (period=filma+filmb)
# Field
{
  #Longitudinal case (effective medium approximation)
  my $grid = zeroes(11,1,1);
  my $x    = $grid->xvals/11;

  ## add two dimensions at the beggining
  my $filma = ($x<=$da)->dummy(0)->dummy(0);
  my $filmb = ($x>$da)->dummy(0)->dummy(0);

  my $epsilonl = $filmb*($eb*$ident5) + $filma*($ea*$ident5);
  my $mul = $filmb*$mub + $filma*$mua;

  my $glf=Photonic::Geometry::FromB->new(B=>$filmb); #only used to get the correct filling fraction
  my $gl=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilonl, L=>pdl([1,1,1])); #long
  my $ml=Photonic::WEM::ST::Metric->new(mu=>$mul, geometry=>$gl, epsilon=>pdl(1),
  				      wavenumber=>pdl(0.01), wavevector=>pdl([0.01,0,0]));
  my $haydock=Photonic::WEM::ST::Haydock->new(
      metric=>$ml, nh=>10, keepStates=>1, polarization=>pdl([1,0,0])->r2C,
      epsilon=>$epsilonl);
  my $flo=Photonic::WEM::ST::Field->new(haydock=>$haydock, nh=>10);
  my $flv=($flo->field)->squeeze;
  my $fla = (1/$ea);
  my $flb = (1/$eb);
  my $fproml=$fla*(1-$glf->f)+$flb*($glf->f);
  ($fla, $flb)=map {$_/$fproml} ($fla, $flb);
  my $flx=pdl([($fla*$filma+$flb*$filmb),0,0])->mv(-1,0)->squeeze;
  ok(Cagree($flv, $flx), "1D long field");
}

{
  #View 2D from 1D superlattice. Long wavelength transverse case
  my $grid = zeroes(1,1,11);
  my $z    = $grid->zvals/11;

  ## add two dimensions at the beggining
  my $filma = ($z<=$da)->dummy(0)->dummy(0);
  my $filmb = ($z>$da)->dummy(0)->dummy(0);

  my $epsilont = $filmb*($eb*$ident5) + $filma*($ea*$ident5);
  my $mut = $filmb*$mub + $filma*$mua;

  my $gt=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilont, L=>pdl([1,1,1])); #long
  my $mt=Photonic::WEM::ST::Metric->new(mu=>$mut, geometry=>$gt, epsilon=>pdl(1),
  				      wavenumber=>pdl(0.0001), wavevector=>pdl([0,0,0.0001]));
  my $nt=Photonic::WEM::ST::Haydock->new(
      metric=>$mt, nh=>10, keepStates=>1, polarization=>pdl([1,0,0])->r2C,
      epsilon=>$epsilont);
  my $fto=Photonic::WEM::ST::Field->new(haydock=>$nt, nh=>10);
  my $ftv=($fto->field);
  my $ftx=r2C(pdl [1,0,0])->dummy(3,11);
  ok(Cagree($ftv, $ftx), "1D transversal field");
}

# Rawfields
{
  #Longitudinal case (effective medium approximation)
  my $grid = zeroes(11,1,1);
  my $x    = $grid->xvals/11;

  ## add two dimensions at the beggining
  my $filma = ($x<=$da)->dummy(0)->dummy(0);
  my $filmb = ($x>$da)->dummy(0)->dummy(0);

  my $epsilonl = $filmb*($eb*$ident5) + $filma*($ea*$ident5);
  my $mul = $filmb*$mub + $filma*$mua;

  my $gl=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilonl, L=>pdl([1,1,1])); #long
  my $ml=Photonic::WEM::ST::Metric->new(mu=>$mul, geometry=>$gl, epsilon=>pdl(1),
  				      wavenumber=>pdl(1), wavevector=>pdl([0.01,0,0]));
  my $haydock=Photonic::WEM::ST::Haydock->new(
      metric=>$ml, nh=>10, keepStates=>1, polarization=>pdl([1,0,0])->r2C,
      epsilon=>$epsilonl);
  my $flo=Photonic::WEM::ST::Field->new(haydock=>$haydock, nh=>10);
  my $flv = ($flo->rawfield)->squeeze;
  my $fla = (1/$ea);
  my $flb = (1/$eb);
  my $flx=pdl(-4*PI*($fla*$filma+$flb*$filmb),0,0)->mv(-1,0)->squeeze;
  ok(Cagree($flv, $flx), "1D long rawfield");
}

{
  #View 2D from 1D superlattice. Long wavelength transverse case
  my $grid = zeroes(1,1,11);
  my $z    = $grid->zvals/11;

  ## add two dimensions at the beggining
  my $filma = ($z<=$da)->dummy(0)->dummy(0);
  my $filmb = ($z>$da)->dummy(0)->dummy(0);

  my $epsilont = $filmb*($eb*$ident5) + $filma*($ea*$ident5);
  my $mut = $filmb*$mub + $filma*$mua;

  my $q=pdl(0.00001);
  my $k=pdl([0,0,0.0002]);

  my $glf=Photonic::Geometry::FromB->new(B=>$filmb); #only used to get the correct filling fraction
  my $gt=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilont, L=>pdl([1,1,1])); #long
  my $mt=Photonic::WEM::ST::Metric->new(mu=>$mut, geometry=>$gt, epsilon=>pdl(1),
  				      wavenumber=>$q, wavevector=>$k);
  my $nt=Photonic::WEM::ST::Haydock->new(
      metric=>$mt, nh=>10, keepStates=>1, polarization=>pdl([1,0,0])->r2C,
      epsilon=>$epsilont);
  my $fto=Photonic::WEM::ST::Field->new(haydock=>$nt, nh=>10);
  my $ftv=($fto->rawfield)->squeeze;
  my $epsM = (1-$glf->f)*$ea+($glf->f)*$eb;
  my $ftx=-4*PI*r2C(pdl [1, 0,0])->dummy(3,11)*$q**2/($epsM*$q**2-$k->inner($k));
  ok(Cagree($ftv, $ftx), "1D trans rawfield");
}