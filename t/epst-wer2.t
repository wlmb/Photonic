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
use PDL::NiceSlice;
use Photonic::WE::R2::Metric;
use Photonic::WE::R2::Haydock;
use Photonic::WE::R2::Green;
use Photonic::WE::R2::GreenP;

use Test::More;
use lib 't/lib';
use TestUtils;

#Check epsilontensor for simple 1D system
#Non-retarded limit
my ($ea, $eb)=(pdl(1),3+4*i);
my $f=6/11;
my $B=zeroes(11,1)->xvals>=5;
my $g=Photonic::Geometry::FromB->new(B=>$B);
my $m=Photonic::WE::R2::Metric->new(
    geometry=>$g, epsilon=>$ea, wavenumber=>pdl(2e-5),
    wavevector=>pdl([1,0])*2.1e-5);
my $et=Photonic::WE::R2::Green->new(nh=>10, metric=>$m, epsB=>$eb);
my $etv=$et->epsilonTensor;
ok(Cagree($etv->((0),(0)), 1/((1-$f)/$ea+$f/$eb)),
			     "Long. perp. non retarded");
ok(Cagree($etv->((1),(1)), (1-$f)*$ea+$f*$eb),
			     "Trans. parallel non retarded");
ok(Cagree($etv->((0),(1)), 0), "xy k perp");
ok(Cagree($etv->((1),(0)), 0), "yx k perp");

$m=Photonic::WE::R2::Metric->new(
    geometry=>$g, epsilon=>$ea, wavenumber=>pdl(2e-5),
    wavevector=>pdl([0,1])*2.1e-5);
$et=Photonic::WE::R2::Green->new(nh=>10, metric=>$m, epsB=>$eb);
$etv=$et->epsilonTensor;
ok(Cagree($etv->((0),(0)), 1/((1-$f)/$ea+$f/$eb)),
			     "Trans. perp. non retarded");
ok(Cagree($etv->((1),(1)), (1-$f)*$ea+$f*$eb),
			     "Long. parallel non retarded");
ok(Cagree($etv->((0),(1)), 0), "xy k parallel");
ok(Cagree($etv->((1),(0)), 0), "yx k parallel");

#Compare to epsilon from transfer matrix.
#Construct normal incidence transfer matrix
($ea, $eb)=(r2C(1),r2C(2));
$g=Photonic::Geometry::FromB->new(B=>$B, L=>pdl(1,1));
my ($na, $nb)=map {sqrt($_)} ($ea, $eb);
my $q=1.2;
my ($ka,$kb)=map {$q*$_} ((1-$f)*$na, $f*$nb); #Multiply by length also
my $ma=pdl([cos($ka), -sin($ka)/$na],[$na*sin($ka), cos($ka)]);
my $mb=pdl([cos($kb), -sin($kb)/$nb],[$nb*sin($kb), cos($kb)]);
my $mt=($ma->(:,*1)*$mb->transpose->(:,:,*1))->sumover;
#Solve exact dispersion relation
my $cospd=($mt->((0),(0))+$mt->((1),(1)))/2;
my $sinpd=sqrt(1-$cospd**2);
my $pd=log($cospd+i()*$sinpd)/i;
warn "Bloch vector not real, $pd" unless $pd->im->abs < 1e-7;
$pd=$pd->re;
#epsilon from transfer matrix
my $epstm=($pd/$q)**2;
#epsilon from photonic
$m=Photonic::WE::R2::Metric->new(
    geometry=>$g, epsilon=>$ea->re, wavenumber=>pdl($q),
    wavevector=>pdl([$pd,0]));
$et=Photonic::WE::R2::Green->new(nh=>1000, metric=>$m,
  reorthogonalize=>1, epsB=>$eb);
$etv=$et->epsilonTensor->((1),(1));
ok(Cagree($etv, $epstm), "Epsilon agrees with transfer matrix");

my $h=Photonic::WE::R2::Haydock->new(nh=>1000, metric=>$m,
   polarization=>r2C(pdl [0,1]), reorthogonalize=>1);
$et=Photonic::WE::R2::GreenP->new(nh=>1000, haydock=>$h, epsB=>$eb);
$etv=$et->epsilon;
ok(Cagree($etv, $epstm, 1e-4), "Projected eps agrees with trans mat. Complex case.");

done_testing;
