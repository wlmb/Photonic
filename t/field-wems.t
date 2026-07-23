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
use Photonic::WEM::S::Haydock;
use Photonic::WEM::S::Metric;
use Photonic::WEM::S::Field;

use Test::More tests => 4;
use lib 't/lib';
use TestUtils;

my $ea=1+2*i;
my $eb=3+4*i;

#Check field for simple 1D system.
{
    #Longitudinal case
    my $B=zeroes(11,1,1)->xvals<5; #1D system
    my $epsilon=$ea*(1-$B)+$eb*$B;
    my $mu=ones(11,1,1)->r2C;
    my $gl=Photonic::Geometry::FromB->new(B=>$B); #long
    my $ml=Photonic::WEM::S::Metric->new(mu=>$mu, geometry=>$gl, epsilon=>pdl(1),
					wavenumber=>pdl(1), wavevector=>pdl([0.01,0,0]));
    my $haydock=Photonic::WEM::S::Haydock->new(
	metric=>$ml, nh=>10, keepStates=>1, polarization=>pdl([1,0,0])->r2C,
	epsilon=>$epsilon);
    my $flo=Photonic::WEM::S::Field->new(haydock=>$haydock, nh=>10);
    my $flv=$flo->field;
    my $fla=1/$ea;
    my $flb=1/$eb;
    my $fproml=$fla*(1-$gl->f)+$flb*($gl->f);
    ($fla, $flb)=map {$_/$fproml} ($fla, $flb);
    my $flx=pdl(($fla*(1-$B)+$flb*$B), 0,0)->mv(-1,0);
    ok(Cagree($flv, $flx), "1D long field");
}

{
    #View 2D from 1D superlattice. Long wavelength transverse case
    my $Bt=zeroes(1,1,11)->yvals<5; #2D flat system
    my $epsilont=$ea*(1-$Bt)+$eb*$Bt;
    my $mu=ones(1,1,11)->r2C;
    my $gt=Photonic::Geometry::FromB->new(B=>$Bt); #trans
    my $mt=Photonic::WEM::S::Metric->new(mu=>$mu, geometry=>$gt, epsilon=>pdl(1),
					wavenumber=>pdl(0.001), wavevector=>pdl([0,0,0.0001]));
    my $nt=Photonic::WEM::S::Haydock->new(
	metric=>$mt, nh=>10, keepStates=>1, polarization=>pdl([1,0,0])->r2C,
	epsilon=>$epsilont);
    my $fto=Photonic::WEM::S::Field->new(haydock=>$nt, nh=>10);
    my $ftv=$fto->field;
    my $ftx=r2C(pdl [1, 0,0]);
    ok(Cagree($ftv, $ftx), "1D trans field");
}

#Check rawfields for simple 1D system.
# Longitudinal case
{
    my $B=zeroes(11,1,1)->xvals<5; #1D system
    my $epsilon=$ea*(1-$B)+$eb*$B;
    my $mu=ones(11,1,1)->r2C;
    my $gl=Photonic::Geometry::FromB->new(B=>$B); #long
    my $ml=Photonic::WEM::S::Metric->new(mu=>$mu, geometry=>$gl, epsilon=>pdl(1),
					wavenumber=>pdl(1), wavevector=>pdl([0.01,0,0]));
    my $haydock=Photonic::WEM::S::Haydock->new(
	metric=>$ml, nh=>10, keepStates=>1, polarization=>pdl([1,0,0])->r2C,
	epsilon=>$epsilon);
    my $flo=Photonic::WEM::S::Field->new(haydock=>$haydock, nh=>10);
    my $flv=$flo->rawfield;
    my $fla=1/$ea;
    my $flb=1/$eb;
    my $flx=pdl(-4*PI*($fla*(1-$B)+$flb*$B), 0,0)->mv(-1,0);
    ok(Cagree($flv, $flx), "1D long rawfield");
}

#View 2D from 1D superlattice. Long wavelength transverse case
{
    my $Bt=zeroes(1,1,11)->yvals<5; #2D flat system
    my $epsilont=$ea*(1-$Bt)+$eb*$Bt;
    my $mu=ones(1,1,11)->r2C;
    my $gt=Photonic::Geometry::FromB->new(B=>$Bt); #trans
    my $q=pdl(0.0001);
    my $k=pdl([0,0,0.0002]);
    my $mt=Photonic::WEM::S::Metric->new(mu=>$mu, geometry=>$gt, epsilon=>pdl(1),
					wavenumber=>$q, wavevector=>$k);
    my $nt=Photonic::WEM::S::Haydock->new(
	metric=>$mt, nh=>10, keepStates=>1, polarization=>pdl([1,0,0])->r2C,
	epsilon=>$epsilont);
    my $fto=Photonic::WEM::S::Field->new(haydock=>$nt, nh=>10);
    my $ftv=$fto->rawfield;
    my $f=$gt->f;
    my $epsM=(1-$f)*$ea+$f*$eb;
    my $ftx=-4*PI*r2C(pdl [1, 0,0])->dummy(3,11)*$q**2/($epsM*$q**2-$k->inner($k));
    ok(Cagree($ftv, $ftx), "1D trans rawfield");
}
