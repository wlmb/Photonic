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
use Photonic::WE::R2::Haydock;
use Photonic::WE::R2::Metric;
use Photonic::WE::R2::Field;

use Test::More tests => 4;
use lib 't/lib';
use TestUtils;

my $ea=pdl(2);
my $eb=3+4*i;

#Check haydock coefficients for simple 1D system.
{
    #Longitudinal case
    my $B=zeroes(11)->xvals<5; #1D system
    my $gl=Photonic::Geometry::FromB->new(B=>$B);
    my $ml=Photonic::WE::R2::Metric->new(geometry=>$gl, epsilon=>$ea,
       wavenumber=>pdl(1), wavevector=>pdl([0.01]));
    my $haydock=Photonic::WE::R2::Haydock->new(metric=>$ml, nh=>10, keepStates=>1,
				       polarization=>pdl([1])->r2C);
    my $flo=Photonic::WE::R2::Field->new(haydock=>$haydock, nh=>10, epsB=>$eb);
    my $flv=$flo->field;
    my $fla=1/$ea;
    my $flb=1/$eb;
    my $fproml=$fla*(1-$gl->f)+$flb*($gl->f);
    ($fla, $flb)=map {$_/$fproml} ($fla, $flb);
    my $flx=($fla*(1-$B)+$flb*$B)->transpose;
    ok(Cagree($flv, $flx), "1D long field");
}

{
    #View 2D from 1D superlattice. Long wavelength transverse case
    my $Bt=zeroes(1,11)->yvals<5; #2D flat system
    my $gt=Photonic::Geometry::FromB->new(B=>$Bt); #trans
    my $mt=Photonic::WE::R2::Metric->new(geometry=>$gt, epsilon=>pdl(1),
       wavenumber=>pdl(0.001), wavevector=>pdl([0,0.0001]));
    my $nt=Photonic::WE::R2::Haydock->new(metric=>$mt, nh=>10, keepStates=>1,
				       polarization=>pdl([1,0])->r2C);
    my $fto=Photonic::WE::R2::Field->new(haydock=>$nt, nh=>10, epsB=>$eb);
    my $ftv=$fto->field;
    my $ftx=r2C(pdl [1, 0]);
    ok(Cagree($ftv, $ftx), "1D trans field");
}


#Check rawfields for simple 1D system.
# Longitudinal case
{
    my $B=zeroes(11)->xvals<5; #1D system
    my $gl=Photonic::Geometry::FromB->new(B=>$B);
    my $ml=Photonic::WE::R2::Metric->new(geometry=>$gl, epsilon=>$ea,
       wavenumber=>pdl(1), wavevector=>pdl([0.01]));
    my $haydock=Photonic::WE::R2::Haydock->new(metric=>$ml, nh=>10, keepStates=>1,
				       polarization=>pdl([1])->r2C);
    my $flo=Photonic::WE::R2::Field->new(haydock=>$haydock, nh=>10, epsB=>$eb);
    my $flv=$flo->rawfield;
    my $fla=1/$ea;
    my $flb=1/$eb;
    my $flx=-4*PI*($fla*(1-$B)+$flb*$B)->dummy(0);
    ok(Cagree($flv, $flx), "1D long rawfield");
}

#View 2D from 1D superlattice. Long wavelength transverse case
{
    my $Bt=zeroes(1,11)->yvals<5; #2D flat system
    my $gt=Photonic::Geometry::FromB->new(B=>$Bt); #trans
    my $q=pdl(0.000000000001);
    my $k=pdl(0.0000000001);
    my $pol=pdl([1,0]);
    my $mt=Photonic::WE::R2::Metric->new(geometry=>$gt, epsilon=>$ea,
       wavenumber=>$q, wavevector=>pdl([0,$k]));
    my $nt=Photonic::WE::R2::Haydock->new(metric=>$mt, nh=>10, keepStates=>1,
				       polarization=>$pol->r2C);
    my $fto=Photonic::WE::R2::Field->new(haydock=>$nt, nh=>10, epsB=>$eb);
    my $ftv=$fto->rawfield;
    my $f=$gt->f;
    my $epsM=(1-$f)*$ea+$f*$eb;
    my $ftx=-4*PI*$pol->r2C->dummy(2,11)*$q**2/($epsM*$q**2-$k**2);
    ok(Cagree($ftv, $ftx), "1D trans rawfield");
}
