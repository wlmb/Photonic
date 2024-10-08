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
use Photonic::Geometry::FromEpsilonTensor;
use Photonic::LE::ST::EpsTensor;

use Test::More tests => 7;
use lib 't/lib';
use TestUtils;

my $ea=1+2*i;
my $eb=3+4*i;
#Check haydock coefficients for simple 1D system
my $B=zeroes(11)->xvals<5; #1D system
my $f=$B->sum/$B->nelem;
my $epsilon=($ea*(1-$B)+$eb*$B)->slice('*1,*1'); # *identity(1);
my $gl=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilon); #long
my $elo=Photonic::LE::ST::EpsTensor->new(geometry=>$gl, nh=>10);
my $elv=$elo->epsTensor;
my $elx=1/((1-$f)/$ea+$f/$eb);
ok(Cagree($elv, $elx), "1D long epsilon");
is($elo->converged,1, "Converged");
#View 2D from 1D superlattice.
my $Bt=zeroes(1,11)->yvals<5; #2D flat system
my $epsilont=($ea*(1-$Bt)+$eb*$Bt)->slice('*1,*1')*identity(2);
my $gt=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilont); #trans
my $eto=Photonic::LE::ST::EpsTensor->new(geometry=>$gt, nh=>10);
my $etv=$eto->epsTensor;
my $etx=(1-$f)*$ea+$f*$eb;
my $etenx=pdl([$etx, r2C(0)],[r2C(0), $elx]);
ok(Cagree($etv, $etenx), "1D trans epsilon");
is($eto->converged,1, "Converged");
#Extend 1D superlattice into 4D (why not?)
my $Bt4=zeroes(11,1,1,1)->xvals<5; #2D flat system
my $epsilont4=($ea*(1-$Bt4)+$eb*$Bt4)->slice('*1,*1')*identity(4);
my $gt4=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilont4); #trans
my $eto4=Photonic::LE::ST::EpsTensor->new(geometry=>$gt4, nh=>10);
my $etv4=$eto4->epsTensor;
my $etenx4=pdl([
    [$elx,   r2C(0),  r2C(0), r2C(0)],
    [r2C(0),  $etx,  r2C(0), r2C(0)],
    [r2C(0), r2C(0), $etx,   r2C(0)],
    [r2C(0), r2C(0), r2C(0), $etx  ]
    ]);
ok(Cagree($etv4, $etenx4), "4D trans epsilon");
is($eto4->converged,1, "Converged");



#Keller
my $Nk=6;
my $Bk=zeroes(2*$Nk,2*$Nk);
$Bk=((($Bk->xvals<$Nk) & ($Bk->yvals<$Nk))
   | (($Bk->xvals>=$Nk) & ($Bk->yvals>=$Nk)));
my $epsilonk=($ea*(1-$Bk)+$eb*$Bk)->slice('*1,*1')*identity(2);
my $gk=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilonk); #
my $eko=Photonic::LE::ST::EpsTensor->new(
    geometry=>$gk, nh=>1000, reorthogonalize=>1, use_mask=>1);
my $etva=$eko->epsTensor;
my $epsilonkk=($eb*(1-$Bk)+$ea*$Bk)->slice('*1,*1')*identity(2);
my $gkk=Photonic::Geometry::FromEpsilonTensor->new(epsilon=>$epsilonkk); #
my $ekko=Photonic::LE::ST::EpsTensor->new(
    geometry=>$gkk, nh=>1000, reorthogonalize=>1, use_mask=>1);
my $etvb=$eko->epsTensor;
my $etvr=zeroes(2,2)->r2C;
$etvr->((0),(0)).= $etvb->((1),(1));
$etvr->((0),(1)).=-$etvb->((1),(0));
$etvr->((1),(0)).=-$etvb->((0),(1));
$etvr->((1),(1)).= $etvb->((0),(0));
my $etvar=($etva->(*1)*$etvr->(:,:,*1))->mv(2,1)->sumover;
ok(Cagree($etvar,$ea*$eb*identity(2), 1e-3), "Keller");
