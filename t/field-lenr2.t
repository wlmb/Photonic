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
use PDL::Complex;
use Photonic::LE::NR2::AllH;
use Photonic::LE::NR2::Field;
use Photonic::LE::NR2::SHP;
use Photonic::LE::NR2::SH;

use Test::More;
use lib 't/lib';
use TestUtils;

my $ea=1+2*i;
my $eb=3+4*i;
#Check field for simple 1D system
my $B=zeroes(11)->xvals<5; #1D system
my $gl=Photonic::Geometry::FromB->new(B=>$B, Direction0=>pdl([1])); #long
my $nr=Photonic::LE::NR2::AllH->new(geometry=>$gl, nh=>10, keepStates=>1);
my $flo=Photonic::LE::NR2::Field->new(nr=>$nr, nh=>10);
my $flv=$flo->evaluate($ea, $eb);
my $fla=1/$ea;
my $flb=1/$eb;
my $fproml=$fla*(1-$gl->f)+$flb*($gl->f);
($fla, $flb)=map {$_/$fproml} ($fla, $flb);
my $flx=pdl([$fla*(1-$B)+$flb*$B])->complex->mv(1,-1);
ok(Cagree($flv, $flx), "1D long field") or diag "got: $flv\nexpected: $flx";

#View 2D from 1D superlattice.
my $Bt=zeroes(1,11)->yvals<5; #2D flat system
my $gt=Photonic::Geometry::FromB->new(B=>$Bt, Direction0=>pdl([1,0])); #trans
my $nt=Photonic::LE::NR2::AllH->new(geometry=>$gt, nh=>10, keepStates=>1);
my $fto=Photonic::LE::NR2::Field->new(nr=>$nt, nh=>10);
my $ftv=$fto->evaluate($ea, $eb);
my $ftx=pdl(r2C(1), r2C(0))->complex;
ok(Cagree($ftv, $ftx), "1D trans field");

my ($dA, $dB) = (0, 1); # vacuum, then anything as is normalised to dB
my $nrshp=Photonic::LE::NR2::SHP->new(nrf=>$fto, densityA=>$dA, densityB=>$dB);
my $nrsh=Photonic::LE::NR2::SH->new(
  shp=>$nrshp, epsA1=>$ea, epsB1=>$eb,
  epsA2=>$ea*$ea, epsB2=>$eb*$eb, # also very arbitrary!
);
my $got=$nrsh->selfConsistentL_G;
my $expected = pdl(<<'EOF')->complex;
[
 [ [0 0] ]
 [ [2.3385878e-18  1.910126e-18] ]
 [ [2.7283262e-18  2.436351e-19] ]
 [ [2.3719287e-18 6.4418678e-20] ]
 [ [4.2223945e-18 -1.807201e-18] ]
 [ [ 9.5067349e-19 -6.5454932e-18] ]
 [ [5.5266405e-18 3.6340151e-18] ]
 [ [ 4.1237233e-18 -2.0133211e-18] ]
 [ [ 1.4476152e-18 -1.8171605e-18] ]
 [ [ 1.6228705e-18 -2.2408099e-18] ]
 [ [-4.4470906e-20 -3.1472783e-18] ]
]
EOF
ok(Cagree($got, $expected), "self-consistent L_G") or diag "got: $got\nexpected: $expected";
$got=$nrsh->alpha1;
$expected = pdl(<<'EOF')->complex;
[
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
]
EOF
ok(Cagree($got, $expected), "alpha1") or diag "got: $got\nexpected: $expected";
$got=$nrsh->alpha1;
$expected = pdl(<<'EOF')->complex;
[
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
 [ [0.15915494 0.31830989] ]
]
EOF
ok(Cagree($got, $expected), "alpha2") or diag "got: $got\nexpected: $expected";
$got=$nrsh->dipolar;
$expected = pdl(<<'EOF')->complex;
[
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [-7.9823116e-17  1.1403302e-17] ] ]
 [ [ [0 0] [ 4.1596509e-17 -5.9423584e-18] ] ]
 [ [ [0 0] [-2.9756908e-17  4.2509868e-18] ] ]
 [ [ [0 0] [ 2.4722933e-17 -3.5318475e-18] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
]
EOF
ok(Cagree($got, $expected), "dipolar") or diag "got: $got\nexpected: $expected";
$got=$nrsh->quadrupolar;
$expected = pdl(<<'EOF')->complex;
[
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
]
EOF
ok(Cagree($got, $expected), "quadrupolar") or diag "got: $got\nexpected: $expected";
$got=$nrsh->externalVecL;
$expected = pdl(<<'EOF')->complex;
[
 [ [ [0 0] [ 3.9327802e-18 -5.6182574e-19] ] ]
 [ [ [0 0] [-7.5890335e-17  1.0841476e-17] ] ]
 [ [ [0 0] [ 4.5529289e-17 -6.5041841e-18] ] ]
 [ [ [0 0] [-2.5824128e-17  3.6891611e-18] ] ]
 [ [ [0 0] [ 2.8655713e-17 -4.0936733e-18] ] ]
 [ [ [0 0] [ 3.9327802e-18 -5.6182574e-19] ] ]
 [ [ [0 0] [ 3.9327802e-18 -5.6182574e-19] ] ]
 [ [ [0 0] [ 3.9327802e-18 -5.6182574e-19] ] ]
 [ [ [0 0] [ 3.9327802e-18 -5.6182574e-19] ] ]
 [ [ [0 0] [ 3.9327802e-18 -5.6182574e-19] ] ]
 [ [ [0 0] [ 3.9327802e-18 -5.6182574e-19] ] ]
]
EOF
ok(Cagree($got, $expected), "externalVecL") or diag "got: $got\nexpected: $expected";
$got=$nrsh->field2;
$expected = pdl(<<'EOF')->complex;
[-0.6756993 -0.10576923]
EOF
ok(Cagree($got, $expected), "field2") or diag "got: $got\nexpected: $expected";
$got=$nrsh->P2;
$expected = pdl(<<'EOF')->complex;
[
 [ [ [0 0] [-1.2467976e-18 -2.5931346e-19] ] ]
 [ [ [0 0] [8.5108064e-20 2.6781772e-18] ] ]
 [ [ [0 0] [-1.9408651e-18  -1.790065e-18] ] ]
 [ [ [0 0] [-7.5028237e-19  8.3574075e-19] ] ]
 [ [ [0 0] [-1.6593174e-18 -1.1691174e-18] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
 [ [ [0 0] [0 0] ] ]
]
EOF
ok(Cagree($got, $expected), "P2") or diag "got: $got\nexpected: $expected";
$got=$nrsh->u1;
$expected = pdl(<<'EOF')->complex;
[-0.22115385 -0.10576923]
EOF
ok(Cagree($got, $expected), "u1") or diag "got: $got\nexpected: $expected";

done_testing;
