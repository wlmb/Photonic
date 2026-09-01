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
use Photonic::LE::NR2::Haydock;
use Photonic::LE::NR2::Field;
use Photonic::LE::NR2::SHP;
use Photonic::LE::NR2::SH;
use Photonic::LE::NR2::SHChiTensor;

use Test::More;
use lib 't/lib';
use TestUtils;

my $ea=1+2*i;
my $eb=3+4*i;
#Check field for simple 1D system
my $B=zeroes(11)->xvals<5; #1D system
my $gl=Photonic::Geometry::FromB->new(B=>$B, Direction0=>pdl([1])); #long
my $haydock=Photonic::LE::NR2::Haydock->new(geometry=>$gl, nh=>10, keepStates=>1);
my $flo=Photonic::LE::NR2::Field->new(haydock=>$haydock, nh=>10, epsA=>$ea, epsB=>$eb);
my $flv=$flo->field;
my $fle=$flo->epsL;
my $fla=1/$ea;
my $flb=1/$eb;
my $favl=$fla*(1-$gl->f)+$flb*($gl->f);
my $flex=1/$favl;
my ($flan, $flbn)=map {$_/$favl} ($fla, $flb);
my $flx=($flan*(1-$B)+$flbn*$B)->slice("*1");
ok(Cagree($flv, $flx), "1D long field") or diag "got: $flv\nexpected: $flx";
ok(Cagree($fle, $flex), "1D long response") or diag "got: $fle\nexpected: $flex";

#View 2D from 1D superlattice.

my $Bt=zeroes(1,11)->yvals<5; #2D flat system
my $gt=Photonic::Geometry::FromB->new(B=>$Bt, Direction0=>pdl([1,0])); #trans
my $nt=Photonic::LE::NR2::Haydock->new(geometry=>$gt, nh=>10, keepStates=>1);
my $fto=Photonic::LE::NR2::Field->new(haydock=>$nt, nh=>10, filter=>ones(1), epsA=>$ea, epsB=>$eb);
my $ftv=$fto->field;
my $fte=$fto->epsL;
my $ftx=pdl([1, 0])->r2C;
ok(Cagree($ftv, $ftx), "1D trans field") or diag "got: $ftv\nexpected: $ftx";;
my $favt=$ea*(1-$gt->f)+$eb*($gt->f);
ok(Cagree($fte, $favt), "1D trans response") or diag "got: $fte\nexpected: $favt";

# check raw fields
# Longitudinal
# 1D
my $flv_raw=$flo->rawfield;
my $flx_raw=-4*PI*($fla*(1-$B)+$flb*$B)->dummy(0);
ok(Cagree($flv_raw, $flx_raw), "1D long rawfield");
# 2D
my $ftv_raw=$fto->rawfield;
my $ftx_raw=-(4*PI+0*i)/($ea*(1-$gt->f)+$eb*$gt->f)*pdl([1,0]);
ok(Cagree($ftv_raw, $ftx_raw), "1D trans field") or diag "got: $ftv_raw\nexpected: $ftx_raw";

## NonRetarded Second Harmonic Polarization
my ($nA, $nB) = (0, 1); # vacuum, then anything as is normalised to dB
my $nrshp=Photonic::LE::NR2::SHP->new(
  haydock=>$nt, nh=>10, filter=>ones(1),
  densityA=>$nA, densityB=>$nB,
);
my $ea2=$ea; # linear permittivity at 2w. Arbitrary
my $eb2=2*$eb; # linear permittivity at 2w. Arbitrary
my $nrsh=Photonic::LE::NR2::SH->new(
  shp=>$nrshp, epsA1=>$ea, epsB1=>$eb, epsA2=>$ea2, epsB2=>$eb2,
  filterflag => 1
);
my $got=$nrsh->dipolar;
my $expected = pdl(0);
ok(Cagree($got, $expected), "dipolar") or diag "got: $got\nexpected: $expected";
# quadrupolar contrib
$got=$nrsh->quadrupolar;
$expected = pdl(0);
ok(Cagree($got, $expected), "quadrupolar") or diag "got: $got\nexpected: $expected";
# SH ext. polarization in reciprocal space
$got=$nrsh->external_G;
$expected = pdl(0);
ok(Cagree($got, $expected), "external_G") or diag "got: $got\nexpected: $expected";
# SH ext. longitudinal polarization in reciprocal space
$got=$nrsh->externalL_G;
$expected = pdl(0);
ok(Cagree($got, $expected), "externalL_G") or diag "got: $got\nexpected: $expected";
# SH ext. longitudinal polarization proj. in recip. space
$got=$nrsh->externalVecL_G;
$expected = pdl(0);
ok(Cagree($got, $expected, 1e-15), "externalVecL_G") or diag "got: $got\nexpected: $expected";
# SH ext. longitudinal polarization proj. in real space'
$got=$nrsh->externalVecL;
$expected = pdl(0);
ok(Cagree($got, $expected, 1e-15), "externalVecL") or diag "got: $got\nexpected: $expected";
# SH ext. longitudinal polarization in Haydock representation
$got=$nrsh->externalL_n;
$expected = pdl(0);
ok(Cagree($got, $expected), "externalL_n") or diag "got: $got\nexpected: $expected";
# SH self consistent longitudinal polarization in Haydock representation
$got=$nrsh->selfConsistentL_G;
$expected = pdl(0);
ok(Cagree($got, $expected, 1e-15), "selfConsistentL_G") or diag "got: $got\nexpected: $expected";
# SH self consistent longitudinal polarization components in reciprocal space
$got=$nrsh->selfConsistentVecL_G;
$expected = pdl(0);
ok(Cagree($got, $expected, 1e-15), "selfConsistentVecL_G") or diag "got: $got\nexpected: $expected";
# SH self consistent longitudinal polarization vector field in real space
$got=$nrsh->selfConsistentVecL;
$expected = pdl(0);
ok(Cagree($got, $expected), "selfConsistentVecL") or diag "got: $got\nexpected: $expected";
# Linear "atomic" polarizability
$got=$nrsh->alpha1;
$expected = ($eb-1)/(4*PI*$nB);
ok(Cagree($got, $expected), "alpha1") or diag "got: $got\nexpected: $expected";
# SH "atomic" polarizability
$got=$nrsh->alpha2;
$expected = ($eb2-1)/(4*PI*$nB);
ok(Cagree($got, $expected), "alpha2") or diag "got: $got\nexpected: $expected";
# longitudinal field at second harmonic
$got=$nrsh->field2;
$expected = pdl([1,0]);
ok(Cagree($got, $expected), "field2") or diag "got: $got\nexpected: $expected";
# SH self consistent total polarization vector field in real space
$got=$nrsh->P2;
$expected = pdl(0);
ok(Cagree($got, $expected), "P2") or diag "got: $got\nexpected: $expected";
# Spectral variable at fundamental
$got=$nrsh->u1;
$expected = 1/(1-$eb/$ea);
ok(Cagree($got, $expected), "u1") or diag "got: $got\nexpected: $expected";
# second harmonic susceptibility tensor
my $chi=Photonic::LE::NR2::SHChiTensor->new(
  geometry=>$gl,
  densityA=>$nA, densityB=>$nB, nhf=>10, nh=>10,
  epsA1=>$ea, epsB1=>$eb, epsA2=>$ea*$ea, epsB2=>$eb*$eb,
  keepStates=>1,
);
# non retarded SH polarizations and fields
$got = Photonic::LE::NR2::SH->new(
  shp=>$chi->nrshp->[0],
  epsA1=>$ea, epsB1=>$eb,epsA2=>$ea*$ea, epsB2=>$eb*$eb,
  filterflag=>0
)->P2;
$expected = pdl(<<'EOF');
[
 [ 0.0070072648+0.0097445079i ]
 [ -0.0019177555-0.002666887i ]
 [ -2.7383241e-17+8.8425853e-18i ]
 [ 0.0019177555+0.002666887i ]
 [ -0.0070072648-0.0097445079i ]
 [  0.014991596-0.011757369i ]
 [ -0.0083074514+0.0065152354i ]
 [  0.0067611307-0.0053025117i ]
 [ -0.0067611307+0.0053025117i ]
 [  0.0083074514-0.0065152354i ]
 [ -0.014991596+0.011757369i ]
]
EOF
ok(Cagree($got, $expected), "SHChiTensor SHPs P2") or diag "got: $got\nexpected: $expected";
# different chis
$got = $chi->evaluate;
$expected = pdl(<<'EOF');
[ [ [ 0 ] ] ]
EOF
ok(Cagree($got, $expected), "P2") or diag "got: $got\nexpected: $expected";
# chi P2
$got = $chi->evaluate(kind => 'f', mask => pdl(1));
$expected = pdl(<<'EOF');
[ [ [ 0 ] ] ]
EOF
ok(Cagree($got, $expected), "P2") or diag "got: $got\nexpected: $expected";

# chi selfConsistentVecL
$got = $chi->evaluate(kind => 'l');
$expected = pdl(<<'EOF');
[ [ [ 0 ] ] ]
EOF
ok(Cagree($got, $expected), "selfConsistentVecL") or diag "got: $got\nexpected: $expected";
# chi P2LMCalt
$got = $chi->evaluate(kind => 'a');
$expected = pdl(<<'EOF');
[ [ [ 0 ] ] ]
EOF
ok(Cagree($got, $expected), "P2LMCalt") or diag "got: $got\nexpected: $expected";
# chi dipolar
$got = $chi->evaluate(kind => 'd');
$expected = pdl(<<'EOF');
[ [ [ 0 ] ] ]
EOF
ok(Cagree($got, $expected), "dipolar") or diag "got: $got\nexpected: $expected";
# chi quadrupolar
$got = $chi->evaluate(kind => 'q');
$expected = pdl(<<'EOF');
[ [ [ 0 ] ] ]
EOF
ok(Cagree($got, $expected), "quadrupolar") or diag "got: $got\nexpected: $expected";
# chi external
$got = $chi->evaluate(kind => 'e');
$expected = pdl(<<'EOF');
[ [ [ 0 ] ] ]
EOF
ok(Cagree($got, $expected), "external") or diag "got: $got\nexpected: $expected";
# chi externalVecL
$got = $chi->evaluate(kind => 'el');
$expected = pdl(<<'EOF');
[ [ [ 0 ] ] ]
EOF
ok(Cagree($got, $expected), "externalVecL") or diag "got: $got\nexpected: $expected";

done_testing;
