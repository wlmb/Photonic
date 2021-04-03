use strict;
use warnings;
use PDL;
use PDL::Complex;
use Photonic::Utils;
Photonic::Utils->import(@Photonic::Utils::EXPORT_OK);
use Test::More;

# not yet covered: GtoR, RtoG, tile, vectors2Dlist

my $x = zeroes(2, 11)->complex;
$x->slice(':,0') .= 1+0*i;
ok approx(HProd($x, $x), r2C(1));
ok approx(EProd($x, $x), r2C(1));

$x = zeroes(2, 1, 11)->complex;
$x->slice(':,:,0') .= 1+0*i;
ok approx(MHProd($x, $x, ones(1, 1, 11)), r2C(1));

$x = zeroes(2, 2, 11)->complex;
$x->slice(':,:,0') .= 1/sqrt(2)+0*i;
ok approx(SProd($x, $x), r2C(1));

$x = zeroes(2, 1, 2, 11)->complex;
$x->slice(':,:,:,0') .= 1/sqrt(2)+0*i;
ok approx(VSProd($x, $x), r2C(1));

my $got = lentzCF(
  [21/11 + i*32/11, 23/11 + i*34/11],
  [r2C(-1), -8.7603535536828499e-17 - 1.98347107438017*i],
  2,
  1e-7,
);
my $expected = 1.4688427299703299 + 2.6112759643916901*i;
ok approx($got, $expected) or diag "got: $got, expected $expected";

done_testing;
