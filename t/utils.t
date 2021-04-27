use strict;
use warnings;
use PDL;
use PDL::Complex;
use Photonic::Utils;
Photonic::Utils->import(@Photonic::Utils::EXPORT_OK);
use Test::More;

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
  pdl([21/11 + 32*i/11, 23/11 + 34*i/11])->cplx,
  pdl([r2C(-1), -8.7603535536828499e-17 - 1.98347107438017*i])->cplx,
  2,
  1e-7,
);
my $expected = 1.4688427299703299 + 2.6112759643916901*i;
ok approx($got, $expected) or diag "got: $got, expected $expected";

$expected = -1.2045454545454499 - 0.79545454545454497*i;
$got = lentzCF(
  $expected,
  -1 - 0.247933884297521*i,
  1,
  1e-7,
);
ok approx($got, $expected) or diag "got: $got, expected $expected";

$expected = pdl <<'EOF';
[
 [
  [
   [     1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115      1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115      1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115]
   [   -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402    -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402    -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402]
   [  1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16   1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16   1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16]
   [    0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402     0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402     0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402]
   [    -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115     -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115     -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115]
   [     1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115      1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115      1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115]
   [   -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402    -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402    -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402]
   [  1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16   1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16   1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16]
   [    0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402     0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402     0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402]
   [    -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115     -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115     -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115]
   [     1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115      1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115      1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115]
   [   -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402    -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402    -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402]
   [  1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16   1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16   1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16]
   [    0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402     0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402     0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402]
   [    -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115     -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115     -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115]
  ]
  [
   [  0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639   0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639   0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639]
   [  -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357   -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357   -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357]
   [ 6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17  6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17  6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17]
   [   0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357    0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357    0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357]
   [ -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639  -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639  -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639]
   [  0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639   0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639   0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639]
   [  -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357   -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357   -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357]
   [ 6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17  6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17  6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17]
   [   0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357    0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357    0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357]
   [ -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639  -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639  -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639]
   [  0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639   0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639   0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639]
   [  -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357   -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357   -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357]
   [ 6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17  6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17  6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17]
   [   0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357    0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357    0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357]
   [ -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639  -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639  -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639]
  ]
 ]
 [
  [
   [  2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937]
   [  1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453]
   [ -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781  -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781  -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781]
   [  1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453]
   [  2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937]
   [  2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937]
   [  1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453]
   [ -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781  -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781  -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781]
   [  1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453]
   [  2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937]
   [  2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937]
   [  1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453]
   [ -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781  -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781  -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781]
   [  1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453   1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453]
   [  2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937   2.8779937   3.7618536   3.0752377   3.7618536   2.8779937]
  ]
  [
   [ 0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445]
   [ 0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552]
   [ -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999  -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999  -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999]
   [ 0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552]
   [ 0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445]
   [ 0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445]
   [ 0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552]
   [ -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999  -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999  -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999]
   [ 0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552]
   [ 0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445]
   [ 0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445]
   [ 0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552]
   [ -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999  -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999  -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999]
   [ 0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552  0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552]
   [ 0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445  0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445]
  ]
 ]
]
EOF
$got = tile(pdl(<<'EOF'), 3, 3);
[
 [
  [
   [     1.3495115   -0.023957563 -5.9880532e-16    0.023957563     -1.3495115]
   [   -0.10812402    -0.85100218  1.4240728e-17     0.85100218     0.10812402]
   [  1.110223e-16  4.9960036e-16 -2.9610637e-16  2.7755576e-16 -5.5511151e-16]
   [    0.10812402     0.85100218  1.5850844e-16    -0.85100218    -0.10812402]
   [    -1.3495115    0.023957563 -1.2080698e-16   -0.023957563      1.3495115]
  ]
  [
   [  0.0081485639   0.0081573366  1.1177956e-17  -0.0081573366  -0.0081485639]
   [  -0.009254357   0.0041060736  5.5476509e-18  -0.0041060736    0.009254357]
   [ 6.9388939e-18  1.2750218e-16  1.2143128e-17  -1.405126e-16  2.2551405e-17]
   [   0.009254357  -0.0041060736 -4.1711898e-18   0.0041060736   -0.009254357]
   [ -0.0081485639  -0.0081573366  -2.084947e-17   0.0081573366   0.0081485639]
  ]
 ]
 [
  [
   [  2.8779937   3.7618536   3.0752377   3.7618536   2.8779937]
   [  1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453]
   [ -3.7800781  -1.4194227 -0.64872161  -1.4194227  -3.7800781]
   [  1.5120453 -0.55214223 -0.25087687 -0.55214223   1.5120453]
   [  2.8779937   3.7618536   3.0752377   3.7618536   2.8779937]
  ]
  [
   [ 0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445]
   [ 0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552]
   [ -0.031969999 -0.0079121337  -0.015230545 -0.0079121337  -0.031969999]
   [ 0.0067999552  -0.013824409  -0.016607758  -0.013824409  0.0067999552]
   [ 0.0091850445   0.017780476   0.024223031   0.017780476  0.0091850445]
  ]
 ]
]
EOF
ok all(approx($got, $expected)) or diag "got: $got, expected $expected";

$expected = pdl(<<'EOF')->complex;
[
 [ [ [0 0] [-4.3260582e-17  6.1800832e-18] ] ]
 [ [ [0 0] [-5.9528635e-17  2.4920212e-17] ] ]
 [ [ [0 0] [-2.7186406e-17  6.2313278e-17] ] ]
 [ [ [0 0] [-1.2929662e-17  5.3169197e-17] ] ]
 [ [ [0 0] [1.3401775e-17   1.10965e-16] ] ]
 [ [ [0 0] [1.5419056e-16 6.9763677e-17] ] ]
 [ [ [0 0] [ 1.2848911e-16 -1.1014649e-16] ] ]
 [ [ [0 0] [-1.8204495e-17 -1.1027889e-16] ] ]
 [ [ [0 0] [ -2.729985e-17 -4.7422124e-17] ] ]
 [ [ [0 0] [-4.3546668e-17 -5.2208553e-17] ] ]
 [ [ [0 0] [-6.4125149e-17 -7.2553859e-18] ] ]
]
EOF
$got = RtoG(pdl(<<'EOF')->complex, 2, 1);
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
ok all(approx($got, $expected)) or diag "got: $got, expected $expected";

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
$got = GtoR(pdl(<<'EOF')->complex, 2, 1);
[
 [ [ [0 0] [0 0] ] ]
 [ [ [-0 0] [-5.9528635e-17  2.4920212e-17] ] ]
 [ [ [-0 0] [-2.7186406e-17  6.2313278e-17] ] ]
 [ [ [-0 0] [-1.2929662e-17  5.3169197e-17] ] ]
 [ [ [0  0] [1.3401775e-17   1.10965e-16] ] ]
 [ [ [0  0] [1.5419056e-16 6.9763677e-17] ] ]
 [ [ [-0 0] [ 1.2848911e-16 -1.1014649e-16] ] ]
 [ [ [0  0] [-1.8204495e-17 -1.1027889e-16] ] ]
 [ [ [0  0] [ -2.729985e-17 -4.7422124e-17] ] ]
 [ [ [0  0] [-4.3546668e-17 -5.2208553e-17] ] ]
 [ [ [0  0] [-6.4125149e-17 -7.2553859e-18] ] ]
]
EOF
ok all(approx($got, $expected)) or diag "got: $got, expected $expected";

$expected = pdl(<<'EOF');
[
  [0 5 10 0 5 10 0 5 10]
  [0 0 0 5 5 5 10 10 10]
  [0 0 0 0 0 0 0 0 0]
  [0 0 0 0 0 0 0 0 0]
]
EOF
$got = pdl(vectors2Dlist(pdl(<<'EOF'), 0, 5));
[
 [
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
 ]
 [
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
 ]
 [
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
 ]
 [
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
 ]
 [
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
 ]
 [
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
 ]
 [
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
 ]
 [
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
 ]
 [
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
 ]
 [
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
 ]
 [
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
  [    0.5751249    0.81806562]
  [    0.2656728    0.96406326]
  [2.2330872e-15             1]
  [   -0.2656728    0.96406326]
  [   -0.5751249    0.81806562]
 ]
 [
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
  [  -0.34115881    0.94000567]
  [ -0.030418434   -0.99953725]
  [-2.453355e-15            -1]
  [  0.030418434   -0.99953725]
  [   0.34115881    0.94000567]
 ]
 [
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
  [ 2.2231984e-15             -1]
  [ 7.4694778e-15             -1]
  [ 1.2043164e-15             -1]
  [-8.4014982e-15             -1]
  [-6.6911796e-16             -1]
 ]
 [
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
  [   0.34115881    0.94000567]
  [  0.030418434   -0.99953725]
  [8.8386295e-16            -1]
  [ -0.030418434   -0.99953725]
  [  -0.34115881    0.94000567]
 ]
 [
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
  [    -0.5751249     0.81806562]
  [    -0.2656728     0.96406326]
  [-1.9584407e-15              1]
  [     0.2656728     0.96406326]
  [     0.5751249     0.81806562]
 ]
]
EOF
ok all(approx($got, $expected)) or diag "got: $got, expected $expected";

my $data = pdl(<<'EOF')->complex;
[
 [    0.72727273              0]
 [      1733.403 -2.8318751e-28]
 [     3466.0788             -0]
]
EOF
$got = tensor($data, [
  pdl('[ [1 0 0] [0.5 1 0.5] [0 0 1] ]'),
  pdl('[1 2 3]'),
], 2);
$expected = pdl(<<'EOF')->complex;
[
 [
  [    0.72727273              0]
  [   -3.6365e-05 -2.8318751e-28]
 ]
 [
  [   -3.6365e-05 -2.8318751e-28]
  [     3466.0788              0]
 ]
]
EOF
ok all(approx($got, $expected)) or diag "got: $got, expected $expected";

$data = pdl(<<'EOF')->complex;
[
 [
  [    0.67272727    -0.10909091]
  [-1.3253896e-09  1.7631692e-09]
 ]
 [
  [-1.3253896e-09  1.7631692e-09]
  [      -153.825     -1737.5269]
 ]
]
EOF
$got = wave_operator($data, 2);
$expected = pdl(<<'EOF')->complex;
[
 [
  [    1.4483986    0.23487544]
  [1.1625849e-12 1.4461083e-12]
 ]
 [
  [ 1.1625849e-12  1.4461083e-12]
  [-5.0556058e-05  0.00057105487]
 ]
]
EOF
ok all(approx($got, $expected)) or diag "got: $got, expected $expected";

done_testing;
