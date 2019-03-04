use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromB;
use Photonic::WE::R2::Metric;
use Photonic::WE::R2::AllH;

use Machine::Epsilon;
use List::Util;

use Test::More tests => 8;

#my $pi=4*atan2(1,1);

sub agree {    
    my $a=shift;
    my $b=shift//0;
    return (($a-$b)*($a-$b))->sum<=1e-7;
}

#Check haydock coefficients for simple 1D system
my $B=zeroes(11)->xvals<5; #1D system
my $g=Photonic::Geometry::FromB->new(B=>$B); #long
my $m=Photonic::WE::R2::Metric->new(geometry=>$g, epsilon=>pdl(1),
   wavenumber=>pdl(1), wavevector=>pdl([0.01]));
my $a=Photonic::WE::R2::AllH->new(metric=>$m,
   polarization=>pdl([1])->r2C, nh=>10); 
$a->run;
my $as=$a->as;
my $bs=$a->bs;
my $b2s=$a->b2s;
ok(agree(pdl($a->iteration), 2), "Number of iterations 1D longitudinal");
ok(agree(pdl($b2s->[0]), 1), "1D L b_0^2");
ok(agree(pdl($b2s->[1]), $g->f*(1-$g->f)), "1D L b_1^2");
ok(agree(pdl($as->[0]), $g->f), "1D L a_0");
ok(agree(pdl($as->[1]), 1-$g->f), "1D L a_1");
ok(agree(pdl($b2s), pdl($bs)**2), "1D L b2==b^2");
#check reorthogonalize with square array
my $Bs=zeroes(15,15)->rvals<5;
my $gs=Photonic::Geometry::FromB->new(B=>$Bs);
my $ms=Photonic::WE::R2::Metric->new(geometry=>$gs, epsilon=>pdl(1),
   wavenumber=>pdl(.01), wavevector=>pdl([.001,0]));
my $als=Photonic::WE::R2::AllH
    ->new(metric=>$ms, polarization=>r2C(pdl([0,1])), nh=>2*15*15,
    reorthogonalize=>1, accuracy=>machine_epsilon(),
    noise=>1e3*machine_epsilon(), normOp=>1e0, smallH=>1e-7);   
$als->run;
ok($als->iteration <= 15*15, "No more iterations than dimensions");
diag("Actual iterations: " . $als->iteration 
     . " Actual orthogonalizations: ", $als->orthogonalizations);
my $st=$als->states;
#check reorthogonalize again with square array
my $B1s=zeroes(21,21)->rvals>2; #21,21, 2
my $g1s=Photonic::Geometry::FromB->new(B=>$B1s, L=>pdl(1,1));
my $m1s=Photonic::WE::R2::Metric->new(geometry=>$g1s, epsilon=>pdl(10),
   wavenumber=>pdl(3.6), wavevector=>pdl([1.01*3.6,0]));
my $al1s=Photonic::WE::R2::AllH
    ->new(metric=>$m1s, polarization=>r2C(pdl([0,1])), nh=>3*21*21,
    reorthogonalize=>1, accuracy=>machine_epsilon(),
    noise=>1e3*machine_epsilon(), normOp=>1e0, smallH=>1e-7);   
$al1s->run;
ok($al1s->iteration <= 2*21*21, "No more iterations than dimensions");
diag("Actual iterations: " . $al1s->iteration 
     . " Actual orthogonalizations: ", $al1s->orthogonalizations);
#foreach(@$st){
#    my $pr=$als->innerProduct($_, $st->[0]);
#    print "$pr\n";
#}
	