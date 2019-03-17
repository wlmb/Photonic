use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;

use lib qw(/home/gortiz/Metas/Photonic/PhotonicWLM/lib);
use Photonic::Geometry::FromB;
use Photonic::WE::R2::Metric;
use Photonic::WE::R2::AllH;

use Machine::Epsilon;
use List::Util;

use Test::More tests => 9;

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
my @N=(10,20,30,40,50);
foreach(@N){
diag("check reorthogonalize with square array, N=" .($_+1)); 
diag("Particle" . " ");
my $Bs=zeroes(($_+1),($_+1))->rvals<($_+1)/3;
my $gs=Photonic::Geometry::FromB->new(B=>$Bs);
my $ms=Photonic::WE::R2::Metric->new(geometry=>$gs, epsilon=>pdl(2),
   wavenumber=>pdl(.01), wavevector=>pdl([.001,0]));
my $als=Photonic::WE::R2::AllH
    ->new(metric=>$ms, polarization=>r2C(pdl([0,1])), nh=>2*($_+1)*($_+1),
    reorthogonalize=>1, accuracy=>machine_epsilon(),
    noise=>1e0*machine_epsilon(), normOp=>1e0, smallH=>1e-7);   
$als->run;
ok($als->iteration <= ($_+1)*($_+1), "No more iterations than dimensions");
diag("Actual iterations: " . $als->iteration 
     . " Actual orthogonalizations: ", $als->orthogonalizations);
diag("Hole" . " ");
#check reorthogonalize again with square array
my $B1s=zeroes(($_+1),($_+1))->rvals>($_+1)/3; #($_+1),($_+1), 2
my $g1s=Photonic::Geometry::FromB->new(B=>$B1s, L=>pdl(1,1));
my $m1s=Photonic::WE::R2::Metric->new(geometry=>$g1s, epsilon=>pdl(2),
   wavenumber=>pdl(0.01), wavevector=>pdl([0.001,0]));
my $al1s=Photonic::WE::R2::AllH
    ->new(metric=>$m1s, polarization=>r2C(pdl([0,1])), nh=>2*($_+1)*($_+1),
    reorthogonalize=>1, accuracy=>machine_epsilon(),
    noise=>1e1*machine_epsilon(), normOp=>1e0, smallH=>1e-7);   
$al1s->run;
ok($al1s->iteration <= 2*($_+1)*($_+1), "No more iterations than dimensions");
diag("Actual iterations: " . $al1s->iteration 
     . " Actual orthogonalizations: ", $al1s->orthogonalizations);

diag("check reorthogonalize again with square array even number, N=" . $_);
my $B2s=zeroes($_,$_)->rvals>$_/2; 
my $g2s=Photonic::Geometry::FromB->new(B=>$B2s, L=>pdl(1,1));
my $m2s=Photonic::WE::R2::Metric->new(geometry=>$g2s, epsilon=>pdl(2),
   wavenumber=>pdl(3.6), wavevector=>pdl([1.01*3.6,0]));
my $al2s=Photonic::WE::R2::AllH
    ->new(metric=>$m2s, polarization=>r2C(pdl([0,1])), nh=>2*$_*$_,
    reorthogonalize=>1, accuracy=>machine_epsilon(),
    noise=>1e0*machine_epsilon(), normOp=>1e0, smallH=>1e-7, use_mask=>1);   
$al2s->run;
ok($al2s->iteration <= 2*$_*$_, "No more iterations than dimensions");
diag("Actual iterations: " . $al2s->iteration 
     . " Actual orthogonalizations: ", $al2s->orthogonalizations);
diag(" ");
}