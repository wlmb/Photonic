use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromB;
use Photonic::LE::NR2::AllH;
use Photonic::Utils qw(HProd);

use Machine::Epsilon;
use List::Util;

use Test::More tests => 12;

#my $pi=4*atan2(1,1);

sub agree {    
    my $a=shift;
    my $b=shift//0;
    return (($a-$b)*($a-$b))->sum<=1e-7;
}

#Check haydock coefficients for simple 1D system
my $B=zeroes(11)->xvals<5; #1D system
my $g=Photonic::Geometry::FromB->new(B=>$B, Direction0=>pdl([1])); #long
my $a=Photonic::LE::NR2::AllH->new(geometry=>$g, nh=>10);
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

#View 1D system as 2D. Transverse direction
my $Bt=zeroes(1,11)->yvals<5; #2D flat system
my $gt=Photonic::Geometry::FromB->new(B=>$Bt, Direction0=>pdl([1,0])); #trans
my $at=Photonic::LE::NR2::AllH->new(geometry=>$gt, nh=>10);
$at->run;
my $ast=$a->as;
my $bst=$a->bs;
my $b2st=$a->b2s;
ok(agree(pdl($at->iteration), 1), "Number of iterations 1D trans");
ok(agree(pdl($b2st->[0]), 1), "1D T b_0^2");
ok(agree(pdl($ast->[0]), $g->f), "1D T a_0");
ok(agree(pdl($b2st), pdl($bs)**2), "1D T b2==b^2");

{
    #check reorthogonalize with square array
    my $Bs=zeroes(15,15)->rvals<5;
    my $gs=Photonic::Geometry::FromB->new(B=>$Bs, Direction0=>pdl([1,0]));
    my $als=Photonic::LE::NR2::AllH
	->new(geometry=>$gs, nh=>2*15*15, reorthogonalize=>1,
	      accuracy=>machine_epsilon(), noise=>machine_epsilon(),
	      normOp=>1);  
    $als->run;
    ok($als->iteration <= 15*15, "No more iterations than dimensions");
    diag("Actual iterations: " . $als->iteration 
	 . " Actual orthogonalizations: ", $als->orthogonalizations);
}
{
    #check reorthogonalize with square array
    my $Bs=zeroes(15,15)->rvals<5;
    my $gs=Photonic::Geometry::FromB->new(B=>$Bs, Direction0=>pdl([1,0]));
    my $als=Photonic::LE::NR2::AllH
	->new(geometry=>$gs, nh=>2*15*15, reorthogonalize=>1,
	      accuracy=>machine_epsilon(), noise=>machine_epsilon(),
	      normOp=>1, stateFN=>"scratch/rem.dat");  
    $als->run;
    ok($als->iteration <= 15*15, "No more iterations than dimensions");
    diag("Actual iterations: " . $als->iteration 
	 . " Actual orthogonalizations: ", $als->orthogonalizations);
}
