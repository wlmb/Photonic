use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromB;
use Photonic::LE::NR2::AllH;
use Photonic::LE::NR2::EpsTensor;

use Machine::Epsilon;
use List::Util;

use Test::More tests => 4;

#my $pi=4*atan2(1,1);

sub Cagree {    
    my $a=shift;
    my $b=shift//0;
    return (($a-$b)->Cabs2)->sum<=1e-7;
}

my $ea=1+2*i;
my $eb=3+4*i;
#Check haydock coefficients for simple 1D system
my $B=zeroes(11)->xvals<5; #1D system
my $gl=Photonic::Geometry::FromB->new(B=>$B); #long
my $elo=Photonic::LE::NR2::EpsTensor->new(geometry=>$gl, nh=>10);
my $elv=$elo->evaluate($ea, $eb);
my $elx=1/((1-$gl->f)/$ea+$gl->f/$eb);
ok(Cagree($elv, $elx), "1D long epsilon");
is($elo->converged,1, "Converged");

#View 2D from 1D superlattice.
my $Bt=zeroes(1,11)->yvals<5; #2D flat system
my $gt=Photonic::Geometry::FromB->new(B=>$Bt); #trans
my $eto=Photonic::LE::NR2::EpsTensor->new(geometry=>$gt, nh=>10);
my $etv=$eto->evaluate($ea, $eb);
my $etx=(1-$gt->f)*$ea+$gt->f*$eb;
my $etenx=pdl([$etx, 0+0*i],[0+0*i, $elx])->complex;
ok(Cagree($etv, $etenx), "1D trans epsilon");
is($eto->converged,1, "Converged");
