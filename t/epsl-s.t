use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromEpsilon;
use Photonic::LE::S::AllH;
use Photonic::LE::S::EpsL;
use Photonic::Utils qw(HProd);

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
my $f=$B->sumover/$B->nelem;
my $epsilon=$ea*(1-$B)+$eb*$B;
my $gl=Photonic::Geometry::FromEpsilon->new(epsilon=>$epsilon,
					    Direction0=>pdl([1])); #long  
my $al=Photonic::LE::S::AllH->new(geometry=>$gl, nh=>10, epsilon=>$epsilon);
my $elo=Photonic::LE::S::EpsL->new(nr=>$al, nh=>10);
my $elv=$elo->epsL;
my $elx=1/((1-$f)/$ea+$f/$eb);
ok(Cagree($elv, $elx), "1D long epsilon");
is($elo->converged,1, "Converged");

#View 1D system as 2D. Transverse direction
my $Bt=zeroes(1,11)->yvals<5; #2D flat system
my $epsilont=$ea*(1-$Bt)+$eb*$Bt;
my $gt=Photonic::Geometry::FromEpsilon->new(epsilon=>$epsilont,
					    Direction0=>pdl([1,0])); #trans  
my $at=Photonic::LE::S::AllH->new(geometry=>$gt, nh=>10);
my $eto=Photonic::LE::S::EpsL->new(nr=>$at, nh=>10);
my $etv=$eto->epsL;
my $etx=(1-$f)*$ea+$f*$eb;
ok(Cagree($etv, $etx), "1D trans epsilon");
is($eto->converged,1, "Converged");
