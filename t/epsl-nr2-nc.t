#Test non cartesian primitive vectors.

use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromB;
use Photonic::LE::NR2::AllH;
use Photonic::LE::NR2::EpsL;

use Machine::Epsilon;
use List::Util;

use Test::More tests => 6;

#my $pi=4*atan2(1,1);

sub Cagree {
    my $a=shift;
    my $b=shift//0;
    my $prec=shift//1e-7;
    return (($a-$b)->Cabs2)->sum<=$prec;
}

my $ea=1+2*i;
my $eb=3+4*i;
#Check haydock coefficients for simple 1D system
my $B=zeroes(1,11)->yvals<5; #1D system
my $gl=Photonic::Geometry::FromB->new(B=>$B,
   primitive=>pdl([1,1],[1,-1]), Direction0=>pdl([1,-1])); #long
my $al=Photonic::LE::NR2::AllH->new(geometry=>$gl, nh=>10);
my $elo=Photonic::LE::NR2::EpsL->new(nr=>$al, nh=>10);
my $elv=$elo->evaluate($ea, $eb);
my $elx=1/((1-$gl->f)/$ea+$gl->f/$eb);
ok(Cagree($elv, $elx), "1D long epsilon");
#diag($elv);
#diag($elx);
is($elo->converged,1, "Converged");

#View 1D system as 2D. Transverse direction
my $Bt=zeroes(1,11)->yvals<5; #2D flat system
my $gt=Photonic::Geometry::FromB->new(B=>$Bt,
   primitive=>pdl([1,1],[1,-1]), Direction0=>pdl([1,1])); #trans
my $at=Photonic::LE::NR2::AllH->new(geometry=>$gt, nh=>10);
my $eto=Photonic::LE::NR2::EpsL->new(nr=>$at, nh=>10);
my $etv=$eto->evaluate($ea, $eb);
my $etx=(1-$gt->f)*$ea+$gt->f*$eb;
ok(Cagree($etv, $etx), "1D trans epsilon");
#diag($etv);
#diag($etx);
is($eto->converged,1, "Converged");

#Test chess nonorthogonal basis.
my $N=10;
my $Bc=zeroes(2*$N,2*$N);
#use (1,0) and (1,1) as primitive vectors.
$Bc=((((($Bc->xvals+$Bc->yvals) % (2*$N)) <$N) & ($Bc->yvals<$N))
     | (((($Bc->xvals+$Bc->yvals) % (2*$N)) >= $N) & ($Bc->yvals>=$N)));
my $gc=Photonic::Geometry::FromB->new(
    B=>$Bc, L=>pdl(1, sqrt(2)), primitive=>pdl([1,0],[1,1]),
    Direction0=>pdl([1,0]));
my $ac=Photonic::LE::NR2::AllH->new(geometry=>$gc, nh=>2000,
				    reorthogonalize=>1);
my $eco=Photonic::LE::NR2::EpsL->new(nr=>$ac, nh=>10000);
my $ecv=$eco->evaluate($ea, $eb);
my $ecx=sqrt($ea*$eb);
#diag($ecv);
#diag($ecx);
#diag($Bc);
ok(Cagree($ecv, $ecx, 1e-4), "Chess board");
#diag("O: ". $ac->orthogonalizations. " I: ".$ac->iteration);
is($eco->converged,1, "Converged");
