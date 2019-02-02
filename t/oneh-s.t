use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromEpsilon;
use Photonic::LE::S::OneH;

use Test::More tests => 4;

#my $pi=4*atan2(1,1);

sub Cagree {    
    my $a=shift;
    my $b=shift//0;
    return (($a-$b)->Cabs2)->sum<=1e-7;
}

#Check haydock coefficients for simple 1D system
#1D system e=1 or 2
my ($ea, $eb)=(1+2*i, 3+4*i);
my $eps=$ea*(zeroes(11)->xvals<5)+ $eb*(zeroes(11)->xvals>=5)+0*i; 
my $g=Photonic::Geometry::FromEpsilon
    ->new(epsilon=>$eps, Direction0=>pdl([1])); 
my $o=Photonic::LE::S::OneH->new(geometry=>$g);
$o->iterate;
my $f=6/11;
ok(Cagree(pdl($o->current_a), $ea*(1-$f)+$eb*$f), "1D a_0");
ok(Cagree(pdl($o->next_b2), ($eb-$ea)**2*$f*(1-$f)), "1D b_1^2");
$o->iterate;
ok(Cagree(pdl($o->current_a), $ea*$f+$eb*(1-$f)), "1D a_1");
ok(Cagree(pdl($o->next_b2), 0), "1D b_2^2");
