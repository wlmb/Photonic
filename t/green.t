use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromEpsilon;
use Photonic::WE::S::Metric;
use Photonic::WE::S::AllH;
use Photonic::WE::S::Green;

use Machine::Epsilon;
use List::Util;

use Test::More tests => 13;

sub agree {    
    my $a=shift;
    my $b=shift//0;
    return (($a-$b)*($a-$b))->sum<=1e-7;
}

sub Cagree {    
    my $a=shift;
    my $b=shift//0;
    return (($a-$b)->Cabs2)->sum<=1e-7;
}

#Check green for simple 1D system
my ($ea, $eb)=(r2C(1), r2C(2));
my $f=6/11;
my $eps=$ea*(zeroes(11,1)->xvals<5)+ $eb*(zeroes(11)->xvals>=5)+0*i; 
my $g=Photonic::Geometry::FromEpsilon
    ->new(epsilon=>$eps); 
my $m=Photonic::WE::S::Metric->new(
    geometry=>$g, epsilon=>pdl(1), wavenumber=>pdl(2e-5),
    wavevector=>pdl([1,0])*1e-8);  
my $gr=Photonic::WE::S::Green->new(nh=>10, metric=>$m);
my $grv=$gr->greenTensor;
print $grv;