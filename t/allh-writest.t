use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromB;
use Photonic::LE::NR2::AllH;

use Test::More; # tests => 11;

#Check haydock coefficients for simple 1D system
my $B=zeroes(11)->xvals<5; #1D system
my $g=Photonic::Geometry::FromB->new(B=>$B, Direction0=>pdl([1])); #long
my $a=Photonic::LE::NR2::AllH->new(geometry=>$g, nh=>10, stateFN=>"rem.txt");
$a->run;
is(1,1,"No test yet");
done_testing;