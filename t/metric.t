use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromB;
use Photonic::Metric::R2;
use Test::More tests => 7;

#my $pi=4*atan2(1,1);

sub agree {    
    my $a=shift;
    my $b=shift//0;
    return (($a-$b)*($a-$b))->sum<=1e-7;
}

my $B=zeroes(11)->xvals<5; #1D system
my $g=Photonic::Geometry::FromB->new(B=>$B);
my $gGG=Photonic::Metric::R2->new(geometry=>$g, epsilon=>pdl(2),
   wavenumber=>pdl(1), wavevector=>pdl([1])); 
my $v=$gGG->value;
ok($v->ndims==3,"Number of dimensions of metric for 1d");
ok(agree(pdl($v->dims),pdl(1,1,11)), "Actual dimensions of metric for 1d");
ok(agree($v, ones(1,1,11)), "Actual metric for 1d");

$B=zeroes(1,11)->xvals<5; #2D system
$g=Photonic::Geometry::FromB->new(B=>$B);
$gGG=Photonic::Metric::R2->new(geometry=>$g, epsilon=>pdl(2),
   wavenumber=>pdl(1), wavevector=>pdl([0,1])); 
$v=$gGG->value;
ok($v->ndims==4,"Number of dimensions of metric for 2d");
ok(agree(pdl($v->dims),pdl(2,2,1,11)), "Actual dimensions of metric for 2d");
ok(agree($v->((1),(1)), ones(1,11)),
			"Longitudinal component of metric in 2D");
ok(agree($v->((0),(0)), 2/(2-((pdl([0,1])+$g->G)**2)->sumover)),
			"Transverse component of metric in 2D");
