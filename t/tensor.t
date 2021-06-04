use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;
use Photonic::Geometry::FromB;
use Photonic::Utils qw(any_complex tensor);
use Test::More;

foreach my $D (1..5){
    my $B=zeroes((1)x$D); # D-dimentional characteristic function
    my $g=Photonic::Geometry::FromB->new(B=>$B);
    my $T=random($D,$D)+i()*random($D,$D); #random 2nd rank complex tensor
    $T=($T+$T->mv(1,2))/2; #symmetrize
    my $pairs=$g->unitPairs; # indices: xy,n
    my $projections=($pairs*($T*$pairs->(:,*1))->sumover)->sumover; #indices: ri,dyad
    my $lu=$g->unitDyadsLU;
    my $newT=tensor($projections, $lu, $D, 2); #reconstructed tensor
    ok all(approx($newT, $T)), "$D-D 2nd rank complex tensor" or diag "got: $newT, expected $T";
}

foreach my $D(1..5){
    my $B=zeroes((1)x$D); # $_-dimentional characteristic function
    my $g=Photonic::Geometry::FromB->new(B=>$B);
    my $T=random($D,$D,$D)+i()*random($D,$D,$D)+1+i; #random 3d rank complex tensor
    $T=($T+$T->mv(2,3))/2; #symmetrize on last two indices, leave first alone. ri,xy,xy,xy
    my $pairs=$g->unitPairs; # indices: xy,n
    my $projections=($pairs->(:,*1)*($T->mv(1,-1)*$pairs->(:,*1,*1))->sumover)
	->sumover->mv(1,-1); #indices: ri,dyad,xy
    my $lu=$g->unitDyadsLU;
    my $newT=tensor($projections, $lu, $D, 3); #reconstructed tensor
    ok all(approx($newT, $T)), "$D-D 3rd rank complex tensor" or diag "got: $newT, expected $T";
}

done_testing;
