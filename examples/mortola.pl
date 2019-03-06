#! /usr/bin/env perl
# Proof of NP module using Mortola results
use strict;
use warnings;
use feature qw(say);

#use DBI;

use Getopt::Long;
use List::Util;

use lib qw(../../Photonic/lib);
use Photonic::Geometry::FromEpsilon;
use Photonic::LE::S::AllH;
use Photonic::LE::S::EpsL;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;

set_autopthread_targ(4);;
set_autopthread_size(4);;


my $nh=10;
my $N=50;
my $l=1;
my $epsA=r2C 2+3*i;
my $epsB=r2C 3+4*i;
my $epsC=r2C 2+3*i;
my $epsD=r2C 3+4*i;
my $e=FourPhases($N,$epsA,$epsB,$epsC,$epsD);

my $filename="epsM_S_eA${epsA}_eB${epsB}_N${N}_Nh${nh}"; $filename=~s/\./_/g; $filename.=".dat";
open(OUT, ">", "../data/$filename") or die "Couldn't open $filename for writing. $!";
print OUT "#   hnu       epsNR-L_re  epsNR-L_im   epsR-T_re    epsR-T_im \n";  

my $gmtnr=Photonic::Geometry::FromEpsilon->new(epsilon=>$e,L=>pdl($l,$l),Direction0=>pdl(1,0));
my $allh=Photonic::LE::S::AllH->new(geometry=>$gmtnr, nh=>$nh);
my $nr=Photonic::LE::S::EpsL->new(nr=>$allh, nh=>$nh);
    
my $enr=$nr->epsL;
    
say OUT join " ", $enr->re, $enr->im; 



sub checkerboard {
    my $N=shift;
    my $B=zeroes(2*$N+1,2*$N+1);
    my $c=$B->ndcoords-pdl($N,$N);
    my $i=whereND($B, $c->((0))->abs > $c->((1))->abs);
    $i.=1;
    #$i=wheLreND($B, $c->((0))->abs == $c->((1))->abs);
    #$i.=0.5; #esto vuelve isotrópico al sistema, pero ya no es binario
    return $B;
}

sub ajedrez { #checkboard
    my $N=shift;
    my $z=zeroes(2*$N+1,2*$N+1); #change to admit arbitrary lattice
    my $r=$z->ndcoords-pdl($N,$N);
    my $signopos=$r->((1))*$r->((0))>0;
#    my $t=$q<(.75+0*cos(3*$theta))*$N;
    return $signopos;
}

sub FourPhases { #checkboard
    my $N=shift;
    my $eA=shift;
    my $eB=shift;
    my $eC=shift;
    my $eD=shift;
#    my $z=zeroes(2*$N+1,2*$N+1); #change to admit arbitrary lattice
    $eB=($eB*ones($N,$N))->mv(0,-1);
    my $z=$eB->glue(1,($eA*ones($N,$N+1))->mv(0,-1));
    $eC=($eC*ones($N+1,$N))->mv(0,-1);
    $z=$z->glue(0,($eC->glue(1,($eD*ones($N+1,$N+1))->mv(0,-1))));
    return $z->mv(-1,0);
}


sub eps{
    my $epsi=shift;
    my $epsilonIn;
    if($epsi eq "au"){
	$epsilonIn="/home/gortiz/Metas/data/eps_au.d";
    }
    open(fileEpsIn,'<',"$epsilonIn") || die "can't open $epsilonIn:\n";
    (my $h_nu,my $eps_re,my $eps_im) = rcols *fileEpsIn;
    close(fileEpsIn);
    my $eps=$eps_re+i*$eps_im;
    return ($h_nu, $eps);
}

    

