#! /usr/bin/env perl
# Retarded vs non retarded MG.
use strict;
use warnings;
use feature qw(say);

#use DBI;

use Getopt::Long;
use List::Util;
use constant PI=>4*atan2(1,1);

use lib qw(../../Photonic/lib);
use Photonic::Geometry::FromB;
use Photonic::NonRetarded::EpsL;
use Photonic::Retarded::Metric;
use Photonic::Retarded::AllH;
use Photonic::Retarded::EpsilonP;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;

set_autopthread_targ(4);;
set_autopthread_size(4);;


my $nh=10;
my $N=50;
my $r0=0;
my $r1= 0.45;
my $B=!ring($N,$r0,$r1);
my $l=1;
my $epsA=pdl(4);
my ($hnu_all,$epsBall)=eps("au");
my $k=0.01;
my $q=2*PI/5;
my $elem=$epsBall->dim(1);
my $filename="Au_cyl_r${r1}_k${k}_N${N}_Nh${nh}"; $filename=~s/\./_/g; $filename.=".dat";
open(OUT, ">", "../data/$filename") or die "Couldn't open $filename for writing. $!";
print OUT "#   hnu       epsNR-L_re  epsNR-L_im   epsR-T_re    epsR-T_im \n";  

my $gmtnr=Photonic::Geometry::FromB->new(B=>$B,L=>pdl($l,$l),Direction0=>pdl(0,1));
my $gmtr=Photonic::Geometry->new(B=>$B,L=>pdl($l,$l));
my $allh=Photonic::NonRetarded::AllH->new(geometry=>$gmtnr, nh=>$nh);
my $nr=Photonic::NonRetarded::EpsL->new(nr=>$allh, nh=>$nh);
    
for(my $j=0;$j<$elem;$j++){
    my $epsB=$epsBall(,($j));
    my $hnu=$hnu_all(($j));
    
    my $enr=$nr->evaluate($epsA->r2C,$epsB);
    
    my $m=Photonic::Retarded::Metric->new(geometry=>$gmtr, epsilon=>$epsA,
					  wavenumber=>pdl($q),
					  wavevector=>pdl([$k,0]));
    my $h=Photonic::Retarded::AllH->new(metric=>$m,
					polarization=>pdl(0,1)->r2C, nh=>$nh);
    my $e=Photonic::Retarded::EpsilonP->new(haydock=>$h, nh=>$nh);
    my $er=$e->evaluate($epsB);
    say OUT join " ", $hnu, $enr->re, $enr->im, $er->re, $er->im; 
}





sub ring {
    my $N=shift;
    my $ra=(shift)*(2*$N+1);  #Note: it was $N
    my $rb=(shift)*(2*$N+1);
    my $z=zeroes(2*$N+1, 2*$N+1);
    my $B=($z->rvals<=$rb) & ($z->rvals >= $ra);
    return $B;
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

    

