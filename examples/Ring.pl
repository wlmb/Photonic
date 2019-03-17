#! /usr/bin/env perl
# Retarded vs non retarded MG.
use strict;
use warnings;
use feature qw(say);

#use DBI;

use Getopt::Long;
use List::Util;
use constant PI=>4*atan2(1,1);

use lib qw(lib);
use Photonic::Geometry::FromB;
use Photonic::LE::NR2::EpsL;
use Photonic::WE::R2::Metric;
use Photonic::WE::R2::AllH;
use Photonic::WE::R2::EpsilonP;
use PDL;
use PDL::NiceSlice;
use PDL::Complex;

set_autopthread_targ(4);;
set_autopthread_size(4);;


my $nh=30;
my $N=30;
my $r0=0;
my $r1= 0.45;
my $B=!ring($N,$r0,$r1);
my $l=1;
my $a=10; # Try with large value of $a to check non retarded results
my $epsA=pdl(4);
my ($hnu_all,$epsBall)=eps("au");
my $k=0.01;
my $q=2*PI/$a;
my $elem=$epsBall->dim(1);
my $filename="Au_cyl_r${r1}_k${k}_N${N}_Nh${nh}"; $filename=~s/\./_/g; $filename.=".dat";
open(OUT, ">", "data/$filename") or die "Couldn't open $filename for writing. $!";
print OUT "#   hnu       epsNR-L_re  epsNR-L_im   epsR-T_re    epsR-T_im \n";  

my $gmtnr=Photonic::Geometry::FromB->new(B=>$B,L=>pdl($l,$l),Direction0=>pdl(0,1));
my $gmtr=Photonic::Geometry::FromB->new(B=>$B,L=>pdl($l,$l));
my $allh=Photonic::LE::NR2::AllH->new(geometry=>$gmtnr, nh=>$nh);
my $nr=Photonic::LE::NR2::EpsL->new(nr=>$allh, nh=>$nh);
my @out=();
    
for(my $j=0;$j<$elem;$j++){
    my $epsB=$epsBall(,($j));
    my $hnu=$hnu_all(($j));
    
    my $enr=$nr->evaluate($epsA->r2C,$epsB);
    
    my $m=Photonic::WE::R2::Metric->new(geometry=>$gmtr, epsilon=>$epsA,
					  wavenumber=>pdl($q),
					  wavevector=>pdl([$k,0]));
    my $h=Photonic::WE::R2::AllH->new(metric=>$m,
					polarization=>pdl(0,1)->r2C, nh=>$nh);
    my $e=Photonic::WE::R2::EpsilonP->new(haydock=>$h, nh=>$nh);
    my $er=$e->evaluate($epsB);
    say OUT join " ", $hnu, $enr->re, $enr->im, $er->re, $er->im; 
    my $linea=pdl($hnu, $enr->re, $enr->im, $er->re, $er->im);
    push @out, $linea;
}

# If you have installed Gnuplot and Gnuplot pdl module, try uncomment lines before
#use PDL::Graphics::Gnuplot;
#my $o=pdl(@out)->mv(-1,0);
#print $o->info,"\n";
#my $x=$o(,0);
#my $y=$o->(:,1:-1);
#my $wgp=gpwin('qt');
#$wgp->plot(xr=>[1.2,3],with=>'lines',$x,$y);


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
	$epsilonIn="data/eps_au.d";
    }
    open(fileEpsIn,'<',"$epsilonIn") || die "can't open $epsilonIn:\n";
    (my $h_nu,my $eps_re,my $eps_im) = rcols *fileEpsIn;
    close(fileEpsIn);
    my $eps=$eps_re+i*$eps_im;
    return ($h_nu, $eps);
}

    

