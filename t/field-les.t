=head1 COPYRIGHT NOTICE

Photonic - A perl package for calculations on photonics and
metamaterials.

Copyright (C) 2016 by W. Luis Mochán

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 1, or (at your option)
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA  02110-1301 USA

    mochan@fis.unam.mx

    Instituto de Ciencias Físicas, UNAM
    Apartado Postal 48-3
    62251 Cuernavaca, Morelos
    México

=cut

use strict;
use warnings;
use PDL;
use PDL::Constants qw(PI);
use Photonic::LE::S::Haydock;
use Photonic::LE::S::Field;
use Photonic::Geometry::FromB;
use Photonic::Geometry::FromEpsilon;
use Test::More tests => 6;
use lib 't/lib';
use TestUtils;

my $ea=1+2*i;
my $eb=3+4*i;

#Check field for simple 1D system
my $B=zeroes(11)->xvals<5; #1D system
my $epsilon=$ea*(1-$B)+$eb*$B;                                  # microscopic 1D diel. func.
my $gl=Photonic::Geometry::FromB->new(B=>$B, Direction0=>pdl([1])); #longitudinal
my $haydock=Photonic::LE::S::Haydock->new(geometry=>$gl, nh=>10,
   keepStates=>1, epsilon=>$epsilon);
my $flo=Photonic::LE::S::Field->new(haydock=>$haydock, nh=>10); # field longitudinal object
my $flv=$flo->field;                                            # field longitudinal value
my $fle=$flo->epsL;                                             # long. diel. func. from field
my $fla=1/$ea;                                                  # long. field. unnormalized
my $flb=1/$eb;
my $favl=$fla*(1-$gl->f)+$flb*($gl->f);                         # inv. long. response exact.
my $flex=1/$favl;                                               # exact long. response
($fla, $flb)=map {$_/$favl} ($fla, $flb);                       # normalized long. fields.
my $flx=($fla*(1-$B)+$flb*$B)->slice("*1");
ok(Cagree($flv, $flx), "1D long field") or diag "got: $flv\nexpected: $flx";
ok(Cagree($fle, $flex), "1D long response") or diag "got: $fle\nexpected: $flex";

#2D View from 1D superlattice.
my $Bt=zeroes(1,11)->yvals<5; #2D flat system
my $epsilont=$ea*(1-$Bt)+$eb*$Bt;                               # microscopic 2D diel. func.
my $gt=Photonic::Geometry::FromB->new(B=>$Bt, Direction0=>pdl([1,0]));
my $nt=Photonic::LE::S::Haydock->new(geometry=>$gt, nh=>10,
				  keepStates=>1, epsilon=>$epsilont);
my $fto=Photonic::LE::S::Field->new(haydock=>$nt, nh=>10);      # field transverse object
my $ftv=$fto->field;                                            # field transverse value
my $fte=$fto->epsL;                                             # trans. diel. func. from field.
my $ftx=pdl([1, 0])->r2C;                                       # exact trans. field.
ok(Cagree($ftv, $ftx), "1D trans field") or diag "got: $ftv\nexpected: $ftx";;
my $favt=$ea*(1-$gt->f)+$eb*($gt->f);                         # exact trans. response
ok(Cagree($fte, $favt), "1D trans response") or diag "got: $fte\nexpected: $favt";

# check raw fields
my $flv_raw=$flo->rawfield;
my $flx_raw=-4*PI*((1-$B)/$ea+$B/$eb)->dummy(0);
ok(Cagree($flv_raw, $flx_raw), "1D long raw field");

my $ftv_raw=$fto->rawfield;
my $ftx_raw=-(4*PI+0*i)/($ea*(1-$gt->f)+$eb*$gt->f)*pdl([1,0]);
ok(Cagree($ftv_raw, $ftx_raw), "1D transverse raw field");
