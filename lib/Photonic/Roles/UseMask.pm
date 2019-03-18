package Photonic::Roles::UseMask;
$Photonic::Roles::UseMask::VERSION = '0.011';
use namespace::autoclean;
use PDL::Lite;
use PDL::NiceSlice;
use PDL::Complex;
use Moose::Role;

has 'use_mask'=>(is=>'ro', default=>1, documentation=>'Use mask if present');
has 'mask'=>(is=>'ro', lazy=>1, builder=>'_build_mask',
    documentation=>'Mask in reciprocal space');
requires qw(ndims);

sub _build_mask { #default mask kills G_max for even dims.
    my $self=shift;
    my $ndims=$self->ndims;
    my $dims=$self->dims;
    my $mask=PDL->ones(@$dims);
    my $masked=0; 
    foreach(0..$ndims-1){
	my $N=$dims->[$_]; 
	next unless $N%2==0; #ignore odd dimensions.
	$mask->mv($_,0)->(($N/2)).=0; #zero terms corresponding to \pm G_max
	$masked=1;
    }
    return $mask if $masked;
    return undef;
}

no Moose::Role;

1;
