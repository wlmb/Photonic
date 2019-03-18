=head1 NAME

Photonic::LE::NR2::AllH

=head1 VERSION

version 0.011

=head1 SYNOPSIS

   use Photonic::LE::NR2::AllH;
   my $iter=Photonic::LE::NR2::AllH->new(geometry=>$geometry,nh=>$Nh); 
   $iter->run;
   my $haydock_as=$iter->as;
   my $haydock_bs=$iter->bs;
   my $haydock_b2s=$iter->b2s;

=head1 DESCRIPTION

Implements the calculation of Haydock coefficients and saves them for
later retrieval for the longitudinal dielectric function non retarded
calculation in a binary metamaterial.

=head1 METHODS

=over 4

=item * new(geometry=>$g, nh=>$nh) 

Initializes an Ph::LE::NR2::AllH object. $nh is the maximum number of desired
coefficients. All other arguments are as in Photonic::LE::NR2::OneH.

=item * run

Runs the iteration to completion

=item * All the Photonic::LE::NR2::OneH methods

=item * All the Photonic::Roles::AllH methods
=back

=head1 ACCESORS (read only)

=over 4

=item * nh

Maximum number of desired Haydock 'a' coefficients and states. The
number of b coefficients is one less.

=item * as

Array of Haydock a coefficients

=item * bs

Array of Haydock b coefficients

=item * b2s

Array of Haydock b coefficients squared

=back

=cut

package Photonic::LE::NR2::AllH;
$Photonic::LE::NR2::AllH::VERSION = '0.011';
use namespace::autoclean;
use Machine::Epsilon;
use PDL::Lite;
use PDL::NiceSlice;
#use Photonic::Utils qw(HProd);
use Moose;
use MooseX::StrictConstructor;
extends 'Photonic::LE::NR2::OneH';
with 'Photonic::Roles::AllH', 'Photonic::Roles::ReorthogonalizeR';

__PACKAGE__->meta->make_immutable;
    
1;
