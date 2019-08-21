package Photonic::Iterator;
$Photonic::Iterator::VERSION='0.011';
use base Exporter;
@EXPORT_OK=qw(iterator nextval);
%EXPORT_TAGS=(all=>\@EXPORT_OK);

sub nextval ($) { $_[0]->(); }
sub iterator(&) { return $_[0] }


=head1 NAME

Photonic::Iterator

=head1 VERSION

version 0.011

=head1 SYNOPSIS

    use Photonic::Iterator qw(iterator nextval);
    
