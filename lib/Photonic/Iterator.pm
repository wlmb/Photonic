package Photonic::Iterator;
$Photonic::Iterator::VERSION='0.011';
use base Exporter;
@EXPORT_OK=qw(iterator nextval);
%EXPORT_TAGS=(all=>\@EXPORT_OK);

sub nextval ($) { $_[0]->(); }
sub iterator(&) { return $_[0] }
