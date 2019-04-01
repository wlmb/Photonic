=head1 NAME

Photonic::Iterator

=head1 VERSION

version 0.011

=head1 SYNOPSYS

     use Photonic::Iterator qw(iterator nextval);
     sub makeiterator {
         my $a=shift;
         my $n=0;
         my $it = iterator {
              return $a->[$n++];
         };
         return $it;
     }
     my @a=(1..100);
     my $it=makeiterator(\@a);
     while(my $x=nextval($it)){
        #do something with $x
     }

=head1 DESCRIPTION

Facilitate the creation of iterators. Iterators are objects that
produce data items, one at a time, until depleted. They are
implemented as functions (closures) that have state (implemented as
lexical variables) to know which is the next item. The example above
creates an iterator for walking over the elements of an array
reference. 

=head1 EXPORTABLE FUNCTIONS 

=over 4

=item iterator

    my $i=iterator {...}

creates an iterator from the function {...}

=item nextval

   my $x=nextval($i)

returns the next value from iterator $i.


=back

=cut

package Photonic::Iterator;
$Photonic::Iterator::VERSION='0.011';
use base Exporter;
@EXPORT_OK=qw(iterator nextval);
%EXPORT_TAGS=(all=>\@EXPORT_OK);

sub nextval ($) { $_[0]->(); }
sub iterator(&) { return $_[0] }
