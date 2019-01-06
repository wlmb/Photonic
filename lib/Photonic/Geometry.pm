package Photonic::Geometry;
$Photonic::Geometry::VERSION='0.010';
use Module::PluginFinder;

my $finder=Module::PluginFinder->new(
    search_path=>'Photonic::Geometry',
    filter=>\&filter
    );
sub filter {
    my ($class, $searchkey)=@_;
    $class->_forme($searchkey);
};

sub new {
    my $self=shift;
    my %args=@_;
    my $module=$finder->construct(\%args, %args);
}
    
1;
