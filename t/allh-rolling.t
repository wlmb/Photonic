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

# Haydock calculators without keepStates hold only 3 states in a
# rotating buffer. Check they give exactly the same results as keeping
# all states, and check storing/loading with state files.

use strict;
use warnings;
use PDL;
# No PDL::NiceSlice: it would turn the code-ref calls ->(...) into slices
use PDL::IO::FastRaw;
use List::Util ();
use Photonic::Geometry::FromB;
use Photonic::Geometry::FromEpsilon;
use Photonic::Geometry::FromEpsilonTensor;
use Photonic::WE::R2::Metric;
use Photonic::WE::R2::Green;
use Photonic::WE::S::Metric;
use Photonic::WE::S::Green;
use Photonic::WEM::S::Metric;
use Photonic::WEM::S::Green;
use Photonic::WEM::ST::Metric;
use Photonic::WEM::ST::Green;

use Test::More;
use lib 't/lib';
use TestUtils;

my @coeffs = qw(as bs b2s cs bcs gs);
# Same arithmetic in the same order: results must agree exactly
sub same { my ($x, $y) = @_; ($x-$y)->abs->max == 0 }
sub slots { $_[0]->_state_pdl->dim(-1) }

# 2D system that needs all 10 iterations, so the 3-slot buffer rotates
my $B2 = zeroes(11,11)->rvals<3;
my $nh = 10;

{ # 1. LE::NR2: same coefficients, 3 slots, states() refused
    my ($k1) = make_store($B2, [1,0], $nh, {keepStates=>1});
    my ($k0) = make_store($B2, [1,0], $nh, {});
    ok($k1->iteration > 3, "LE::NR2 iterates beyond the 3-slot buffer");
    is($k0->iteration, $k1->iteration, "LE::NR2 same iterations");
    ok(same($k0->$_, $k1->$_), "LE::NR2 same $_ without keepStates") for @coeffs;
    is(slots($k0), 3, "LE::NR2 3 slots without keepStates");
    is(slots($k1), $nh+1, "LE::NR2 nh+1 slots with keepStates");
    eval { $k0->states };
    like($@, qr/states were not kept/, "LE::NR2 states() refused without keepStates");
    my ($sf) = make_store($B2, [1,0], $nh, {stateFN=>make_fn()});
    ok($sf->keepStates, "stateFN implies keepStates");
}

# 2. Green's tensors, built with keepStates 0 and 1
my ($ea, $eb) = (r2C(1), r2C(2));
my $eps2 = $ea + ($eb-$ea)*$B2;
my %make_green = (
    'WE::R2' => sub {
	my $g=Photonic::Geometry::FromB->new(B=>$B2);
	my $m=Photonic::WE::R2::Metric->new(geometry=>$g, epsilon=>pdl(1),
	    wavenumber=>pdl(2e-5), wavevector=>pdl([1,0])*1e-8);
	Photonic::WE::R2::Green->new(nh=>$nh, metric=>$m, epsB=>$eb, @_);
    },
    'WE::S' => sub {
	my $g=Photonic::Geometry::FromEpsilon->new(epsilon=>$eps2);
	my $m=Photonic::WE::S::Metric->new(geometry=>$g, epsilon=>pdl(1),
	    wavenumber=>pdl(2e-5), wavevector=>pdl([1,0])*1e-8);
	Photonic::WE::S::Green->new(nh=>$nh, metric=>$m, @_);
    },
    'WEM::S' => sub {
	my $B3 = zeroes(7,7,7)->rvals<2;
	my $g=Photonic::Geometry::FromEpsilon->new(epsilon=>$ea+($eb-$ea)*$B3);
	my $m=Photonic::WEM::S::Metric->new(geometry=>$g, epsilon=>pdl(1),
	    mu=>ones(7,7,7)->r2C,
	    wavenumber=>pdl(0.01), wavevector=>pdl([0.01,0,0]));
	Photonic::WEM::S::Green->new(nh=>$nh, metric=>$m, @_);
    },
    'WEM::ST' => sub {
	my $B3 = (zeroes(7,7,7)->rvals<2)->dummy(0)->dummy(0);
	my $id = identity(3)->slice(":,:,*1,*1,*1")->r2C;
	my $g=Photonic::Geometry::FromEpsilonTensor->new(
	    epsilon=>$ea*$id+($eb-$ea)*$B3*$id);
	my $m=Photonic::WEM::ST::Metric->new(geometry=>$g, epsilon=>pdl(1),
	    mu=>$id*ones(1,1,7,7,7),
	    wavenumber=>pdl(0.01), wavevector=>pdl([0.01,0,0]));
	Photonic::WEM::ST::Green->new(nh=>$nh, metric=>$m, @_);
    },
);
for my $name (sort keys %make_green) {
    my $g1 = $make_green{$name}->(keepStates=>1);
    my $g0 = $make_green{$name}->();
    ok(same($g0->greenTensor, $g1->greenTensor),
       "$name same greenTensor without keepStates");
    my @h0 = @{$g0->haydock};
    my @h1 = @{$g1->haydock};
    ok(List::Util::max(map $_->iteration, @h0) > 3,
       "$name iterates beyond the 3-slot buffer");
    ok((List::Util::all { slots($_) == 3 } @h0),
       "$name Haydock calculators hold 3 states without keepStates");
    ok((List::Util::all { slots($_) == $nh+1 } @h1),
       "$name Haydock calculators hold nh+1 states with keepStates");
    for my $c (@coeffs) {
	ok((List::Util::all { same($h0[$_]->$c, $h1[$_]->$c) } 0..$#h0),
	   "$name same $c without keepStates");
    }
}

# 3. WEM Greens pass smallH to their Haydock calculators
for my $name (qw(WEM::S WEM::ST)) {
    my $gr = $make_green{$name}->(smallH=>0.123);
    is($gr->haydock->[0]->smallH, 0.123, "$name passes smallH");
}

# Storing and loading. Reference run keeping all states:
my ($full) = make_store($B2, [1,0], $nh, {keepStates=>1});
{ # 4. load all states into an object built without keepStates
    my $fn = make_fn();
    make_store($B2, [1,0], $nh, {keepStates=>1, storeAllFN=>$fn});
    my ($l) = make_store($B2, [1,0], $nh, {loadAllFN=>$fn});
    ok($l->keepStates, "keepStates restored from file");
    is(slots($l), $nh+1, "3-slot buffer grown to nh+1 on load");
    ok(same($l->states, $full->states), "states restored into grown buffer");
}
{ # 5. stored without states: coefficients load, resuming refused
    my $fn = make_fn();
    make_store($B2, [1,0], 2, {storeAllFN=>$fn});
    my ($c) = eval { make_store($B2, [1,0], 2, {loadAllFN=>$fn}) };
    ok(!$@ && same($c->as, $full->as->slice('0:1')),
       "coefficients load without states for the same nh") or diag $@;
    eval { make_store($B2, [1,0], $nh, {loadAllFN=>$fn}) };
    like($@, qr/Cannot resume Haydock iteration/,
	 "resuming to a larger nh without states refused");
}
{ # 6. stale state file of a different system is recreated, not reused
    my $fn = make_fn();
    make_store(zeroes(11)->xvals<5, [1], $nh, {stateFN=>$fn});
    my ($n) = eval { make_store($B2, [1,0], $nh, {stateFN=>$fn}) };
    ok(!$@, "stale state file of another system does not stop the run")
	or diag $@;
    ok($n && same($n->states, $full->states), "states correct in recreated file");
}
{ # 7. resume into the same file with a larger nh: the file grows
    my $fn = make_fn();
    make_store($B2, [1,0], 4, {keepStates=>1, storeAllFN=>$fn});
    my ($r) = make_store($B2, [1,0], $nh,
	{keepStates=>1, loadAllFN=>$fn, storeAllFN=>$fn});
    is($r->iteration, $full->iteration, "resumed into same file: iterations");
    ok(same($r->$_, $full->$_), "resumed into same file: $_") for @coeffs;
    ok(same($r->states, $full->states), "resumed into same file: states");
    is(mapfraw("$fn.states", {ReadOnly=>1})->dim(-1), $nh+1,
       "state file grown to nh+1 slots");
    ok(!-e "$fn.states.old", "temporary .old file removed");
}
{ # 8. resume into the same file whose states file is missing
    my $fn = make_fn();
    make_store($B2, [1,0], 4, {keepStates=>1, storeAllFN=>$fn});
    unlink "$fn.states", "$fn.states.hdr";
    eval { make_store($B2, [1,0], $nh,
	{keepStates=>1, loadAllFN=>$fn, storeAllFN=>$fn}) };
    like($@, qr/to resume from does not exist/,
	 "missing state file refused when resuming into it");
}
{ # 9. resume into the same file holding states of a different system
    my $fn = make_fn();
    make_store(zeroes(11)->xvals<5, [1], 4, {keepStates=>1, storeAllFN=>$fn});
    eval { make_store($B2, [1,0], $nh,
	{keepStates=>1, loadAllFN=>$fn, storeAllFN=>$fn}) };
    like($@, qr/belongs to a different system/,
	 "state file of a different system refused when resuming into it");
}

done_testing;
