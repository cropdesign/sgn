# perl bin/pedigreeplausibility.pl -H localhost -D fixture -p 1 -o sink.txt -f bin/ped.txt
# psql -U postgres
# \c fixture
# select * from stock join cvterm on (type_id = cvterm_id) where cvterm.name = 'accession';
# ssh production@cassava-devel.sgn.cornell.edu pw #cp_/root/*_?
# perl bin/pedigreeplausibility.pl -H db5.sgn.cornell.edu -D cxgn_cassava -p 3 -o sink.txt -f bin/ped.txt
# postgres pw Eo0vair1

use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;
use Bio::Chado::Schema;
use CXGN::DB::InsertDBH;
use CXGN::Chado::Stock;
use CXGN::Genotype;
use CXGN::Genotype::Search;

our ($opt_H, $opt_D, $opt_p, $opt_o, $opt_f); # host, database, genotyping protocol_id, out, in
getopts('H:D:p:o:f:');

if (!$opt_p) {
    print STDERR "Need -p with genotyping protocol id.\n";
    exit();
}

my $protocol_id = $opt_p;

my $filename = $opt_f;

my %pedigreehash;
my $childname;
my @pedigreearray;

open(IN, $filename) or die "Could not open file $filename $!";

while (my $row = <IN>){
@pedigreearray = split(/\s/ ,$row);
$childname = $pedigreearray[0];
$pedigreehash{$childname} = "1";
}

my $dbh = CXGN::DB::InsertDBH->new( {
    dbhost => $opt_H,
    dbname => $opt_D,
    dbuser => "postgres",
				    }
    );

my $OUT;
my $is_stdin =0;
}

my $schema = Bio::Chado::Schema->connect(sub { $dbh });

my $accession_cvterm_id = $schema->resultset("Cv::Cvterm")->find({ name=> "accession" })->cvterm_id();

my $stock_rs = $schema->resultset("Stock::Stock")->search( { type_id => $accession_cvterm_id });

my @scores;

while (my $row = $stock_rs->next()) {
    print STDERR "working on accession ".$row->uniquename()."\n";
    unless($pedigreehash{$row->uniquename()}){next;}
    my $stock = CXGN::Chado::Stock->new($schema, $row->stock_id());
    my @parents = $stock->get_direct_parents();

    if (@parents == 2) {

	  my $gts = CXGN::Genotype::Search->new( {
	    bcs_schema => $schema,
	    accession_list => [ $row->stock_id ],
	    protocol_id => $protocol_id,
							    });

  my @self_gts = $gts->get_genotype_info_as_genotype_objects();
	$gts = CXGN::Genotype::Search->new( {
	    bcs_schema => $schema,
	    accession_list => [ $parents[0]->[0]],
	    protocol_id => $protocol_id,
							    });

	my @mom_gts = $gts->get_genotype_info_as_genotype_objects();

	$gts = CXGN::Genotype::Search->new( {
	    bcs_schema => $schema,
	    accession_list => [ $parents[1]->[0]],
	    protocol_id => $protocol_id,
							    });

	my (@dad_gts) = $gts->get_genotype_info_as_genotype_objects();

	if (! (@self_gts)) {
	    print STDERR "Genotype of accession ".$row->uniquename()." not available. Skipping...\n";
	    next;
	}
	if (!@mom_gts) {
	    print STDERR "Genotype of female parent missing. Skipping.\n";
	    next;
	}
	if (! @dad_gts) {
	    print STDERR "Genotype of male parent missing. Skipping.\n";
	    next;
	}
  if ($opt_o) {
      open($OUT, '>>', $opt_o);


  #else {
      #$OUT =  *STDIN;
      #$is_stdin = 1;
  print $OUT "at mom genos".@mom_gts;
  print $OUT "at dad genos".@dad_gts;
  print $OUT "at child genos".@self_gts;
  #print $OUT "d child genos".$self_gts;
  #print $OUT "d father genos".$dad_gts;
  #print $OUT "d mother genos".$mom_gts;
}
##check length
##index of array for
	foreach my $s (@self_gts) {
	    foreach my $m (@mom_gts) {
		foreach my $d (@dad_gts) {
		    my ($concordant, $discordant, $non_informative) =
			$s->compare_parental_genotypes($m, $d);
		    my $score = $concordant / ($concordant + $discordant);
		    push @scores, $score;

		    print $OUT join "\t", map { ($_->name(), $_->id()) } ($s, $m, $d);
		    print $OUT "\t$score\n";
		}
	    }
	}
    }
}

$dbh->disconnect();

print STDERR "Done.\n";
