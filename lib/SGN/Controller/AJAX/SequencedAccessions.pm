
package SGN::Controller::AJAX::SequencedAccessions;

use Moose;

use Data::Dumper;
use DateTime;
use List::MoreUtils qw(uniq);
use CXGN::Stock::SequencingInfo;

__PACKAGE__->config(
    default   => 'application/json',
    stash_key => 'rest',
    map       => { 'application/json' => 'JSON', 'text/html' => 'JSON' },
   );


BEGIN { extends 'Catalyst::Controller::REST' };

sub get_all_sequenced_stocks :Path('/ajax/genomes/sequenced_stocks') {
    my $self = shift;
    my $c = shift;

    my $user_id = $c->user() ? $c->user->get_object()->get_sp_person_id() : undef;
    my $schema = $c->dbic_schema("Bio::Chado::Schema");
    my @sequenced_stocks = CXGN::Stock::SequencingInfo->all_sequenced_stocks($schema);

    my @info = $self->retrieve_sequencing_infos($schema, $user_id, @sequenced_stocks);

    $c->stash->{rest} = { data => \@info };
}

sub get_sequencing_info_for_stock :Path('/ajax/genomes/sequenced_stocks') Args(1) {
    my $self = shift;
    my $c = shift;
    my $stock_id = shift;

    my $user_id = $c->user() ? $c->user->get_object()->get_sp_person_id() : undef;
    print STDERR "Retrieving sequencing info for stock $stock_id...\n";
    my @info = $self->retrieve_sequencing_infos($c->dbic_schema("Bio::Chado::Schema"), $user_id, $stock_id);

    # print STDERR "SEQ INFOS: ".Dumper(\@info);

    $c->stash->{rest} = { data => \@info };
}

sub retrieve_sequencing_infos {
    my $self = shift;
    my $schema = shift;
    my $user_id = shift;
    my @stock_ids = @_;

    print STDERR Dumper(\@stock_ids);

    my @data = ();

    foreach my $stock_id (uniq @stock_ids) {
	print STDERR "retrieving data for stock stock_id...\n";
	my $infos = CXGN::Stock::SequencingInfo->get_sequencing_project_infos($schema, $stock_id);

	# print STDERR "INFO = ".Dumper($infos);

	if ($infos) {
	    foreach my $info (@$infos) {

		my $website = "";
		if ($info->{website}) {
		    $website = qq | <a href="https://$info->{website}">$info->{website}</a> |;
		}

		my $blast_link = "BLAST";
		if ($info->{blast_link}) {
		    $blast_link = qq | <a href="$info->{blast_link}">BLAST</a> |;
		}

		my $jbrowse_link = "Jbrowse";
		if ($info->{jbrowse_link}) {
		    $jbrowse_link = qq | <a href="$info->{jbrowse_link}">JBrowse</a> |;
		}

		my $delete_link_js = "window.jsMod['sequenced_accessions'].delete_sequencing_info(".$info->{stockprop_id}.");";

		my $edit_link_js = "window.jsMod['sequenced_accessions'].edit_sequencing_info(".$info->{stockprop_id}.");";
		my $edit_delete_html = "Edit | Delete";

		if ($user_id && ($info->{sp_person_id} == $user_id)) {
		    $edit_delete_html = '<a href="javascript:'.$edit_link_js.'">Edit</a> | <a href="javascript:'.$delete_link_js.'">Delete</a>';
		}


		push @data, [
		    "<a href=\"/stock/$stock_id/view\">".$info->{uniquename}."</a>",
		    $info->{sequencing_year},
		    $info->{organization},
		    $website,
		    $blast_link." | ".$jbrowse_link,
		    $edit_delete_html,
		];
	    }
	}
    }

    print STDERR "Sequencing Data for this stock: ".Dumper(\@data);
    return @data;
}

sub get_sequencing_info :Path('/ajax/genomes/sequencing_info') Args(1) {
    my $self = shift;
    my $c = shift;
    my $stockprop_id = shift;

    my $si = CXGN::Stock::SequencingInfo->new( { schema => $c->dbic_schema("Bio::Chado::Schema"), stockprop_id => $stockprop_id });

    my $hashref = $si->to_hashref();

    $c->stash->{rest} = { data => $hashref };
}


sub store_sequencing_info :Path('/ajax/genomes/store_sequencing_info') Args(0) {
    my $self =shift;
    my $c = shift;

    my $params = $c->req->params();

    # print STDERR "Params for store: ".Dumper($params);

    if (!$c->user()) {
	$c->stash->{rest} = { error => "You need to be logged in to add sequencing inforaiton" };
	return;
    }
    if (!$c->user()->check_roles("curator")) {
	$c->stash->{rest} = { error => "You need to be logged in as a curator to submit this information" };
	return;
    }

    if (!$params->{stockprop_id}) { $params->{stockprop_id} = undef; } # force it to undef if it is a ""

    my $si = CXGN::Stock::SequencingInfo->new( { schema => $c->dbic_schema("Bio::Chado::Schema") });

    my $timestamp = DateTime->now();
    $params->{sp_person_id} = $c->user()->get_object()->get_sp_person_id();
    $params->{timestamp} = $timestamp->ymd()." ".$timestamp->hms();
    $params->{stock_id} = $params->{sequencing_status_stock_id};
    $si->from_hash($params);

    eval {
	$si->store();
    };
    if ($@) {
	$c->stash->{rest} = { error => "Error. The operation could not be completed. ($@)" };
	return;
    }

    $c->stash->{rest} = { success => 1};
}


sub delete_sequencing_info :Path('/ajax/genomes/sequencing_info/delete') Args(1) {
    my $self = shift;
    my $c = shift;

    my $stockprop_id = shift;

    print STDERR "delete_sequencing_info...\n";

    if (!$c->user() && !$c->user()->check_roles("curator")) {
	$c->stash->{rest} = { error => "Log in required for sequencing info dele
tion." };
	return;
    }

    my $si = CXGN::Stock::SequencingInfo->new(
	{
	    schema => $c->dbic_schema("Bio::Chado::Schema"),
	    stockprop_id => $stockprop_id,
	});

    if ($si->sp_person_id() && $si->sp_person_id() != $c->user()->get_object()->get_sp_person_id()) {
	$c->stash->{rest} = { error => "You don't own this entry so it cannot be deleted." };
	return;
    }

    print STDERR "Starting delete of stockprop_id $stockprop_id...(in object: ".$si->stockprop_id()."), type_id =". $si->type_id()."\n";

    my $success;
    if ($si->stockprop_id()) {
	eval {
	    $success = $si->delete();
	};
	if ($@) {
	    print STDERR "An error occurred during deletion. Sorry.\n";
	    $c->stash->{rest} = { error => "An error occurred while deleting sequencing info. ($@)" };
	    return;
	}
    }

    if ($success) {
	$c->stash->{rest} = { success => 1 };
    }
    else {
	$c->stash->{rest} = { error => "An error occurred during deletion." }
    }
}



1;
