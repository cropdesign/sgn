use strict;

package SGN::Controller::AJAX::Stability;

use Moose;
use Data::Dumper;
use File::Temp qw | tempfile |;
use File::Slurp;
use File::Spec qw | catfile|;
use File::Basename qw | basename |;
use File::Copy;
use CXGN::Dataset;
use CXGN::Dataset::File;
use CXGN::Tools::Run;
use CXGN::Page::UserPrefs;
use CXGN::Tools::List qw/distinct evens/;
use CXGN::Blast::Parse;
use CXGN::Blast::SeqQuery;
use Cwd qw(cwd);

BEGIN { extends 'Catalyst::Controller::REST' }

__PACKAGE__->config(
    default   => 'application/json',
    stash_key => 'rest',
    map       => { 'application/json' => 'JSON', 'text/html' => 'JSON' },
    );


sub shared_phenotypes: Path('/ajax/stability/shared_phenotypes') : {
    my $self = shift;
    my $c = shift;
    my $dataset_id = $c->req->param('dataset_id');
    my $people_schema = $c->dbic_schema("CXGN::People::Schema");
    my $schema = $c->dbic_schema("Bio::Chado::Schema", "sgn_chado");
    my $ds = CXGN::Dataset->new(people_schema => $people_schema, schema => $schema, sp_dataset_id => $dataset_id);
    my $traits = $ds->retrieve_traits();
    my @trait_info;
    foreach my $t (@$traits) {
        my $tobj = CXGN::Cvterm->new({ schema=>$schema, cvterm_id => $t });
        push @trait_info, [ $tobj->cvterm_id(), $tobj->name()];
    }

    
    $c->tempfiles_subdir("stability_files");
    my ($fh, $tempfile) = $c->tempfile(TEMPLATE=>"stability_files/trait_XXXXX");
    $people_schema = $c->dbic_schema("CXGN::People::Schema");
    $schema = $c->dbic_schema("Bio::Chado::Schema", "sgn_chado");
    my $temppath = $c->config->{basepath}."/".$tempfile;
    my $ds2 = CXGN::Dataset::File->new(people_schema => $people_schema, schema => $schema, sp_dataset_id => $dataset_id, file_name => $temppath);
    my $phenotype_data_ref = $ds2->retrieve_phenotypes();

    print STDERR Dumper(@trait_info);
    $c->stash->{rest} = {
        options => \@trait_info,
        tempfile => $tempfile."_phenotype.txt",
#        tempfile => $file_response,
    };
}


sub extract_trait_data :Path('/ajax/stability/getdata') Args(0) {
    my $self = shift;
    my $c = shift;

    my $file = $c->req->param("file");
    my $trait = $c->req->param("trait");

    $file = basename($file);

    my $temppath = File::Spec->catfile($c->config->{basepath}, "static/documents/tempfiles/stability_files/".$file);
    print STDERR Dumper($temppath);

    my $F;
    if (! open($F, "<", $temppath)) {
  $c->stash->{rest} = { error => "Can't find data." };
  return;
    }

    my $header = <$F>;
    chomp($header);
    print STDERR Dumper($header);
    my @keys = split("\t", $header);
    print STDERR Dumper($keys[1]);
    for(my $n=0; $n <@keys; $n++) {
        if ($keys[$n] =~ /\|CO\_/) {
        $keys[$n] =~ s/\|CO\_.*//;
        }
    }
    my @data = ();

    while (<$F>) {
  chomp;

  my @fields = split "\t";
  my %line = {};
  for(my $n=0; $n <@keys; $n++) {
      if (exists($fields[$n]) && defined($fields[$n])) {
    $line{$keys[$n]}=$fields[$n];
      }
  }
    print STDERR Dumper(\%line);
  push @data, \%line;
    }

    $c->stash->{rest} = { data => \@data, trait => $trait};
}

sub generate_results: Path('/ajax/stability/generate_results') : {
    my $self = shift;
    my $c = shift;
    my $dataset_id = $c->req->param('dataset_id');
    my $trait_id = $c->req->param('trait_id');
    print STDERR $dataset_id;
    print STDERR $trait_id;
    $c->tempfiles_subdir("stability_files");
    my $stability_tmp_output = $c->config->{cluster_shared_tempdir}."/stability_files";
    mkdir $stability_tmp_output if ! -d $stability_tmp_output;
    my ($tmp_fh, $tempfile) = tempfile(
      "stability_download_XXXXX",
      DIR=> $stability_tmp_output,
    );

    my $pheno_filepath = $tempfile . "_phenotype.txt";
    

    my $people_schema = $c->dbic_schema("CXGN::People::Schema");
    my $schema = $c->dbic_schema("Bio::Chado::Schema", "sgn_chado");

    my $temppath = $stability_tmp_output . "/" . $tempfile;

    my $ds = CXGN::Dataset::File->new(people_schema => $people_schema, schema => $schema, sp_dataset_id => $dataset_id, file_name => $temppath);

    my $phenotype_data_ref = $ds->retrieve_phenotypes($pheno_filepath);

    my $newtrait = $trait_id;
    $newtrait =~ s/\s/\_/g;
    $newtrait =~ s/\//\_/g;
    $trait_id =~ tr/ /./;
    $trait_id =~ tr/\//./;

    my $AMMIFile = $tempfile . "_" . "AMMIFile.png";
    my $figure1file = $tempfile . "_" . "figure1.png";
    my $figure2file = $tempfile . "_" . "figure2.png";
    my $figure3file = $tempfile . "_" . "figure3.png";


    # $trait_id =~ tr/ /./;
    # $trait_id =~ tr/\//./;


    my $cmd = CXGN::Tools::Run->new({
            backend => $c->config->{backend},
            temp_base => $c->config->{cluster_shared_tempdir} . "/stability_files",
            queue => $c->config->{'web_cluster_queue'},
            do_cleanup => 0,
            # don't block and wait if the cluster looks full
            max_cluster_jobs => 1_000_000_000,
        });

        print STDERR Dumper $pheno_filepath;

    # my $job;
    $cmd->run_cluster(
            "Rscript ",
            $c->config->{basepath} . "/R/stability/ammi_script.R",
            $pheno_filepath,
            $trait_id,
            $figure1file,
            $figure2file,
            $figure3file,
            $AMMIFile
    );
    $cmd->alive;
    $cmd->is_cluster(1);
    $cmd->wait;

    my $resultfile = $AMMIFile.".results";
    my $error;
    my $lines;

    if (! -e $resultfile) { 
    $error = "The analysis could not be completed. The factors may not have sufficient numbers of levels to complete the analysis. Please choose other parameters."
    }
    else { 
    $lines = read_file($temppath.".results");
    }
    my $figure_path = $c->{basepath} . "./documents/tempfiles/stability_files/";
    copy($AMMIFile, $figure_path);
    copy($figure1file, $figure_path);
    copy($figure2file, $figure_path);
    copy($figure3file, $figure_path);

    my $AMMIFilebasename = basename($AMMIFile);
    my $AMMIFile_response = "/documents/tempfiles/stability_files/" . $AMMIFilebasename;
    
    my $figure1basename = basename($figure1file);
    my $figure1_response = "/documents/tempfiles/stability_files/" . $figure1basename;
    
    my $figure2basename = basename($figure2file);
    my $figure2_response = "/documents/tempfiles/stability_files/" . $figure2basename;

    my $figure3basename = basename($figure3file);
    my $figure3_response = "/documents/tempfiles/stability_files/" . $figure3basename;

    print $AMMIFile_response;
        
    $c->stash->{rest} = {
        AMMITable => $AMMIFile_response,
        figure1 => $figure1_response,
        figure2 => $figure2_response,
        figure3 => $figure3_response,
        dummy_response => $dataset_id
        # dummy_response2 => $trait_id,
    };
}

1

