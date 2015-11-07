
package SGN::Controller::InSilicoPcr;

use Moose;
use File::Temp qw | tempfile |;
use File::Basename qw | basename |;
use Storable qw | nstore retrieve |;
use Bio::Graphics::Gel;
use File::NFSLock qw | uncache |;
use File::Copy qw | copy |;
use Data::Dumper;

use Config::Any;
use Try::Tiny;
use Tie::UrlEncoder; our %urlencode;
use File::Spec qw | catfile |;
use File::Slurp qw | read_file write_file |;
use CXGN::Tools::Run;
use CXGN::BlastDB;
use Bio::Seq;
use CXGN::Page::FormattingHelpers qw/ html_break_string/;


BEGIN { extends 'Catalyst::Controller'; }

sub AUTO { 
    my $self = shift;
    my $c = shift;
    SGN::Schema::BlastDb->dbpath($c->config->{blast_db_path});
}

sub index :Path('/tools/in_silico_pcr/') :Args(0) { 
  my $self = shift;
  my $c = shift;

  my $db_id = $c->req->param('db_id');

  my $seq = $c->req->param('seq');
  my $schema = $c->dbic_schema("SGN::Schema");
  
  my $group_rs = $schema->resultset("BlastDbGroup")->search( name => "Genome Sequences");
  # my $group_rs = $schema->resultset("BlastDbGroup")->search( undef, { order_by=>'ordinal' });

  my $databases = {};
  my $dataset_groups = [];
  
  my $preselected_database = $c->config->{preselected_blastdb};
  my $preselected_category = '';
  
  if ($db_id) { 
    my $rs = $schema->resultset("BlastDb")->search( { blast_db_id => $db_id }, { join => 'blast_db_group' });
    
    if ($rs == 0) {
      $c->throw( is_error => 0, message => "The blast database with id $db_id could not be found.");
    }
    
    $preselected_database = $rs->first()->blast_db_id(); # first database of the category
    $preselected_category = $rs->first()->blast_db_group_id();
  }
    
  foreach my $g ($group_rs->all()) { 
    my @blast_dbs = $g->blast_dbs();
    push @$dataset_groups, [ $g->blast_db_group_id, $g->name() ];
    
    my @dbs_AoA;

    foreach my $db (@blast_dbs) {
      push @dbs_AoA, [ $db->blast_db_id(), $db->title(), $db->type() ];
    }

    my @arr = sort {$a->[1] cmp $b->[1]} @dbs_AoA;
    $databases->{ $g->blast_db_group_id } = \@arr;
  }

  $c->stash->{preselected_database} = $preselected_database;
  $c->stash->{preselected_category} = $preselected_category;
  $c->stash->{preload_id} = $c->req->param('preload_id');
  $c->stash->{preload_type} = $c->req->param('preload_type');

  $c->stash->{databases} = $databases;
  $c->stash->{dataset_groups} = $dataset_groups;
  
  $c->stash->{template} = '/tools/in_silico_pcr/in_silico_pcr.mas';
}


sub _reverse_complement{
	my $seq = shift;
	my $rev_seq = reverse $seq;
  $rev_seq =~ tr/ACGTacgt/TGCAtgca/;
  
	return $rev_seq;
}


sub run_pcr_blast :Path('/tools/pcr_results') :Args(0) {
  my ($self, $c) = @_;
  
  my @errors; #to store erros as they happen
  
  my $fprimer = $c->req->param("fprimer");
  my $rprimer = $c->req->param("rprimer");
  my $productLength = $c->req->param("productLength");
  my $allowedMismatches = $c->req->param('allowedMismatches');
  my $frevcom = $c->req->param('frevcom'); #forward primer reverse complement
  my $rrevcom = $c->req->param('rrevcom'); #reverse primer reverse complement
  
  my $params = $c->req->params();

  #getting the length of the primers
  my $flength = length($fprimer);
  my $rlength = length($rprimer);
  $params->{flength} = $flength;
  $params->{rlength} = $rlength;

  #reverse complement if checked
  if ($frevcom){
  	$fprimer = _reverse_complement($fprimer);
  }
  if ($rrevcom){
  	$rprimer = _reverse_complement($rprimer);
  }

  # print STDERR "fprimer: $fprimer, rprimer: $rprimer\n";
  # print STDERR "DB id: $blast_db_id\n";

  
  #validating the primers input
  if  (!$fprimer){
      push ( @errors , "Forward Primer was not provided!\n");
  }
  elsif ($fprimer =~ /[^a-zA-Z]/g){
       push (  @errors , "Forward Primer Can only hold letters (no numbers are allowed)\n");
  }
  if (!$rprimer){
      push (  @errors , "Reverse Primer was not provided!\n");
  }
  elsif ($rprimer =~ /[^a-zA-Z]/g){
       push (  @errors , "Reverse Primer Can only hold letters (no numbers are allowed)\n");
  }

  #validating productLength
  push (@errors , "Max Product Length should be a positive digit\n")
  	if ($productLength <= 0 or $productLength !~ /^[\d]*$/g);

  #validating AllowedMismatches
  push (@errors , "Allowed mismatches should be a positive digit\n")
  	if ($allowedMismatches < 0 or $allowedMismatches !~ /^[\d]*$/g);
  
  
  # return errors
  if (scalar (@errors) > 0){
    $c->stash->{errors} = join("<BR>" , @errors);
    $c->stash->{template} = '/tools/in_silico_pcr/insilicopcr_output.mas';
  }
  
  # giving the primers a fasta format
  $fprimer = ">FORWARD-PRIMER\n$fprimer";
  $rprimer = ">REVERSE-PRIMER\n$rprimer";

  my $sequence = "$fprimer\n$rprimer\n";
  $sequence =~ s/^\s+|\s+$|\n\s*\n//g; #< trim out leading and trailing whitespace and blank lines
  
  # create blast input and output file
  my $blast_tmp_output = $c->config->{cluster_shared_tempdir}."/blast";
  mkdir $blast_tmp_output if ! -d $blast_tmp_output;

  my ($seq_fh, $seqfile) = tempfile( 
  	"blast_XXXXXX",
  	DIR=> $blast_tmp_output,
	);
  
  print STDERR "Opening file for sequence ($seqfile)... ";
  open(my $FH, ">", $seqfile) || die "Can't open file for query ($seqfile)\n";
  print $FH $sequence if $sequence;
  close($FH);

  my $jobid = basename($seqfile);
  print STDERR "JOB ID CREATED: $jobid\n";

  print STDERR "Done.\n";



  # get matrix
  my $m = $params->{matrix};
  $m =~ /^(BLOSUM|PAM)\d+$/ or $c->throw( is_error => 0, message => "invalid matrix '$m'" );
           
  # get evalue
  my $evalue = $params->{evalue} ? $params->{evalue} : 1;
	
  # get filter
  my $filter = $params->{filterq} ? 'T' : 'F';
	
  # get blast database
  my $schema = $c->dbic_schema("SGN::Schema");
  my $bdb = $schema->resultset("BlastDb")->find($params->{database} ) or die "could not find bdb with file_base '$params->{database}'";
  my $basename = File::Spec->catfile($c->config->{blast_db_path},$bdb->file_base());
  
  
  my $result_file =  $seqfile.".out";
  # my $result_file = File::Spec->catfile($c->config->{basepath}, $c->tempfiles_subdir('blast'), $jobid.".out");
  
  # create blast command
  my $blast_cmd = "blastall -p blastn -i $seqfile -d $basename -m 8 -M $m -F $filter -e $evalue -o $result_file";
	
  print STDERR "COMMAND: $blast_cmd\n";
  my $blast_error = system($blast_cmd);
  
  if ($blast_error) { 
    print STDERR "An error occurred! $blast_error\n";
    $c->stash->{rest} = { error => $blast_error };
  }
  
  
  my $pcr_img_path = $c->config->{'tempfiles_subdir'};
  my ($pcr_seq,$gel_url) = _blast_to_pcr($result_file,$params,$pcr_img_path,$basename);
  
  
  

      
  
  
  
  $c->stash->{pcr_seq} = $pcr_seq;
  $c->stash->{gel_url} = $gel_url;
  $c->stash->{template} = '/tools/in_silico_pcr/insilicopcr_output.mas';
}

sub _blast_to_pcr {
  my $blast_file = shift;
  my $params = shift;
  my $gel_img_tempdir = shift;
  my $basename = shift;
  
  #Parsing the blast m8 result file

  #Report parsing method was taken from /util/blast_result_handle.pl by Aure 

  my (%query_id, %subject_id, %identity, %align_length, %mismatches, %gap_openings, %q_start, %q_end, %s_start, %s_end, 
        %e_value, %hit_score, %orientation);

  my (@fprimer_ids ,@rprimer_ids);
  
  open(my $res_fh, "<", $blast_file) or die "$! opening $blast_file for reading";
  
  my $line=0;

  while (<$res_fh>) {
      $line++;
    
      my @data=split(/\t/, $_);
    
      #separating forward primers from reverse primers using 2 arrays
      push (@fprimer_ids , $line) if ($data[0] eq 'FORWARD-PRIMER');
      push (@rprimer_ids , $line) if ($data[0] eq 'REVERSE-PRIMER');

      # print STDERR "$data[1]\n";
      
      $query_id{$line}=$data[0];
      $subject_id{$line}=$data[1];
      $identity{$line}=$data[2];
      $align_length{$line}=$data[3];
      $mismatches{$line}=$data[4];
      $gap_openings{$line}=$data[5];
      $q_start{$line}=$data[6];
      $q_end{$line}=$data[7];
      $s_start{$line}=$data[8];
      $s_end{$line}=$data[9];
      $e_value{$line}=$data[10];
      $hit_score{$line}=$data[11];

      #finding the orientation of the strand "+" is  5'->3' and "-" 3'->5'
      $orientation{$line}= ($s_end{$line}-$s_start{$line} > 0)? '+' : '-';
  }
 
  close $res_fh;
  
  
  #Finding Results
  my @pcr_results; #is an array of array references [forward Primer # , reverse primer #, + or - for parent orientation] 

  foreach my $forward (@fprimer_ids){
    
    print STDERR "fwd: ".$subject_id{$forward}."\t".$s_start{$forward}."\t".$align_length{$forward}."\t".$mismatches {$forward}."\n";
    print STDERR "params: product length: ".$params->{productLength}."\tmm: ".$params->{allowedMismatches}."\tfwd length: ".$params->{flength}."\trev length: ".$params->{rlength}."\n";
    
    foreach my $reverse (@rprimer_ids){
      print STDERR "rev: ".$subject_id{$reverse}."\t".$s_start{$reverse}."\t".$align_length{$reverse}."\t".$mismatches {$reverse}."\n";
	    
  		if ($subject_id{$forward} eq $subject_id{$reverse}    #both on the same subject seq
        and  $s_start{$reverse}- $s_start{$forward}<= $params->{productLength} #product Length is within user's choice
        and  $mismatches {$forward} <= $params->{allowedMismatches}  #Allowed mismatches by user
        and  $mismatches {$reverse} <= $params->{allowedMismatches}
        and  $align_length{$forward} == $params->{flength}  #primers match exact length
        and  $align_length {$reverse} == $params->{rlength}
  		) {
        print STDERR "\ninside the if\n";
        print STDERR "orientation fwd: ".$orientation{$forward}." s_end fwd: ".$s_end{$forward}."\n";
        print STDERR "orientation rev: ".$orientation{$reverse}." s_end rev: ".$s_end{$reverse}."\n";
        
          #if the product is in the + starnd of parent seq add a + sign in the array
        
          if ( $orientation{$forward} eq '+'     #forward is on the + strand
          and  $orientation{$reverse} eq '-'     #reverse is on the - strand b/c its a complement
          and  $s_end{$forward} < $s_end{$reverse}  #end of forward located upstream of beginning of reverse 
          ){
            print STDERR "first if";
          	push (@pcr_results , [$forward,$reverse, '+']) ;
           }
        	
          #if the product is in the - strand of the parent seq add a - sign in the array	
          elsif ( $orientation{$forward} eq '-'     #forward is on the - strand (complemet here)
             and  $orientation{$reverse} eq '+'     #reverse is on the + strand 
             # and  $s_end{$forward} < $s_end{$reverse}  #end of forward located downstream of beginning of reverse
             and  $s_end{$forward} > $s_end{$reverse}  #end of forward located downstream of beginning of reverse
            )
             {
               print STDERR "scond if\n";
                push (@pcr_results , [$forward,$reverse, '-']);
             }	
             else {
               print STDERR "\nin else\n";
             }
       }
    }#end of the 4each loop

  }
  
  
  
  
  ##############################################################################################################################
  
  my $fs = Bio::BLAST::Database->open(full_file_basename => "$basename",);
			
  
  

  # my $bdb = CXGN::BlastDB->from_id($params->{database});
  # die "No such database" if (!$bdb);
  #
  # my $basename = $bdb->full_file_basename;

  # print STDERR "bdb: $bdb\n";
  print STDERR "basename: $basename\n";

  my $find_seq;
  my $find_subseq;
  my $find_id;
  my $pcr_seq;
  
  # my $report_download_link = qq|[ <a href="$raw_report_url">BLAST OUTPUT</a> ]|;


  my @product_sizes; #for the gel




  # print info_section_html( title => 'PCR Report',
  #           subtitle => $report_download_link,
  #           contents => "\n");



  if (scalar(@pcr_results) ==  0 ){
    print STDERR "No PCR Product Found\n";
  }
  else{

    foreach my $result (@pcr_results){
	
  	#finding parent sequence
    $find_id = $subject_id{$result->[0]};
    $find_seq = $fs->get_sequence($find_id);
    # $find_seq = $bdb->get_sequence($find_id);
    print STDERR "id: $find_id\n";
    print STDERR "seq: $find_seq\n";
    
    print STDERR "s_sstart1: ".$s_start{$result->[1]}." s_start0: ".$s_start{$result->[0]}." res2: ".$result->[2]."\n";

  	#finding the pcr result sequence
  	$find_subseq = $find_seq->subseq($s_start{$result->[0]},$s_start{$result->[1]}) if $result->[2] eq '+';
    $find_subseq = $find_seq->subseq($s_start{$result->[1]},$s_start{$result->[0]}) if $result->[2] eq '-';
	
  	######################################################################################
	
  	#generating sequence object for the result to be able to find its molecular weight
  	my $seq_obj = Bio::Seq->new(-seq       => $find_subseq ,
                                  -alphabet  => 'dna' 
                                  );


  	my $seq_size = $seq_obj->length;
     	push (@product_sizes , $seq_size);
    
      #finding the ID of the sequence and adding + sign if it is on the plus strand and - if its on minus strand and some coordinates
  	$find_id = $find_seq->id();
  	$find_id .= $result->[2] eq '+' ? ' strand = plus, ' : ' strand = minus, ';
  	$find_id .= " start = $s_start{$result->[0]}, end = $s_start{$result->[1]}, size = $seq_size bp";
	 
  	#######################################################################################

  	#reverse complementing $find_subseq if the orientation is '-'
  	$pcr_seq = $seq_obj->revcom->seq if $result->[2] eq '-'; 
	
    # $find_subseq = html_break_string($find_subseq , 90);
    # for (my $i=0; $i<length($find_subseq); $i=$i+100) {
    #   $pcr_seq = $pcr_seq.substr($find_subseq,$i,100)."\n";
    # }
    # chomp($pcr_seq);

  	print STDERR " $find_id\n $pcr_seq\n";

  }

  ##############################################################################################################################
  #Generating a gel of the results 

      my $gel = Bio::Graphics::Gel->new('pcr' => \@product_sizes,             	
              	      -lane_length => 200,
              	      -bandcolor => [0xff,0xc6,0x00]);
    
      my $gel_img = $gel->img->png;
    

      #saving the gel img in a temp file 
      # my $gel_img_tempdir = $c->path_to( $c->tempfiles_subdir('temp_images') );
    
      my ($fh ,$temp_file) = tempfile( DIR => $gel_img_tempdir, TEMPLATE=>"gel_XXXXXX", SUFFIX => ".png");
      print $fh $gel_img;

      my $pcr_img_name = basename ($temp_file);

      #generating the url
      my $img_url = $temp_file;

      print STDERR "img_url: $img_url\n";
      print STDERR "base_temp: $pcr_img_name\n";
      
      return ($pcr_seq,$pcr_img_name);
  } 
  
  
  
}

1;
