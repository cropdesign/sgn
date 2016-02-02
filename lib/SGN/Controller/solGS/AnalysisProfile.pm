package SGN::Controller::solGS::AnalysisProfile;

use Moose;
use namespace::autoclean;
use File::Path qw / mkpath  /;
use File::Spec::Functions qw / catfile catdir/;
use File::Slurp qw /write_file read_file :edit prepend_file/;
use JSON;
use CXGN::Tools::Run;
use Try::Tiny;


BEGIN { extends 'Catalyst::Controller' }


sub check_user_login :Path('/solgs/check/user/login') Args(0) {
  my ($self, $c) = @_;

  my $user = $c->user();
  my $ret->{loggedin} = 0;

  if ($user) 
  { 
      my $salutation = $user->get_salutation();
      my $first_name = $user->get_first_name();
      my $last_name  = $user->get_last_name();
          
      $self->get_user_email($c);
      my $email = $c->stash->{user_email};

      $ret->{loggedin} = 1;
      my $contact = { 'name' => $first_name, 'email' => $email};
     
      $ret->{contact} = $contact;
  }
   
  $ret = to_json($ret);
       
  $c->res->content_type('application/json');
  $c->res->body($ret);  
  
}


sub save_analysis_profile :Path('/solgs/save/analysis/profile') Args(0) {
    my ($self, $c) = @_;
   
    my $analysis_profile = $c->req->params;
    $c->stash->{analysis_profile} = $analysis_profile;
   
    my $analysis_page = $analysis_profile->{analysis_page};
    $c->stash->{analysis_page} = $analysis_page;
   
    my $ret->{result} = 0;
   
    $self->save_profile($c);
    my $error_saving = $c->stash->{error};
    
    if (!$error_saving) 
    {
	$ret->{result} = 1;	
    }

    $ret = to_json($ret);
       
    $c->res->content_type('application/json');
    $c->res->body($ret);  
    
}


sub save_profile {
    my ($self, $c) = @_;
        
    $self->analysis_log_file($c);
    my $log_file = $c->stash->{analysis_log_file};

    $self->add_headers($c);

    $self->format_profile_entry($c);
    my $formatted_profile = $c->stash->{formatted_profile};
    
    write_file($log_file, {append => 1}, $formatted_profile);
   
}


sub add_headers {
  my ($self, $c) = @_;

  $self->analysis_log_file($c);
  my $log_file = $c->stash->{analysis_log_file};

  my $headers = read_file($log_file);
  
  unless ($headers) 
  {  
      $headers = 'User_name' . 
	  "\t" . 'User_email' . 
	  "\t" . 'Analysis_name' . 
	  "\t" . "Analysis_page" . 	 
	  "\t" . "Status" .
	  "\t" . "Submitted on" .
	  "\t" . "Arguments" .
	  "\n";

      write_file($log_file, $headers);
  }
  
}


sub index_log_file_headers {
   my ($self, $c) = @_;
   
   $self->analysis_log_file($c);
   my $log_file = $c->stash->{analysis_log_file};
   
   my @headers = split(/\t/, (read_file($log_file))[0]);
   
   my $header_index = {};
   my $cnt = 0;
   
   foreach my $header (@headers)
   {
       $header_index->{$header} = $cnt;
       $cnt++;
   }
  
   $c->stash->{header_index} = $header_index;

}


sub format_profile_entry {
    my ($self, $c) = @_; 
    
    my $profile = $c->stash->{analysis_profile};
    my $time    = POSIX::strftime("%m/%d/%Y %H:%M", localtime);
    my $entry   = join("\t", 
		       (
			$profile->{user_name}, 
			$profile->{user_email}, 
			$profile->{analysis_name}, 
			$profile->{analysis_page},
			'Submitted',
			$time,
			$profile->{arguments},
		       )
	);

    $entry .= "\n";
	
    $c->stash->{formatted_profile} = $entry; 

}


sub run_saved_analysis :Path('/solgs/run/saved/analysis/') Args(0) {
    my ($self, $c) = @_;
   
    my $analysis_profile = $c->req->params;
    $c->stash->{analysis_profile} = $analysis_profile;

    $self->parse_arguments($c);
    
    $self->run_analysis($c);  
     
    $self->structure_output_details($c); 
    
    my $output_details = $c->stash->{bg_job_output_details};
      
    $c->stash->{r_temp_file} = 'analysis-status';
    $c->controller('solGS::solGS')->create_cluster_acccesible_tmp_files($c);
    my $out_temp_file = $c->stash->{out_file_temp};
    my $err_temp_file = $c->stash->{err_file_temp};
   
    my $temp_dir = $c->stash->{solgs_tempfiles_dir};
    my $status;
   
    try 
    { 
        my $job = CXGN::Tools::Run->run_cluster_perl({
           
            method        => ["solGS::AnalysisReport" => "check_analysis_status"],
    	    args          => [$output_details],
    	    load_packages => ['solGS::AnalysisReport'],
    	    run_opts      => {
    		              out_file    => $out_temp_file,
			      err_file    => $err_temp_file,
    		              working_dir => $temp_dir,
			      max_cluster_jobs => 1_000_000_000,
	    },
	    
         });
	
    }
    catch 
    {
	$status = $_;
	$status =~ s/\n at .+//s;           
    };

    if (!$status) 
    { 
	$status = $c->stash->{status}; 
    }
   
    my $ret->{result} = $status;	

    $ret = to_json($ret);
       
    $c->res->content_type('application/json');
    $c->res->body($ret);  

} 


sub parse_arguments {
  my ($self, $c) = @_;
 
  my $analysis_data =  $c->stash->{analysis_profile};
  my $arguments = $analysis_data->{arguments};

  if ($arguments) 
  {
      my $json = JSON->new();
      $arguments = $json->decode($arguments);
      
      foreach my $k ( keys %{$arguments} ) 
      {
	  if ($k eq 'population_id') 
	  {
	      my @pop_ids = @{ $arguments->{$k} };
	      $c->stash->{pop_ids} = \@pop_ids;
	      
	      if (scalar(@pop_ids) == 1) 
	      {		  
		  $c->stash->{pop_id}  = $pop_ids[0];
	      }
	  }

	  if ($k eq 'combo_pops_id') 
	  {
	      my $combo_pops_id = $arguments->{$k};
	      $c->stash->{combo_pops_id} = $arguments->{$k};
	  }

	  if ($k eq 'combo_pops_list') 
	  {
	      my @pop_ids = @{ $arguments->{$k} };
	      $c->stash->{combo_pops_list} = \@pop_ids;
	      
	      if (scalar(@pop_ids) == 1) 
	      {		  
		  $c->stash->{pop_id}  = $pop_ids[0];
	      }
	  }

	  if ($k eq 'trait_id') 
	  {
	      my @selected_traits = @{ $arguments->{$k} };
	      $c->stash->{selected_traits} = \@selected_traits;
	  } 
	  
	  if ($k eq 'analysis_type') 
	  {
	      $c->stash->{analysis_type} = $arguments->{$k};
	  }	 

	  if ($k eq 'data_set_type') 
	  {
	      $c->stash->{data_set_type} =  $arguments->{$k};
	  }	 
      }
  }
	    
}


sub structure_output_details {
    my ($self, $c) = @_;

    my $analysis_data =  $c->stash->{analysis_profile};
    my $arguments = $analysis_data->{arguments};
 
    $self->parse_arguments($c);
   
    my @traits_ids;
    
    if ($c->stash->{selected_traits}) 
    {
	@traits_ids = @{$c->stash->{selected_traits}};
    }
   
    my $pop_id        = $c->stash->{pop_id}; 
    my $combo_pops_id = $c->stash->{combo_pops_id};
    
    my $base = $c->req->base;    
    if ( $base !~ /localhost/)
    {
	$base =~ s/:\d+//;    
    } 
           
    my %output_details = ();

    my $solgs_controller = $c->controller('solGS::solGS');
    my $analysis_page = $analysis_data->{analysis_page};

    $analysis_page =~ s/$base/\//;

    if ($analysis_page =~ m/(solgs\/analyze\/traits\/|solgs\/trait\/|solgs\/model\/combined\/trials\/)/) 
    {	
	foreach my $trait_id (@traits_ids)
	{	    
	    $c->stash->{cache_dir} = $c->stash->{solgs_cache_dir};

	    $solgs_controller->get_trait_details($c, $trait_id);	    
	    $solgs_controller->gebv_kinship_file($c);
	 
	    my$trait_abbr = $c->stash->{trait_abbr};
	    my $trait_page;
	    my $referer = $c->req->referer;   
	    
	    if ( $referer =~ m/solgs\/population\// ) 
	    {
		$trait_page = $base . "solgs/trait/$trait_id/population/$pop_id";
	    }
	    
	    if ( $referer =~ m/solgs\/search\/trials\/trait/ ) 
	    {
		$trait_page = $base . "solgs/trait/$trait_id/population/$pop_id";
	    }
	    
	    if ( $referer =~ m/solgs\/populations\/combined\// ) 
	    {
		$trait_page = $base . "solgs/model/combined/trials/$pop_id/trait/$trait_id";
	    }

	    if ( $analysis_page =~ m/solgs\/model\/combined\/trials\// ) 
	    {
		$trait_page = $base . "solgs/model/combined/trials/$pop_id/trait/$trait_id";

		$c->stash->{combo_pops_id} = $pop_id;
		$solgs_controller->cache_combined_pops_data($c);		
	    }
	    
	    $output_details{$trait_abbr} = {
		'trait_id'       => $trait_id, 
		'trait_name'     => $c->stash->{trait_name}, 
		'trait_page'     => $trait_page,
		'gebv_file'      => $c->stash->{gebv_kinship_file},
		'pop_id'         => $pop_id,
		'phenotype_file' => $c->stash->{trait_combined_pheno_file},
		'genotype_file'  => $c->stash->{trait_combined_geno_file},
		'data_set_type'  => $c->stash->{data_set_type},
	    };
	}

    }
    elsif ( $analysis_page =~ m/solgs\/population\// ) 
    {
	my $population_page = $base . "solgs/population/$pop_id";

	$c->stash->{pop_id} = $pop_id;

	$solgs_controller->phenotype_file($c);	
	$solgs_controller->genotype_file($c);
	$solgs_controller->get_project_details($c, $pop_id);

	$output_details{$pop_id} = {
		'population_page' => $population_page,
		'population_id'   => $pop_id,
		'population_name' => $c->stash->{project_name},
		'phenotype_file'  => $c->stash->{phenotype_file},
		'genotype_file'   => $c->stash->{genotype_file},  
		'data_set_type'   => $c->stash->{data_set_type},
	};		
    }
    elsif ($analysis_page =~ m/solgs\/populations\/combined\//) 
    {
	my $combined_pops_page = $base . "solgs/populations/combined/$combo_pops_id";

	$solgs_controller->multi_pops_pheno_files($c);	
	$solgs_controller->multi_pops_geno_files($c);

	my $ph_files = $c->stash->{multi_pops_pheno_files};
	my @pheno_files = split(/\t/, \$ph_files);

	my $match_status = $c->stash->{pops_with_no_genotype_match};
	
	my @combined_pops_ids = @{$c->stash->{combo_pops_list}};

	foreach my $pop_id (@combined_pops_ids) 
	{	    
	    $solgs_controller->get_project_details($c, $pop_id);
	    my $population_name = $c->stash->{project_name};
	    
	    $c->stash->{pop_id} = $pop_id;
	    $solgs_controller->phenotype_file($c);
	    my $pheno_file = $c->stash->{phenotype_file};

	    $solgs_controller->genotype_file($c);
	    my $geno_file = $c->stash->{genotype_file};

	    my $population_page = $base . "solgs/population/$pop_id";

	    $output_details{$pop_id} = {
		'population_page'   => $population_page,
		'population_id'     => $pop_id,
		'population_name'   => $population_name,
		'combo_pops_id'     => $combo_pops_id,	
		'phenotype_file'    => $pheno_file,
		'genotype_file'     => $geno_file,
		'data_set_type'     => $c->stash->{data_set_type},
	    };	    
	}
	
	$output_details{no_match}           = $match_status;
	$output_details{combined_pops_page} = $combined_pops_page;
    }

    $self->analysis_log_file($c);
    my $log_file = $c->stash->{analysis_log_file};

    $output_details{analysis_profile}  = $analysis_data;
    $output_details{r_job_tempdir}     = $c->stash->{r_job_tempdir};
    $output_details{contact_page}      = $base . 'contact/form';
    $output_details{data_set_type}     = $c->stash->{data_set_type};
    $output_details{analysis_log_file} = $log_file;
   
    $c->stash->{bg_job_output_details} = \%output_details;
   
}


sub run_analysis {
    my ($self, $c) = @_;
 
    my $analysis_profile = $c->stash->{analysis_profile};
    my $analysis_page    = $analysis_profile->{analysis_page};

    my $base =   $c->req->base;
    $analysis_page =~ s/$base/\//;
   
    $c->stash->{background_job} = 1;
    
    my @selected_traits = @{$c->stash->{selected_traits}} if $c->stash->{selected_traits};
 
    if ($analysis_page =~ /solgs\/analyze\/traits\//) 
    {  	
	$c->controller('solGS::solGS')->build_multiple_traits_models($c);	
    } 
    elsif ($analysis_page =~  /solgs\/models\/combined\/trials\// )
    {
	if ($c->stash->{data_set_type} =~ /combined populations/)
	{
	    $c->stash->{combo_pops_id} = $c->stash->{pop_id};
	   
	    foreach my $trait_id (@selected_traits)		
	    {		
		$c->controller('solGS::solGS')->get_trait_details($c, $trait_id);   	
		$c->controller('solGS::combinedTrials')->combine_data_build_model($c);
	    }
	}
    }
    elsif ($analysis_page =~ /solgs\/model\/combined\/trials\// )	  
    {
	$c->stash->{combo_pops_id} = $c->stash->{pop_id};
	my $trait_id = $c->stash->{selected_traits}->[0];		
	
	$c->controller('solGS::solGS')->get_trait_details($c, $trait_id);
	$c->controller('solGS::combinedTrials')->combine_data_build_model($c);
       
    }
    elsif ($analysis_page =~ /solgs\/trait\//) 
    {
	$c->stash->{trait_id} = $selected_traits[0];
	$c->controller('solGS::solGS')->build_single_trait_model($c);
    }
    elsif ($analysis_page =~ /solgs\/population\//)
    {
	$c->controller('solGS::solGS')->phenotype_file($c);
	$c->controller('solGS::solGS')->genotype_file($c);
    }
    elsif ($analysis_page =~ /solgs\/populations\/combined\//)
    {
	my $combo_pops_id = $c->stash->{combo_pops_id};

	$c->controller("solGS::combinedTrials")->prepare_multi_pops_data($c);	
	
	$c->controller('solGS::solGS')->multi_pops_geno_files($c);
	my $g_files = $c->stash->{multi_pops_geno_files};
	my @geno_files = split(/\t/, $g_files);
	$c->controller('solGS::solGS')->compare_genotyping_platforms($c, \@geno_files);
    }
    else 
    {
	$c->stash->{status} = 'Error';
	print STDERR "\n I don't know what to analyze.\n";
    }

    my @error = @{$c->error};
    
    if ($error[0]) 
    {
	$c->stash->{status} = 'Failed submitting';
    }
    else 
    {    
	$c->stash->{status} = 'Submitted';
    }
 
    $self->update_analysis_progress($c);
 
}


sub update_analysis_progress {
    my ($self, $c) = @_;
     
    my $analysis_data =  $c->stash->{analysis_profile};
    my $analysis_name= $analysis_data->{analysis_name};
    my $status = $c->stash->{status};
    
    $self->analysis_log_file($c);
    my $log_file = $c->stash->{analysis_log_file};
  
    my @contents = read_file($log_file);
   
    map{ $contents[$_] =~ m/\t$analysis_name\t/
	     ? $contents[$_] =~ s/error|submitted/$status/ig 
	     : $contents[$_] } 0..$#contents; 
   
    write_file($log_file, @contents);

}


sub get_user_email {
    my ($self, $c) = @_;
   
    my $user = $c->user();

    my $private_email = $user->get_private_email();
    my $public_email  = $user->get_contact_email();
     
    my $email = $public_email 
	? $public_email 
	: $private_email;

    $c->stash->{user_email} = $email;

}


sub analysis_log_file {
    my ($self, $c) = @_;
      
    $self->create_analysis_log_dir($c);   
    my $log_dir = $c->stash->{analysis_log_dir};
    
    $c->stash->{cache_dir} = $log_dir;

    my $cache_data = {
	key       => 'analysis_log',
	file      => 'analysis_log',
	stash_key => 'analysis_log_file'
    };

    $c->controller('solGS::solGS')->cache_file($c, $cache_data);

}


sub confirm_request :Path('/solgs/confirm/request/') Args(0) {
    my ($self, $c) = @_;
    
    my $referer = $c->req->referer;
    
    $c->stash->{message} = "<p>Your analysis is running.</p>
                            <p>You will receive an email when it is completed.</p>
                            <p><a href=\"$referer\">[ Go back ]</a></p>";

    $c->stash->{template} = "/generic_message.mas"; 

}


sub display_analysis_status :Path('/solgs/display/analysis/status') Args(0) {
    my ($self, $c) = @_;
    
    my @panel_data = $self->solgs_analysis_status_log($c);

    my $ret->{data} = \@panel_data;
    
    $ret = to_json($ret);
       
    $c->res->content_type('application/json');
    $c->res->body($ret);  
    
}


sub solgs_analysis_status_log {
    my ($self, $c) = @_;
    
    $self->analysis_log_file($c);
    my $log_file = $c->stash->{analysis_log_file};
 
    my $ret = {};
    my @panel_data;
   
    if ($log_file)
    {    
	my @user_analyses = grep{$_ !~ /User_name\s+/i }
	                    read_file($log_file);

	$self->index_log_file_headers($c);
	my $header_index = $c->stash->{header_index};
	
	foreach my $row (@user_analyses) 
	{
	    my @analysis = split(/\t/, $row);
	    
	    my $analysis_name   = $analysis[$header_index->{'Analysis_name'}];
	    my $result_page     = $analysis[$header_index->{'Analysis_page'}];
	    my $analysis_status = $analysis[$header_index->{'Status'}];
	    my $submitted_on    = $analysis[$header_index->{'Submitted on'}];

	    if ($analysis_status =~ /(Failed|Submitted)/i) 
	    {
		$result_page = 'N/A';
	    }
	    else 
	    {
		$result_page = qq | <a href=$result_page>[ View ]</a> |;
	    }

	    push @panel_data, [$analysis_name, $submitted_on, $analysis_status, $result_page];
	}		
    }
 
    return \@panel_data;
}


sub create_analysis_log_dir {
    my ($self, $c) = @_;
        
    my $user_id = $c->user->id;
      
    $c->controller('solGS::solGS')->get_solgs_dirs($c);

    my $log_dir = $c->stash->{analysis_log_dir};

    $log_dir = catdir($log_dir, $user_id);
    mkpath ($log_dir, 0, 0755);

    $c->stash->{analysis_log_dir} = $log_dir;
  
}


sub begin : Private {
    my ($self, $c) = @_;

    $c->controller("solGS::solGS")->get_solgs_dirs($c);
  
}




__PACKAGE__->meta->make_immutable;


####
1;
####
