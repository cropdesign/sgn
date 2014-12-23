
=head1 NAME

CXGN::Trial - helper class for trials

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=head1 METHODS

=cut

package CXGN::Trial;


use Moose;
use Try::Tiny;

=head2 bcs_schema()

accessor for bcs_schema. Needs to be set when calling the constructor.

=cut

has 'bcs_schema' => ( isa => "Ref",
		      is => 'rw',
		      required => 1,
    );


=head2 trial_id()

accessor for the trial_id. Needs to be set when calling the constructor.

=cut

has 'trial_id' => (isa => 'Int',
		   is => 'rw',
		   reader => 'get_trial_id',
		   writer => 'set_trial_id',
    );


has 'layout' => (isa => 'CXGN::Trial::TrialLayout',
		 is => 'rw',
		 reader => 'get_layout',
		 writer => 'set_layout',
		 predicate => 'has_layout',
		 

    );


=head2 get_year()

getter for the year property.

=cut

sub get_year { 
    my $self = shift;

    my $type_id = $self->get_year_type_id();

    my $rs = $self->bcs_schema->resultset('Project::Project')->search( { 'me.project_id' => $self->get_trial_id() })->search_related('projectprops', { type_id => $type_id } );

    if ($rs->count() == 0) { 
	return undef;
    }
    else { 
	return $rs->first()->value();
    }
}

=head2 set_year()

setter for the year property.

=cut

sub set_year { 
    my $self = shift;
    my $year = shift;
    
    my $type_id = $self->get_year_type_id();
    
    my $row = $self->bcs_schema->resultset('Project::Projectprop')->find( { project_id => $self->get_trial_id(), type_id => $type_id  });

    if ($row) { 
	$row->value($year);
	$row->update();
    }
    else { 
	$row = $self->bcs_schema->resultset('Project::Projectprop')->create(
	    { type_id => $type_id,
	    value => $year,
	    } );
    }
}

=head2 get_description()

getter for the description

=cut

sub get_description { 
    my $self = shift;

    my $rs = $self->bcs_schema->resultset('Project::Project')->search( { project_id => $self->get_trial_id() });

    return $rs->first()->description();

}


sub set_description { 
    my $self = shift;
    my $description = shift;
    
    my $row = $self->bcs_schema->resultset('Project::Project')->find( { project_id => $self->get_trial_id() });

    print STDERR "Setting new description $description for trial ".$self->get_trial_id()."\n";

    $row->description($description);

    $row->update();

}


sub get_location { 
    my $self = shift;

    if ($self->get_location_type_id()) { 
	my $row = $self->bcs_schema->resultset('Project::Projectprop')->find( { project_id => $self->get_trial_id() , type_id=> $self->get_location_type_id() });
	
	if ($row) { 
	    my $loc = $self->bcs_schema->resultset('NaturalDiversity::NdGeolocation')->find( { nd_geolocation_id => $row->value() });
	    
	    return [ $row->value(), $loc->description() ];
	}
	else { 
	    return [];
	}
    }
    

}


sub add_location { 
    my $self = shift;
    my $location_id = shift;

    my $row = $self->bcs_schema->resultset('Project::Projectprop')->create( 
	{ 
	    project_id => $self->get_trial_id(),
	    type_id => $self->get_location_type_id(),
	    value => $location_id,
	});
    
    
}

sub remove_location { 
    my $self = shift;
    my $location_id = shift;
    
    my $row = $self->bcs_schema->resultset('Project::Projectprop')->find( 
	{ 
	    project_id => $self->get_trial_id(),
	    type_id => $self->get_location_type_id(),
	    value => $location_id,
	});
    if ($row) { 
	print STDERR "Removing location $location_id from trail ".$self->get_trial_id()."\n";
	$row->delete();
    }

}

# sub get_project_type { 
#     my $self = shift;
#     my $row = $self->bcs_schema->resulset('Project::Projectprop')->find( { project_id => $self->get_trial_id() , type_id=> $self->get_location_type_id() });
    
#     return $row->value();
    

# }


sub set_project_type { 
    

}

sub get_project_type { 
    my $self = shift;
    my $row = $self->bcs_schema->resultset('Cv::Cv')->find( { name => 'project_types' } );

    my @types;
    if ($row) { 
	my $rs = $self->bcs_schema->resultset('Project::Projectprop')->search( { project_id => $self->get_trial_id() })->search_related('type', { cv_id => $row->cv_id() });
	foreach my $r ($rs->next()) { 
	    push @types, [ $r->cvterm_id(), $r->name() ];
	}
	
	return @types;
    }
	
    return ();

}

sub get_location_type_id { 
    my $self = shift;
    my $rs = $self->bcs_schema->resultset('Cv::Cvterm')->search( { name => 'project location' });

    if ($rs->count() > 0) { 
	return $rs->first()->cvterm_id();
    }

}

sub get_year_type_id { 
    my $self = shift;

    my $rs = $self->bcs_schema->resultset('Cv::Cvterm')->search( { name => 'project year' });

    return $rs->first()->cvterm_id();
}

sub get_name { 
    my $self = shift;
    my $row = $self->bcs_schema->resultset('Project::Project')->find( { project_id => $self->get_trial_id() });
    
    if ($row) { 
	return $row->name();
    }
}
 
sub set_name { 

}   

sub delete_phenotype_data { 
    my $self = shift;

    my $trial_id = $self->get_trial_id();

    eval { 
	$self->bcs_schema->txn_do( 
	    sub { 
		# first, delete metadata entries
		#
		$self->delete_metadata($trial_id);
		
		# delete phenotype data associated with trial
		#
		my $trial = $self->bcs_schema()->resultset("Project::Project")->search( { project_id => $trial_id });
		
		my $nd_experiment_rs = $self->bcs_schema()->resultset("NaturalDiversity::NdExperimentProject")->search( { project_id => $trial_id });
		my @nd_experiment_ids = map { $_->nd_experiment_id } $nd_experiment_rs->all();
		
		$self->_delete_phenotype_experiments(@nd_experiment_ids); # cascading deletes should take care of everything (IT DOESNT????)
		
	    });
    };
    if ($_) { 
	return "Error deleting phenotype data for trial $trial_id. $_\n";
    }
    return '';
    
}
    
sub delete_field_layout { 
    my $self = shift;

    my $trial_id = $self->get_trial_id();
    
    # Note: metadata entries need to be deleted separately using delete_metadata()
    #
    my $error = '';
    eval { 
	$self->bcs_schema()->txn_do( 
	    sub {
		my $trial = $self->bcs_schema()->resultset("Project::Project")->search( { project_id => $trial_id });
		
		my $nd_experiment_rs = $self->bcs_schema()->resultset("NaturalDiversity::NdExperimentProject")->search( { project_id => $trial_id });
		my @nd_experiment_ids = map { $_->nd_experiment_id } $nd_experiment_rs->all();
		
		$self->_delete_field_layout_experiment($trial_id); 
	    }
	    );
    };
    if ($_) { 
	print STDERR "ERROR $_\n";
	return "An error occurred: $_\n";
    }
    
    return '';
}

sub delete_metadata { 
    my $self = shift;
    my $metadata_schema = shift;

    if (!$metadata_schema) { die "Need metadata schema parameter\n"; }

    my $trial_id = $self->get_trial_id();

    # first, deal with entries in the md_metadata table, which may reference nd_experiment (through linking table)
    my $q = "SELECT distinct(metadata_id) FROM nd_experiment_project JOIN phenome.nd_experiment_md_files using(nd_experiment_id) JOIN metadata.md_files using(file_id) JOIN metadata.md_metadata using(metadata_id) WHERE project_id=?";
    my $h = $self->bcs_schema->storage()->dbh()->prepare($q);
    $h->execute($trial_id);

    while (my ($md_id) = $h->fetchrow_array()) { 
	my $mdmd_row = $self->metadata_schema->resultset("MdMetadata")->find( { metadata_id => $md_id } );
	if ($mdmd_row) { 
	    $mdmd_row -> update( { obsolete => 1 });
	}
    }

    # delete the entries from the linking table...
    $q = "SELECT distinct(file_id) FROM nd_experiment_project JOIN phenome.nd_experiment_md_files using(nd_experiment_id) JOIN metadata.md_files using(file_id) JOIN metadata.md_metadata using(metadata_id) WHERE project_id=?";
    $h = $self->bcs_schema->storage()->dbh()->prepare($q);
    $h->execute($trial_id);
    
    while (my ($file_id) = $h->fetchrow_array()) { 
	my $ndemdf_rs = $self->phenome_schema->resultset("NdExperimentMdFiles")->search( { file_id=>$file_id });

	foreach my $row ($ndemdf_rs->all()) { 
	    $row->delete();
	}
    }
}


sub _delete_phenotype_experiments { 
    my $self = shift;
    my @nd_experiment_ids = @_;

    # retrieve the associated phenotype ids (they won't be deleted by the cascade)
    #
    my $phenotypes_deleted = 0;
    my $nd_experiments_deleted = 0;

    my $phenotype_rs = $self->bcs_schema()->resultset("NaturalDiversity::NdExperimentPhenotype")->search( { nd_experiment_id=> { -in => [ @nd_experiment_ids ] }}, { join => 'phenotype' });
    if ($phenotype_rs->count() > 0) { 
	foreach my $p ($phenotype_rs->all()) { 
	    $p->delete();
	    $phenotypes_deleted++;
	}
    }
    
    # delete the experiments
    #
    my $delete_rs = $self->bcs_schema()->resultset("NaturalDiversity::NdExperiment")->search({ nd_experiment_id => { -in => [ @nd_experiment_ids] }});
    $nd_experiments_deleted = $delete_rs->count();
    $delete_rs->delete_all();

    return { phenotypes_deleted => $phenotypes_deleted, 
	     nd_experiments_deleted => $nd_experiments_deleted
    };
}

=head2 _delete_field_layout_experiment

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub _delete_field_layout_experiment { 
    my $self = shift;
    my $trial_id = shift;

    # check if there are still associated phenotypes...
    #
    if ($self->trial_has_phenotype_data()) { 
	print STDERR "Attempt to delete field layout that still has associated phenotype data.\n";
	return { error => "Trial still has associated phenotyping experiment, cannot delete." };
    }

    my $field_layout_type_id = $self->bcs_schema->resultset("Cv::Cvterm")->find( { name => "field layout" })->cvterm_id();
    #print STDERR "Field layout type id = $field_layout_type_id\n";

    my $plot_type_id = $self->bcs_schema->resultset("Cv::Cvterm")->find( { name => 'plot' })->cvterm_id();
    #print STDERR "Plot type id = $plot_type_id\n";

    my $q = "SELECT stock_id FROM nd_experiment_project JOIN nd_experiment USING (nd_experiment_id) JOIN nd_experiment_stock ON (nd_experiment.nd_experiment_id = nd_experiment_stock.nd_experiment_id) JOIN stock USING(stock_id) WHERE nd_experiment.type_id=? AND project_id=? AND stock.type_id=?";
    my $h = $self->bcs_schema->storage()->dbh()->prepare($q);
    $h->execute($field_layout_type_id, $trial_id, $plot_type_id);

    my $plots_deleted = 0;
    while (my ($plot_id) = $h->fetchrow_array()) { 
	my $plot = $self->bcs_schema()->resultset("Stock::Stock")->find( { stock_id => $plot_id });
	#print STDERR "Deleting associated plot ".$plot->name()." (".$plot->stock_id().") \n";
	$plots_deleted++;
	$plot->delete();
    }

    $q = "SELECT nd_experiment_id FROM nd_experiment JOIN nd_experiment_project USING(nd_experiment_id) WHERE nd_experiment.type_id=? AND project_id=?";
    $h = $self->bcs_schema->storage()->dbh()->prepare($q);
    $h->execute($field_layout_type_id, $trial_id);
    
    my ($nd_experiment_id) = $h->fetchrow_array();
    if ($nd_experiment_id) { 
	#print STDERR "Delete corresponding nd_experiment entry  ($nd_experiment_id)...\n";
	my $nde = $self->bcs_schema()->resultset("NaturalDiversity::NdExperiment")->find( { nd_experiment_id => $nd_experiment_id });
	$nde->delete();
    }


    #return { success => $plots_deleted };
    #return { success => 1 };
}

sub trial_has_phenotype_data { 
    my $self = shift;
    my $trial_id = shift;

    my $phenotyping_experiment_type_id = $self->bcs_schema->resultset("Cv::Cvterm")->find( { name => 'phenotyping experiment' })->cvterm_id();
    
    my $phenotype_experiment_rs = $self->bcs_schema()->resultset("NaturalDiversity::NdExperimentProject")->search( 
	{ 
	    project_id => $trial_id, 'nd_experiment.type_id' => $phenotyping_experiment_type_id}, 
	{ 
	    join => 'nd_experiment'  
	}
	);
    
    return $phenotype_experiment_rs->count();

}



1;
