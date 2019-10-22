
package CXGN::Analysis::AnalysisMetadata;

use Moose;

extends 'CXGN::JSONProp';

has 'dataset_id' => ( isa => 'Maybe[Int]', is => 'rw');

has 'dataset_data' => (isa => 'Maybe[Str]', is => 'rw');

has 'analysis_protocol' => (isa => 'Maybe[Str]', is => 'rw');

sub BUILD {
    my $self = shift;
    my $args = shift;

    $self->prop_table('projectprop');
    $self->prop_namespace('Project::Projectprop');
    $self->prop_primary_key('projectprop_id');
    $self->prop_type('analysis_metadata');
    $self->cv_name('project_property');
    $self->allowed_fields([ qw | dataset_id dataset_data analysis_protocol | ]);
    $self->parent_table('project');
    $self->parent_primary_key('project_id');

    $self->load();

}

