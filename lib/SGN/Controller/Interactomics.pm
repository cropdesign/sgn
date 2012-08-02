package SGN::Controller::Interactomics;


=head1 NAME

SGN::Controller::Interactomics - Provides and interface for loading of CytoScape via Web Start with interactomic data.


=cut


use Moose;
use URI::FromHash 'uri';


BEGIN { extends 'Catalyst::Controller'; }


=head2

   Loads SGN page for Web Start app.

=cut


sub interactomics :Path("/tools/interactomics") :Args(0) { 
    my ($self, $c) = @_;

    $c->stash->{template} = '/interactomics/index.mas';

}



=head2 create_jnlp

  Usage: In an html form, pass it the file (.cys or other cytoscape compatible) location with name "Network" and call this URL with the form action.  
  Desc : Uses the current URL base and a file passed as a parameter to create the body of a jnlp file for loading CytoScape via Web Start with preloaded data.
  Ret  : JNLP file that the user's browser should automatically download and run with Java Web Start.  
  Args : none
  Side Effects: ?

=cut

sub create_jnlp :Path("/tools/interactomics/coffee.jnlp") :Args(0) {
    my ($self, $c) = @_;
    
    my $codebase = ($c->req->base)."/CytoScape/";
    
    my $file = $c->req->param("Network");
    
    #Content of jnlp file ultimately downloaded by user
    #Change this if and only if the CytoScape build at /static/CytoScape is changed in any way.
    #codebase is specified on line 44 of this controller
    #file is specified on line 176 of this controller
    my $doc =  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<jnlp codebase=\"".$codebase."\" href=\"\">   
  <security>
    <all-permissions />
  </security>
  <information>
    <title>Cytoscape Webstart</title>
    <vendor>Cytoscape Collaboration</vendor>
    <homepage href=\"http://cytoscape.org\" />
    <offline-allowed />
  </information>
  <resources>
    <j2se version=\"1.5+\" max-heap-size=\"1024M\" />
    <!--All lib jars that cytoscape requires to run should be in this list-->
    <jar href=\"cytoscape.jar\" />
    <jar href=\"lib/cytoscape-geom-spacial.jar\" />
    <jar href=\"lib/freehep-util-2.0.2.jar\" />
    <jar href=\"lib/biojava-1.4.jar\" />
    <jar href=\"lib/resolver.jar\" />
    <jar href=\"lib/cytoscape-graph-dynamic.jar\" />
    <jar href=\"lib/colt.jar\" />
    <jar href=\"lib/phoebe.jar\" />
    <jar href=\"lib/activation.jar\" />
    <jar href=\"lib/freehep-export-2.1.1.jar\" />
    <jar href=\"lib/tclib.jar\" />
    <jar href=\"lib/i4jruntime.jar\" />
    <jar href=\"lib/FastInfoset.jar\" />
    <jar href=\"lib/jaxws-api.jar\" />
    <jar href=\"lib/jaxb-impl.jar\" />
    <jar href=\"lib/cytoscape-util-intr.jar\" />
    <jar href=\"lib/freehep-graphicsio-svg-2.1.1.jar\" />
    <jar href=\"lib/saaj-api.jar\" />
    <jar href=\"lib/jaxws-tools.jar\" />
    <jar href=\"lib/cytoscape-render-export.jar\" />
    <jar href=\"lib/freehep-graphicsio-2.1.1.jar\" />
    <jar href=\"lib/piccolo.jar\" />
    <jar href=\"lib/itext-2.0.4.jar\" />
    <jar href=\"lib/glf.jar\" />
    <jar href=\"lib/jdom-1.0.jar\" />
    <jar href=\"lib/freehep-io-2.0.2.jar\" />
    <jar href=\"lib/jaxws-rt.jar\" />
    <jar href=\"lib/freehep-jas-plotter-2.2.jar\" />
    <jar href=\"lib/freehep-swing-2.0.3.jar\" />
    <jar href=\"lib/jhall.jar\" />
    <jar href=\"lib/jaxb-api.jar\" />
    <jar href=\"lib/junit.jar\" />
    <jar href=\"lib/jsr173_1.0_api.jar\" />
    <jar href=\"lib/com-nerius-math-xform.jar\" />
    <jar href=\"lib/swing-layout-1.0.1.jar\" />
    <jar href=\"lib/giny.jar\" />
    <jar href=\"lib/undo.support.jar\" />
    <jar href=\"lib/commons-cli-1.x-cytoscape-custom.jar\" />
    <jar href=\"lib/freehep-graphicsio-ps-2.1.1.jar\" />
    <jar href=\"lib/sjsxp.jar\" />
    <jar href=\"lib/http.jar\" />
    <jar href=\"lib/cytoscape-cruft-obo.jar\" />
    <jar href=\"lib/cytoscape-render-stateful.jar\" />
    <jar href=\"lib/cytoscape-render-immed.jar\" />
    <jar href=\"lib/stax-ex.jar\" />
    <jar href=\"lib/ding.jar\" />
    <jar href=\"lib/cytoscape-geom-rtree.jar\" />
    <jar href=\"lib/jsr250-api.jar\" />
    <jar href=\"lib/swingx-2006_10_27.jar\" />
    <jar href=\"lib/cytoscape-graph-fixed.jar\" />
    <jar href=\"lib/saaj-impl.jar\" />
    <jar href=\"lib/cytoscape-task.jar\" />
    <jar href=\"lib/wizard.jar\" />
    <jar href=\"lib/fing.jar\" />
    <jar href=\"lib/coltginy.jar\" />
    <jar href=\"lib/freehep-xml-2.1.1.jar\" />
    <jar href=\"lib/l2fprod-common-all.jar\" />
    <jar href=\"lib/freehep-graphics2d-2.1.1.jar\" />
    <jar href=\"lib/freehep-graphicsio-java-2.1.1.jar\" />
    <jar href=\"lib/jnlp.jar\" />
    <jar href=\"lib/violinstrings-1.0.2.jar\" />
    <jar href=\"lib/concurrent.jar\" />
    <jar href=\"lib/streambuffer.jar\" />
    <jar href=\"lib/looks-2.1.4.jar\" />
    <jar href=\"lib/jsr181-api.jar\" />
    <!--These are the plugins you wish to load, edit as necessary.-->
    <jar href=\"plugins/biopax.jar\" />
    <jar href=\"plugins/cPath.jar\" />
    <jar href=\"plugins/yLayouts.jar\" />
    <jar href=\"plugins/browser.jar\" />
    <jar href=\"plugins/GraphMerge.jar\" />
    <jar href=\"plugins/CytoscapeEditor.jar\" />
    <jar href=\"plugins/psi_mi.jar\" />
    <jar href=\"plugins/filter.jar\" />
    <jar href=\"plugins/filters.jar\" />
    <jar href=\"plugins/ManualLayout.jar\" />
    <jar href=\"plugins/linkout.jar\" />
    <jar href=\"plugins/cpath2.jar\" />
    <jar href=\"plugins/SBMLReader.jar\" />
    <jar href=\"plugins/quick_find.jar\" />
    <jar href=\"plugins/TableImport.jar\" />
    <jar href=\"plugins/AutomaticLayout.jar\" />
  </resources>
  <!--This starts-up Cytoscape, specify your plugins to load, and other command line arguments.  Plugins not specified here will not be loaded.-->
  <application-desc main-class=\"cytoscape.CyMain\">
    <argument>-p</argument>
    <argument>org.mskcc.biopax_plugin.plugin.BioPaxPlugIn</argument>
    <argument>-p</argument>
    <argument>org.cytoscape.coreplugin.cpath.plugin.CPathPlugIn</argument>
    <argument>-p</argument>
    <argument>yfiles.YFilesLayoutPlugin</argument>
    <argument>-p</argument>
    <argument>browser.AttributeBrowserPlugin</argument>
    <argument>-p</argument>
    <argument>GraphMerge.GraphMerge</argument>
    <argument>-p</argument>
    <argument>cytoscape.editor.CytoscapeEditorPlugin</argument>
    <argument>-p</argument>
    <argument>org.cytoscape.coreplugin.psi_mi.plugin.PsiMiPlugIn</argument>
    <argument>-p</argument>
    <argument>filter.cytoscape.CsFilter</argument>
    <argument>-p</argument>
    <argument>cytoscape.filters.FilterPlugin</argument>
    <argument>-p</argument>
    <argument>ManualLayout.ManualLayoutPlugin</argument>
    <argument>-p</argument>
    <argument>linkout.LinkOutPlugin</argument>
    <argument>-p</argument>
    <argument>org.cytoscape.coreplugin.cpath2.plugin.CPathPlugIn2</argument>
    <argument>-p</argument>
    <argument>sbmlreader.SBMLReaderPlugin</argument>
    <argument>-p</argument>
    <argument>csplugins.quickfind.plugin.QuickFindPlugIn</argument>
    <argument>-p</argument>
    <argument>edu.ucsd.bioeng.coreplugin.tableImport.TableImportPlugin</argument>
    <argument>-p</argument>
    <argument>csplugins.layout.LayoutPlugin</argument>
    <argument>-N</argument>
    <argument>".$file."</argument>
  </application-desc>
</jnlp>";


    $c->res->content_type("application/x-java-jnlp-file");
    $c->res->body($doc);
    

}


1;
