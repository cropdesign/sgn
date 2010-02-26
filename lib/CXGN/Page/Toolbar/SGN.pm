package CXGN::Page::Toolbar::SGN;
use Moose;

has 'toolbar_data' => (
    is       => 'ro',
    isa      => 'ArrayRef',
    required => 1,
    default  => \&default_toolbar_data
);

sub default_toolbar_data {
    [
        {
            name => 'search',
            link => '/search/direct_search.pl?search=loci',
            desc => 'Search our database.',
            menu => [

                {
                    name => 'Genes',
                    link => '/search/direct_search.pl?search=loci',
                    desc => 'Gene search.'
                },
                {
                    name => 'Phenotypes',
                    link => '/search/direct_search.pl?search=phenotypes',
                    desc => 'Phenotype search.'
                },
                {
                    name => 'QTLs/Traits[beta]',
                    link => '/search/direct_search.pl?search=cvterm_name',
                    desc => 'A web interface for uploading QTL raw data, on-the-fly QTL mapping and search interface for QTLs.'
                },
                {
                    name => 'Unigenes',
                    link => '/search/direct_search.pl?search=unigene',
                    desc => 'Unigene search.'
                },
                {
                    name => "Families",
                    link => "/search/direct_search.pl?search=family",
                    desc => "Family search."
                },
                {
                    name => 'Markers',
                    link => '/search/direct_search.pl?search=markers',
                    desc => 'Marker search.'
                },
                {
                    name => 'BACs',
                    link => '/search/direct_search.pl?search=bacs',
                    desc => 'BAC (Bacterial Artificial Chromosome) search.'
                },
                {
                    name => 'ESTs and libraries',
                    link => '/search/direct_search.pl?search=est',
                    desc =>
'Find EST (Expressed Sequence Tag) libraries by keyword, e.g. library name, organism, tissue, development stage, or authors.'
                },
                {
                    name => 'Insitu database',
                    link => '/insitu/',
                    desc =>
'A database of in-situ images which can be updated by users, with image add, delete and annotation editing functions. The current images were generated by the floral genome project (FGP). For more information on the floral genome project, visit http://floralgenome.org/ and http://pgn.cornell.edu/.'
                },
                {
                    name => 'People',
                    link => '/search/direct_search.pl?search=directory',
                    desc => 'Search database of researchers who use SGN.'
                },
                {
                    name => 'Glossary',
                    link => '/search/glossarysearch.pl',
                    desc => 'Search a glossary of terms used on SGN.'
                },
                {
                    name => 'Images',
                    link => '/search/direct_search.pl?search=images',
                    desc => 'Search images contained in the SGN database.'
                }
            ]
        },
        {
            name => 'maps',
            link => '/cview/',
            desc => 'View and compare chromosomes from various organisms.',
            menu => [
                {
                    name => 'Tomato-EXPEN 2000',
                    link =>
'/cview/map.pl?map_id=9&amp;show_offsets=1&amp;show_ruler=1',
                    desc =>
'A combined genetic and physical map, showing the association of BACs (Bacterial Artificial Chromosomes) and BAC-Contigs to this genetic map.',
                    style => 'color:#400'
                },
                {
                    name => 'Tomato physical map',
                    link =>
'/cview/map.pl?map_id=p9&amp;show_offsets=1&amp;show_ruler=1',
                    desc =>
'This map is extracted from the S. lycopersicum LA925 x S. pennellii LA716 type F2.2000 genetic map.',
                    style => 'color:#400'

                },
                {
                    name => 'Tomato FPC contig map',
                    link =>
'/cview/map.pl?map_id=c9&amp;show_offsets=1&amp;show_ruler=1',
                    desc =>
'This map shows the positions of the Sanger HindIII/MboI FPC contigs relative to the the S. lycopersicum LA925 x S. pennellii LA716 type F2.2000 genetic map.',
                    style => 'color:#400'

                },

                {
                    name => 'Tomato IL map',
                    link =>
'/cview/map.pl?map_id=il6&amp;show_offsets=1&amp;show_ruler=1',
                    desc =>
'This map shows the locations of the introgressed fragments of S. pennellii in the S. lycopersicum genome, tied to the 1992 map.',
                    style => 'color:#400'
                },

                {
                    name => 'Tomato FISH map',
                    link => '/cview/map.pl?map_id=13',
                    desc =>
                      'Map of FISH (Fluorescent In Situ Hybridization) data.',
                    style => 'color:#400'
                },

                {
                    name => 'Tomato AGP map',
                    link => '/cview/map.pl?map_id=agp',
                    desc =>
'The Accessioned Golden Path (AGP) map for the tomato genome.',
                    style => 'color:#400'
                },
                {
                    name => 'Tomato ITAG map',
                    link => '/cview/map.pl?map_id=itag',
                    desc =>
'The ITAG map represents the position of the contigs used in the ITAG annotation.',
                    style => 'color:#400'
                },
                {
                    name  => 'Pepper COSII map',
                    link  => '/cview/map.pl?map_version_id=58',
                    desc  => 'COSII map',
                    style => 'color:#530'
                },

                {
                    name => 'Pepper-AC99',
                    link => '/cview/map.pl?map_id=11',
                    desc =>
'This map is based on 100 individuals from an inter-specific F2 population of Capsicum. annuum cv. NuMex RNaky and Capsicum chinense var PI159234. 426 molecular markers, including 359 SSR markers 3 specific PCR markers and 68 RFLP markers were used to construct this linkage map, with a total length of 1304.8 cM.',
                    style => 'color:#530'
                },
                {
                    name => 'Pepper-FA03',
                    link => '/cview/map.pl?map_id=10',
                    desc =>
'This map is based on 100 individuals from an inter-specific F2 population of Capsicum. annuum cv. NuMex RNaky and Capsicum frutescens BG 2814-6. 728 molecular markers, including 489 SSR markers 195 AFLP, 8 specific PCR markers and 36 RFLP were used to construct this linkage map, with a total length of 1358.7 cM.',
                    style => 'color:#530'
                },
                {
                    name => 'Potato-TXB 1992',
                    link => '/cview/map.pl?map_id=3',
                    desc =>
'This is a map based on a S. tuberosum x S. berthaultii BC S. tuberosum population reported in Tanksley et al (1992).',
                    style => 'color:#220'
                },
                {
                    name  => 'Eggplant COSII map',
                    link  => '/cview/map.pl?map_version_id=63',
                    desc  => 'COSII map',
                    style => 'color:#202'
                },
                {
                    name => 'Eggplant-LXM 2002',
                    link => '/cview/map.pl?map_id=6',
                    desc =>
'This map is based on 58 F2 plants from an interspecific cross between S. linnaeanum (MM195) and S. melongena (MM738) and contains 233 RFLP markers.',
                    style => 'color:#202'
                },
                {
                    name  => 'Tobacco SSR 2007',
                    link  => '/cview/map.pl?map_id=15',
                    desc  => 'Nicotiana tabacum',
                    style => 'color:#333'
                },
                {
                    name => 'Arabidopsis COSII',
                    link =>
'/cview/map.pl?map_id=7&amp;show_offsets=1&amp;show_ruler=1',
                    desc  => 'Arabidopsis thaliana sequenced-based COSII map.',
                    style => 'color:#002800'
                },

                {
                    name => 'Comparative viewer...',
                    link =>
'/cview/view_chromosome.pl?show_offsets=1&amp;show_ruler=1&amp;show_offsets=1',
                    desc =>
                      'Compare chromosomes from different maps, side by side.',
                    style => 'font-style:italic;'
                },
                {
                    name  => 'See all available maps...',
                    link  => '/cview/index.pl',
                    desc  => 'All chromosome maps that are available on SGN.',
                    style => 'font-style:italic'
                }
            ]
        },
        {
            name => 'sequencing',
            link => '/about/tomato_sequencing/',
            desc => 'Read about sequencing projects.',
            menu => [
                {
                    name => 'International Tomato Sequencing Project',
                    link => '/about/tomato_sequencing.pl',
                    desc =>
'The gene-rich euchromatic portion of the tomato genome is being sequenced by an international consortium. This page gives more information about the progress of the sequencing effort.'
                },
                {
                    name => 'Solanaceae project (SOL)',
                    link => '/solanaceae-project/index.pl',
                    desc =>
'Over the coming decade the International Solanaceae Genome Project (SOL) will create a coordianted network of knowledge about the Solanaceae family aimed at answering two of the most important questions about life and agriculture: How can a common set of genes/proteins give rise to such a wide range of morphologically and ecologically distinct organisms that occupy our planet? How can a deeper understanding of the genetic basis of diversity be harnessed to better meet the needs of society in an environmentally-friendly way? On this page, you will find more information about the strategy that will be used to answer these questions.'
                },
                {
                    name => 'U.S. tomato sequencing project',
                    link => '/about/us_tomato_sequencing.pl',
                    desc => 'An overview of the US Tomato Sequencing Project'
                },
                {
                    name => 'Browse Tomato Genome Data',
                    link => '/gbrowse/',
                    desc => 'Browse tomato genomic sequences and annotations.',
                },
            ]
        },
        {
            name => 'tools',
            link => '/tools/',
            desc => 'Tools for biologists.',
            menu => [

                #Sequence Analysis
                { name => 'Sequence Analysis' },
                {
                    name => 'BLAST',
                    link => '/tools/blast/',
                    desc =>
'Compare a given sequence to known sequences using the Basic Local Alignment Search Tool, NCBI BLAST v. 2.2.9 [May-01-2004]'
                },
                {
                    name => 'Alignment Analyzer',
                    link => '/tools/align_viewer/',
                    desc =>
'Calculate alignments using Muscle v3.6, or view existing alignments using our navigator.'
                },
                {
                    name => 'Tree Browser',
                    link => '/tools/tree_browser/',
                    desc =>
'Visualize and manipulate phylogenetic trees from a newick-formatted string.'
                },
                {
                    name => 'Intron Finder',
                    link => '/tools/intron_detection/find_introns.pl',
                    desc =>
'The SGN Intron Finder works by doing a BLAST search for <em>Arabidopsis thaliana</em> proteins that are similar to the translated protein sequence of the DNA input.'
                },

                #Mapping
                { name => 'Mapping' },
                {
                    name => 'Comparative Viewer',
                    link =>
'/cview/view_chromosome.pl?show_offsets=1&amp;show_ruler=1&amp;show_offsets=1',
                    desc =>
                      'Compare chromosomes from different maps, side by side.',
                },
                {
                    name => 'Fast Mapping',
                    link => '/tools/fastmapping/',
                    desc =>
'A quick mapping program which requires a file in mapmaker format that you can upload.'
                },
                {
                    name => 'CAPS Designer',
                    link => '/tools/caps_designer/caps_input.pl',
                    desc =>
'Designs CAPS (Cleaved Amplified Polymorphic Sequence) assays for up to twelve sequences. Two types of nucleotide inputs are accepted: fasta sequences and clustal aligment. It generates a list of polymorphic enzymes that cut the sequences into different length products.'
                },
                {
                    name => 'Seed BAC Finder',
                    link => '/tools/seedbac/',
                    desc =>
'Lists all anchored BACs for a given chromosome to help identify seed BACs, or suggests a seed BAC given a marker name. '
                },
		{
                    name => 'QTL Mapping',
                    link => '/qtl/',
                    desc =>
'A web interface for uploading QTL raw data, on-the-fly QTL mapping and search interface for QTLs.'
                },
                #Molecular Biology
                { name => "Molecular Biology" },

# 				{
#                     name=>'Virtual DNA Digestor',
#                     link=>'/tools/dna_digestion/digestion.pl',
#                     desc=>'Visualize the digestion sites in a given sequence for a selected enzyme.'
#                 },
#				{
#                    name=>'Primer 3',
#                    link=>'/tools/primer3/',
#                    desc=>'Design primers for PCR reactions using several criteria.'
#                },
                {
                    name => 'Signal Peptide Finder',
                    link => '/tools/sigpep_finder/input.pl',
                    desc =>
'This method of searching for signal sequences is designed to complement SignalP, and has similar success rates.'
                },
                {
                    name => 'In Silico PCR',
                    link => '/tools/insilicopcr',
                    desc => 'In Silico PCR tool based on BLAST'
                },

                #Systems Biology
                { name => "Systems Biology" },
                {
                    name => 'SolCyc Biochemical Pathways',
                    link => '/tools/solcyc/',
                    desc => 'An interactive map of metabolic pathways.'
                },

                #Bulk Query
                { name => 'Bulk Query' },

                {
                    name => 'Unigene and BAC information',
                    link => '/bulk/input.pl',
                    desc =>
'Rather than sifting through SGN page by page, here you can download large quantities of information in a single step.'
                },
                {
                    name => 'FTP Site',
                    link => '/bulk/input.pl?mode=ftp',
                    desc =>
'FTP (File Transfer Protocol) is more efficient than HTTP (HyperText Transfer Protocol) for transferring large files. Many of our larger archive files are accessible here.'
                },
                {
                    name => 'ID Converter (SGN <=> TIGR)',
                    link => '/tools/convert/input.pl',
                    desc =>
'The Institute for Genomic Research and SGN maintain independent unigene databases, entries in which tend to have common member ESTs (Expressed Sequence Tags), although they tend not to correspond completely. This tool uses common members to convert back and forth between the two identifier sets.'
                },

                #Other
                { name => "Other" },
                {
                    name => 'SGN Ontology Browser [beta]',
                    link => '/tools/onto/',
                    desc => 'AJAX ontology browser.'
                },
            ]
        }
    ];
}

=head2 headings

    my @headings=$tb->headings();  

Returns a list of the names of the main toolbar buttons
('search','maps','sequencing','tools'). These can then be used to
access other parts of this object.

=cut

sub headings {
    map { $_->{name} } @{ shift->toolbar_data };
}

=head2 heading_link 

    print '<a href="'.$_->heading_link($heading).qq|">$_</a>|
        for @headings;

Takes a heading (which can be gotten from the headings
function). Returns the stuff that is supposed to go in the href
attribute for this heading.

=cut

sub heading_link {
    my $self = shift;
    my ($heading) = @_;
    for my $heading_hash ( @{ $self->toolbar_data } ) {
        if ( $heading_hash->{name} eq $heading ) {
            return $heading_hash->{link};
        }
    }
}

=head2 heading_desc

    for my $heading(@headings)
    {
        print"$heading description: ".$tb->heading_desc($heading);
    }

Takes a heading (which can be gotten from the headings function). Returns a description of what the linked-to page is supposed to do.

=cut

sub heading_desc {
    my $self = shift;
    my ($heading) = @_;
    for my $heading_hash ( @{ $self->toolbar_data } ) {
        if ( $heading_hash->{name} eq $heading ) {
            return $heading_hash->{desc};
        }
    }
}

=head2 menu_options

    for my $heading(@headings)
    {
        my @menu_options=$tb->menu_options($heading);  
    }

Takes a heading (which can be gotten from the headings function). Returns a list of the names of the options below that heading ('Genes','Unigenes','Markers','BACs','Annotations','ESTs and libraries','People'). Like the list returned by the headings function, these can then be used to access other information from this object.

=cut

sub menu_options {
    my $self = shift;
    my ($heading) = @_;
    for my $heading_hash ( @{ $self->toolbar_data } ) {
        if ( $heading_hash->{name} eq $heading ) {

#            the toolbar links don't work if you add tooltips to the toolbar. fix?
#            my @menu_options=map {CXGN::Page::FormattingHelpers::tooltipped_text($_->{name},$_->{desc},'toolbar_help')} @{$heading_hash->{menu}};
            my @menu_options = map { $_->{name} } @{ $heading_hash->{menu} };
            return @menu_options;
        }
    }
}

=head2 option_link

    for my $heading(@headings)
    {
        my @menu_options=$tb->menu_options($heading);
        for my $option(@menu_options)
        {
            print"<a href=\"".$tb->option_link($heading,$option)."\">$option</a>";            
        } 
    }

Takes a heading and an option. Returns the stuff that goes in the href for that option.

=cut

sub option_link {
    my $self = shift;
    my ( $heading, $option ) = @_;

    #print STDERR "option_link received: heading '$heading' option '$option'\n";
    for my $heading_hash ( @{ $self->toolbar_data } ) {
        if ( $heading_hash->{name} eq $heading ) {
            for my $menu_option_hash ( @{ $heading_hash->{menu} } ) {
                if ( $menu_option_hash->{name} eq $option ) {
                    my $option_link = $menu_option_hash->{link};
                    if ($option_link) { return $option_link }
                }
            }
        }
    }

    #die"Option link not found for heading '$heading' option '$option'";
}

=head2 option_style

    for my $heading(@headings)
    {
        my @menu_options=$tb->menu_options($heading);
        for my $option(@menu_options)
        {
            my $style = $tb->option_style($heading,$option);            		print "<span style='$style'>";
				#print option item
			  print "</span>";
        } 
    }

Takes a heading and an option. Returns the stuff that goes in the href for that option. Override me.


=cut

sub option_style {
    my $self = shift;
    my ( $heading, $option ) = @_;

    #print STDERR "option_link received: heading '$heading' option '$option'\n";
    for my $heading_hash ( @{ $self->toolbar_data } ) {
        if ( $heading_hash->{name} eq $heading ) {
            for my $menu_option_hash ( @{ $heading_hash->{menu} } ) {
                if ( $menu_option_hash->{name} eq $option ) {
                    my $option_style = $menu_option_hash->{style};
                    if ($option_style) { return $option_style }
                }
            }
        }
    }

    #die"Option link not found for heading '$heading' option '$option'";
}

=head2 option_desc

    for my $heading(@headings)
    {
        my @menu_options=$tb->menu_options($heading);
        for my $option(@menu_options)
        {
            print"$option description: ".$tb->option_desc($heading,$option);            
        } 
    }

Takes a heading and an option. Returns the stuff that goes in the href for that option.

=cut

sub option_desc {
    my $self = shift;
    my ( $heading, $option ) = @_;
    for my $heading_hash ( @{ $self->toolbar_data } ) {
        if ( $heading_hash->{name} eq $heading ) {
            for my $menu_option_hash ( @{ $heading_hash->{menu} } ) {
                if ( $menu_option_hash->{name} eq $option ) {
                    return $menu_option_hash->{desc};
                }
            }
        }
    }
}

sub as_html {
    my $self     = shift;
    my @headings = $self->headings();

    #decorate the headings array with its menu alignments
    @headings = map [ 'center', $_ ], @headings;
    if ( @headings > 1 ) {
        $headings[0]->[0] = 'left';

        #      $headings[-1]->[0] = 'right';
    }

    my $index        = 0;
    my $toolbar_html = <<EOH
<table class="toolbar" summary="" cellspacing="0" cellpadding="0">
<tr class="toolbar_content">
  <td class="toolbar_l"><img src="/documents/img/toolbar_l.gif" alt="" /></td>
EOH
      . join(
        '',
        map {
            my ( $align, $heading ) = @$_;
            my $heading_link = $self->heading_link($heading);
            my $popup_menu   = $self->render_menu($heading);
            <<EOH
  <td class="toolbar_content" onmouseover="CXGN.Page.Toolbar.showmenu('$heading')" onmouseout="CXGN.Page.Toolbar.hidemenu()" style="text-align: $align">
    <a class="toolbar_menuname" href="$heading_link">$heading</a>
    <div style="position: relative; text-align: left">
$popup_menu
    </div>
  </td>
EOH
          } @headings
      )
      . <<EOH

  <td class="toolbar_search toolbar_content">
    <form class="quicksearch" action="/search/quick_search.pl">
      <input class="quicksearch field" type="text" name="term" size="14" />
      <input class="quicksearch imgbutton" type="image" src="/documents/img/sol_search_button.gif" value="sol search" />
    </form>
  </td>
  <td class="toolbar_r" height="1">
    <img src="/documents/img/toolbar_r.gif" alt="" />
  </td>
</tr>
</table>
<script language="JavaScript" type="text/javascript">
EOH
      . join( '', map <<EOH, @headings )
  CXGN.Page.Toolbar.addmenu('$_->[1]');
EOH
      . <<EOH
</script>
EOH
}

#return HTML for the popup menu for the  menu heading
sub render_menu {
    my ( $self, $heading ) = @_;

    my @menu_options = $self->menu_options($heading);

    return <<EOH
    <table cellspacing="0" cellpadding="0" id="$heading" class="toolbar_popmenu">
EOH
      . join(
        '',
        map {    ########## FINISH REWRITING THIS
            my $menu_option = $_;

     #menu options for each heading, these will be invisible unless hovered over
            my $option_link = $self->option_link( $heading, $menu_option );
            $option_link &&=
              qq|<a class="toolbar_item" href="$option_link">$menu_option</a>|;
            $option_link ||= qq|<span class="toolbar_item">$menu_option</span>|;
            my $style = $self->option_style( $heading, $menu_option );
            <<EOH
       <tr><td class="toolbar_item">$option_link</td></tr>
EOH
          } @menu_options
      )
      . <<EOH
    </table>
EOH
}

sub index_page {
    my $self         = shift;
    my ($heading)    = @_;
    my $content      = '';
    my $disp_heading = ucfirst($heading);
    for my $option ( $self->menu_options($heading) ) {
        my $link = $self->option_link( $heading, $option );
        my $desc = $self->option_desc( $heading, $option );
        if ($link) {
            $content .=
qq|<div style="margin-left:10px"><a href="$link">$option</a> - $desc<br /><br /></div>|;
        }
        else {
            $content .=
"<br /><span style='font-size:1.00em; font-weight:bold'>$option</span> $desc";
            $content .= ": $desc" if $desc;
            $content .= "<br /><br />";
        }
    }
    return <<END_HEREDOC;
<center>
<table class="boxbgcolor2" width="100%" summary="">
<tr>
<td width="15%">&nbsp;</td>
<td width="70%" class="left">
<div class="boxcontent">
<div class="subheading" style='font-size:1.1em'>
$disp_heading<br />
</div>
<div class="boxsubcontent">
$content
</div>
</div>
</td>
<td width="15%">&nbsp;</td>
</tr>
</table>
</center>
END_HEREDOC
}

__PACKAGE__->meta->make_immutable;

###
1;    #
###
