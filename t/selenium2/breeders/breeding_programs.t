

use lib 't/lib';

use Test::More;

use SGN::Test::WWW::WebDriver;


my $t = SGN::Test::WWW::WebDriver->new();

$t->while_logged_in_as("curator", sub { 
    $t->get_ok('/breeders/manage_programs');
    
    my $new_bp_link = $t->find_element_ok('new_breeding_program_link', 'id', 'new breeding program link');

    $new_bp_link->click();

    sleep(2);

    my $breeding_program_name_input = $t->find_element_ok('new_breeding_program_name', 'id', 'find add breeding program name input');

    $breeding_program_name_input->send_keys('WEBTEST');

    my $breeding_program_desc_input = $t->find_element_ok('new_breeding_program_desc', 'id', 'find add breeding program description input');

    $breeding_program_desc_input->send_keys('Test description.');

    my $ok_button = $t->find_element_ok('new_breeding_program_submit', 'id', 'find add breeding program button');


    $ok_button->click();

    print STDERR "\n\nCLICKED OK... so far so good...\n\n";

    $t->driver()->accept_alert();

    sleep(2); # wait until page is re-loaded

    $t->get_ok('/breeders/manage_programs');
    sleep(2);

    ok($t->driver->get_page_source() =~ m/WEBTEST/, "breeding program addition successful");

#    print STDERR $t->driver->get_page_source();

    my $delete_link = $t->find_element_ok('delete_breeding_program_link_WEBTEST', 'id', 'find breeding program delete link');

    sleep(2); 

     $delete_link->click();

    sleep(2);
    
    $t->driver()->accept_alert();
    
    $t->get_ok('/breeders/manage_programs');
    
    
    print STDERR "Marker 4\n";

#    sleep(2);

#  ok($t->driver->get_page_source() !~ m/WEBTEST/, "breeding program deletion successful");

    done_testing();
		       
		       });
