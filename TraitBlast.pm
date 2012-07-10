
=head1 NAME

TraitBlast.pm - a Perl script to determine genes that are missing from contrast groups

=head1 DESCRIPTION

Synopsis:

mx-run TraitBlast --blast_file <blast-m8-file>

Options:

--blast_file: A blast file in m8 format

Output: A list of gene identifiers

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut


package TraitBlast;

use Moose;
use Config::Any;

##my $conf = Config::Any->read_files(files=> [ 'traitblast.conf' ]);



with 'MooseX::Runnable';
with 'MooseX::Getopt';

has 'blast_file' => (
    is => 'rw',
    isa => 'Str',
    traits => ['Getopt'],
    );


has 'outgroup' => (is=> 'rw',
		   isa => 'ArrayRef', 
    );


has 'ingroup' => (is =>'rw',
		  isa => 'ArrayRef',
    );

has 'contrast_group' => (is=> 'rw',
			 isa => 'ArrayRef',
    );

has 'min_ingroups' => (is=>'rw',
		       isa=> 'Int',
		       default => 7,
    );

has 'min_outgroups' => (is=>'rw',
			isa =>'Int',
			default => 5,
    );


sub is_in_contrast_group { 
    my $self = shift;
    my $species = shift;    
    return $self->is_in_list($species, @{$self->contrast_group()});
}

sub is_in_outgroup { 
    my $self = shift;
    my $species = shift;
    return $self->is_in_list($species, @{$self->outgroup()} );
}

sub is_in_ingroup { 
    my $self = shift;
    my $species = shift;
    return $self->is_in_list($species, @{$self->ingroup()} );
}

sub is_in_list { 
    my $self = shift;
    my $element = shift;
    my @list = @_;
    foreach my $e (@list) { 
	if ($e eq $element) {  
	    return 1;
	    
	}
    }   
    return 0;
}


sub run { 
    my $self = shift;

    $self->outgroup( [ qw | rice maize brachypodium sorghum selaginella | ]);
    $self->ingroup( [ qw | tomato castorbean soybean medicago poplar cassava grape potato papaya cucumber | ]);
    $self->contrast_group ([ qw | arabidopsis lyrata brapa crubella thalophila |]);
    
    my $old_q = "";

    my $outgroup_flag = 0;
    my $contrast_flag =0;
    my $ingroup_flag =0;
    my $candidate_flag = 0;

    my %outgroups = ();
    my %ingroups  = ();

    my %species = ();

    my $min_score = 10000;
    my $max_score = 0;

    my $old_min_score = 0;
    my $old_max_score = 0;

    my $candidate_count = 0;
    my $skip =0;

    my $query_species = "";

    open (my $F, "<", $self->blast_file()) || die "Can't open \"".$self->blast_file()."\". ";
    while (<$F>) { 
	chomp;
	my ($q, $s, $score) = ( split /\t/ ) [0,1,11];

	if ($score > $max_score) { $max_score = $score; }

	if ($score < $min_score) { $min_score = $score; }
	
	
	if ($q ne $old_q) { 

	    if ($candidate_flag) { 
		print "$old_q\n";# $min_score $max_score\n";
		$candidate_count++;
		
	    }
	    else { 
		#print "NOT A CANDIDATE: $old_q\n";
	    }

	    $min_score = 10000;
	    $max_score = 0;

	    $candidate_flag = 0;
	    $outgroup_flag = 0;
	    $contrast_flag = 0;

	    %ingroups = ();
	    %outgroups = ();

	    $skip = 0;

	    $query_species = "";
	}

	if ($skip) {  next; }

	#print STDERR "$q, $s, $score\n";
	my $query_species = $self->id2species($q);

	my $species = $self->id2species($s);

	$species{$species}++;

	if ($self->is_in_contrast_group($species)) { 
	    $contrast_flag =1;
	}

	if ($self->is_in_outgroup($species)) { 
	    #warn "Adding $species to outgroup... ".scalar(keys(%outgroups))."\n";
	    $outgroups{$species}++;
	    $outgroup_flag = 1;
	}
	
	if ($self->is_in_ingroup($species)) { 
	    #warn "Adding $species to ingroup... ".scalar(keys(%ingroups))."\n";
	    $ingroups{$species}++;
	    $ingroup_flag = 1;

	    # don't allow outgroup matches before ingroup matches
	    #
	    if (keys(%ingroups) < $self->min_ingroups() && keys(%outgroups)>0) { 
		$candidate_flag=0;
		$skip = 1;
		next;
		
	    }
	}

	if ($query_species eq $species && ($q ne $s) ) { 
	    print "skipping $q because $q and $s are both in this group and of species $species.\n";
	    next;
	}

	if ((scalar(keys(%outgroups))>=$self->min_outgroups()) && (scalar(keys(%ingroups))>$self->min_ingroups()) && ($contrast_flag ==0) ) { 
	    $candidate_flag= 1;
	}

	$old_q = $q;
	$old_min_score = $min_score;
	$old_max_score = $max_score;
    }
    

    print STDERR "Species detected in dataset: ".(join ", ", (map { "$_ ($species{$_})" } sort(keys(%species))))."\n";

    print STDERR "Candidates detected: $candidate_count\n";
    
}


sub id2species { 
    my $self = shift;
    my $id = shift;

    if ($id =~ /At[1-5CM]g\d+/i) { return "arabidopsis"; }
    if ($id =~ /^Solyc/)   { return "tomato"; }
    if ($id =~ /^GSVIV/)   { return "grape"; }
    if ($id =~ /^GRMZM/ || $id=~ /^\w{2}.*\d+\_FGP\d+/)   { return "maize"; }
    if ($id =~ /evm.model.supercontig/i || /evm.TU.contig/i)  { return "papaya"; }
    if ($id =~ /^LOC_OS/i || $id=~ /^ChrSy.fgenesh.mRNA/ || $id=~ /^ChrUn.fgenesh.mRNA/) { return "rice"; }
    if ($id =~ /^Sb/)      { return "sorghum"; }
    if ($id =~ /Glyma/i)   { return "soybean"; }
    if ($id =~ /Bradi/)    { return "brachypodium"; }
    if ($id =~ /\d{5}\_?\.?m\d{5,6}/) { return "castorbean"; }
    if ($id =~ /^Cucsa/) { return "cucumber"; }
    if ($id =~ /^Sm/)   { return "selaginella"; }
    if ($id =~ /Aly/i)   { return "lyrata"; }
    if ($id =~ /POPTR_/) { return "poplar"; }
    if ($id =~ /^Bra\d{6}/) { return "brapa";}
    if ($id =~ /^Medtr[1-8]g\d{6}/ || $id=~ /^AC/ || $id=~ /^CU/ || $id=~ /^contig_\d+/) { return "medicago";}   
    if ($id =~ /^CM/ || $id=~ /^Lj/ || $id =~ /^chr[1-6]/) { return "lotus";}
    if ($id =~ /^orange/) { return "orange";}
    if ($id =~ /^clementine/) { return "clementine";}
    if ($id =~ /^cassava/) { return "cassava";}
    if ($id =~ /^mgv/) { return "mimulus";}
    if ($id =~ /^PGSC/) { return "potato";}
    if ($id =~ /^Carubv/) { return "crubella";}
    if ($id =~ /^Thhalv100/) { return "thalophila";}
    if ($id =~ /^Pp1s/) { return "moss";}
    if ($id =~ /^Cre/) { return "chlamy";}
    if ($id =~ /^Vocar200/) { return "volvox";}
    if ($id =~ /^Ostta4-/) { return "tauri";}
    if ($id =~ /^Olu-/) { return "lucimarinus";}





    #die "don't know $id\n";
    
    return "unknown";


}



1;
