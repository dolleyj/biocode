#!/usr/bin/env perl

=head1 NAME

fragment_sequence.pl - fragment FASTA sequences according to user's parameters.

Based on Enhancement request:
https://github.com/jorvis/biocode/issues/36


=head1 SYNOPSIS

USAGE: fragment_sequence.pl 
		--input_file=/path/to/sequence_file.fsa 
		--output_file=/path/to/ouput_file.fsa
		--min_overlap_distance=10
		--max_overlap_distance=100
		[ --fragmentation_factor=5 
		--log
		]

=head1 OPTIONS

B<--input_file,-i>
	name the path to the input FASTA file

B<--output_file,-o>
	the path and name of the output FASTA file

B<--min_overlap_distance,-m>
	This is the minimal distance a fragment will overlap.
	* Note: a (-) negative number will cause fragments to be at most N part
	
	min_overlap_distance=50		means the fragment will overlap at least 50 bps
	min_overlap_distance=-200	means the fragment will at most be 200 bps apart

B<--max_overlap_distance,-n>
	This is the maximum distance a fragment will overlap. 
	* Note: MUST be positive(+)

B<--fragmentation_factor,-f>
	optional.  This is the degree to which fragmentation will occur.
	A higher fragmentation factor yields more (and shorter) fragments ( default = 5 ).
	* Note: MUST be positive(+) whole number 

B<--log,-l> 
	Log file

B<--help,-h>
	This help message

=head1  DESCRIPTION

The purpose of this script is to break apart a known nucleotide sequence, 
or sequences, into smaller fragments. A fragment can either contain a(n): 

1) Overlap with the sequence before it 
	

	5'---------------------3'
             	5'----------------------------------3'
               
               
               
2) Gap between the sequence before it
	

	5'---------------------3'
                                          5'----------------------------------3'   
                                     
The length of an overlap is noted by a positive(+) min_overlap_distance
and max_overlap_distance.
The length of a gap is noted by a negative(-) min_overlap_distance. 
	*Note* that max_overlap_distance is still positive(+)
			...in fact, it should ALWAYS be positive.

fragmentation_factor is the degree to which fragments will be made. Essentially,
the sequence length is divided by the fragmentation factor. This should be 
set to at least 2.

Example:
Sequence length of 1000 bases
Fragmentation Factor of 10
				
Produces fragments that are approximately 100 bases long
				
				
This script is an attempt to fullfill the feature request:
	https://github.com/jorvis/biocode/issues/36


=head1  INPUT

One FASTA file enters.
Can be a single sequence or multiple sequences.

=head1  OUTPUT

One FASTA file leaves. 
Sequence(s) are fragmented according to user's parameters

=head1  CONTACT

	Dustin Olley

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
						'input_file|i=s',
						'output_file|o=s',
						'min_overlap_distance|m=s',
						'max_overlap_distance|n=s',
						'fragmentation_factor|f=s',
						'log|l=s',
						'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
	open($logfh, ">$options{log}") || die "can't create log file: $!";
}

my $output_fh;
if (defined $options{output_file}) {
	open($output_fh, ">$options{output_file}") || die "can't create output file: $!";
}

#my $min_overlap;
#if (defined $options{min_overlap_distance}) {
#	$min_overlap = $options{min_overlap_distance};
#}
#
#my $max_overlap;
#if (defined $options{max_overlap_distance}) {
#	$min_overlap = $options{max_overlap_distance};
#}
#
#my $frag_factor;
#if (defined $options{fragmentation_factor}) {
#	$min_overlap = $options{fragmentation_factor};
#}



open (my $inputfh, "<$options{input_file}") || die "can't open input file: $!";

my $current_header;
my @current_seq_lines = ();


while (my $line = <$inputfh>) {
	chomp $line;
	#print "$line\n";
	if ($line =~ /^\>(.+)$/) {
		my $next_header = $1;
		
#		$current_header = $next_header;
#		@current_seq_lines = ();
		
		if ($current_header) {
			process_seq($current_header, \@current_seq_lines);
		}
		$current_header = $next_header;
		@current_seq_lines = ();
		
	} else {
		#push(my @current_sequence, $line);
		push(@current_seq_lines, $line);  
	}
#	if ($current_header) {
#		process_seq($current_header, \@current_seq_lines);
#		print "current_header defined\n";
#	}
}#end while loop

#for single sequence OR last sequence in a multi sequence file
process_seq($current_header, \@current_seq_lines); 

close $output_fh;

exit(0);

sub process_seq {
	my ($header, $lines) = @_;

	#Do overlap & fragmentation arguments go here or above?
	#...they seem to be working here
	my $min_overlap;
	if (defined $options{min_overlap_distance}) {
		$min_overlap = $options{min_overlap_distance};
	}
	my $max_overlap;
	if (defined $options{max_overlap_distance}) {
		$max_overlap = $options{max_overlap_distance};
	}
	my $frag_factor;
	if (defined $options{fragmentation_factor}) {
		$frag_factor = $options{fragmentation_factor};
	}
	

	# join seq lines into 1 string
	my $seq = join('', @$lines);
	$seq =~ s|\s||g;
	
	my $seq_length = length($seq);
	#print "seq_length = $seq_length\n"; #FOR TESTING
	my $min_ol_number = abs($min_overlap);
	
	#print "header: $header\nseq: $seq\n"; #FOR TESTING


	if ($seq_length < $min_ol_number ) {
	#if ($seq_length < abs($min_overlap)) {
		_log(">$header\nSkipped because sequence is shorter than min_overlap_distance\n");
		print $output_fh "$header was skipped because it is shorter than min_overlap_distance\n";
		#print "$seq\n"; # for testing
	} elsif ($seq_length < $max_overlap) {
			_log(">$header\nSkipped because sequence is shorter than max_overlap_distance.\n");
			print "$header was skipped because it is shorter than max_overlap_distance\n";
			#print "$seq\n"; # for testing	
	} else {
		
		my @seq_range = (0 .. ($seq_length -1)); #setup a list that's the length of the sequence
		my $counter;
		my $frag_length = int($seq_length / $frag_factor); ### ASSUMES a higher fragmentation factor means smaller fragments 
		
		for (my $base_start = 0; $base_start < @seq_range; $base_start += $frag_length) {
		#foreach my $base_start(@seq_range) { #makes too many fragments
#			$counter++;

			my $base_location = $base_start + $frag_length + 1; #which number base are we at?
			if ($base_location >= $seq_length) {
##				print "Done!\n";
##				return;
				next;
			} else {
			
				my @overlap_range = ($min_overlap .. $max_overlap); #fragments are a distance apart --> overlapping...
				#print @overlap_range;
				my $ol_range_length = @overlap_range;
				
				my $overlap_increment = int($ol_range_length / $frag_factor);
				for (my $ol_number = $min_overlap; $ol_number < @overlap_range; $ol_number += $overlap_increment) {
				#foreach my $ol_number(@overlap_range) {	#makes too many fragments
					#$counter++;
					
					# if min_overlap is negative, fragments will be a distance apart... 
					if ( $ol_number <= 0 ) {
						$counter++;
#						my $seq_fragment = substr($seq, $base_start, $frag_length);	# get a sequence fragment
#						#push @seq_fragments, $seq_fragment;	
#						
#						print $output_fh ">$header, frag: $counter, overlap_distance: 0\n";
#						while ( $seq_fragment =~ /(.{1,60})/g ) {
#							print $output_fh "$1\n";
#						}
				    
				    #gives a start point to the fragment (distance apart)
						my $overlap_start = $base_start + $frag_length + abs($ol_number); 
						if ($overlap_start >= $seq_length) {
##							print "Done!\n";
##							return;
							next;
						} else { 
							my $overlapping_seq_fragment = substr($seq, $overlap_start, $frag_length);	# get a overlapping sequence
							#push @seq_fragments, $overlapping_seq_fragment;	
						
							print $output_fh ">$header, fragment: $counter, starts_at_base: $base_location, length: $frag_length, overlap_distance: $ol_number, fragment_factor: $frag_factor\n";
							while ( $overlapping_seq_fragment =~ /(.{1,60})/g ) {
								print $output_fh "$1\n";
							}
						}
					} 
					# if min_overlap is positive. fragment will be overlapped...
					elsif ( $ol_number > 0 ) {
						$counter++;
#						my $seq_fragment = substr($seq, $base_start, $frag_length);	# get a sequence fragment
#						#push @seq_fragments, $seq_fragment;	
#						
#						print $output_fh ">$header, frag: $counter, overlap_distance: 0\n";
#						while ( $seq_fragment =~ /(.{1,60})/g ) {
#							print $output_fh "$1\n";
#						}
			
						my $overlap_start = $base_start + $frag_length - $ol_number; #gives a start point to the overlapping fragment 
						my $overlapping_seq_fragment = substr($seq, $overlap_start, $frag_length);	# get a overlapping sequence
						#push @seq_fragments, $overlapping_seq_fragment;	
						
						#my $header_new = ">$header, frag: $counter, overlap_distance: $ol_number\n";
						
						print $output_fh ">$header, fragment: $counter, starts_at_base: $base_location, length: $frag_length, overlap_distance: $ol_number, fragment_factor: $frag_factor\n";
						while ( $overlapping_seq_fragment =~ /(.{1,60})/g ) {
							print $output_fh "$1\n";
						}		
					}
				} #end for ol_number loop 
			} #end else 
		} #end for base_start loop
	} #end else
} #end sub
	



sub _log {
	my $msg = shift;

	print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
	my $options = shift;
    
	## make sure required arguments were passed
	my @required = qw( input_file output_file min_overlap_distance max_overlap_distance );
	for my $option ( @required ) {
		unless  ( defined $$options{$option} ) {
			die "--$option is a required option";
		}
	}

	##
	## you can do other things here, such as checking that files exist, etc.
	## check input is FASTA - .fasta .fas .fa .fsa .fna .ffn .faa .frn
	if ( ! $$options{input_file} =~ m/\.(fasta|fas|fa|fsa|fna|ffn|faa|frn)/ ) {
		die "--input_file not a FASTA file."; 
 	}
 	## no negative fragmentation factor allowed
 	if ( $$options{fragmentation_factor} < 0 ) {
 		die "--fragmentation_factor cannot be negative number.";
 	}
 	
 	## fragmentation of 2 is min. cuts sequences in half
 	if ( $$options{fragmentation_factor} >= 0 && $$options{fragmentation_factor} < 2 ) {
 		die "If you want fragments, --fragmentation_factor should be at least 2.";	
 	}

	## handle some defaults
	$options{fragmentation_factor}   = '5'  unless ($options{fragmentation_factor}); #default frag_factor to 5 
}
