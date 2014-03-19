#!/gsc/bin/perl
#
#   THIS SCRIPT IS A WRAPPER AROUND THE BasoVac.pm SOMATIC VARIANT CALLER
#
#   Copyright (C) 2010 Washington University
#
#   mwendl@wustl.edu
#
#   This script is free software; you can use it for free, even
#   redistribute it under the terms of the GNU Library General Public License as
#   published by the Free Software Foundation.
#
#
#   PROGRAMMER NOTES
#   ------------------------------
#
#   (1) typical command line:
#
#
##      DNMT3A Mutation in AML 31
#
#   ./basovac.pl --rfile /gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa --tfile /gscmnt/sata423/info/medseq/dlarson/AML31_tumor_DNMT3A.bam --nfile /gscmnt/sata423/info/medseq/dlarson/AML31_normal_DNMT3A.bam --normal_contam_by_tumor 0.05 --tumor_contam_by_normal 0.03 --tumor_background_mutation_rate 0.000001 --normal_variant_freq 0.5 --tumor_variant_freq 0.5
#
#   (2) description of samtools output
#
#       ordered by chromosome (low->hi) then by position (low->hi)
#
#       column 0: chromosome
#       column 1: base position
#       column 2: reference base
#       column 3: total number of covering reads
#       column 4: symbolic string-based read summary
#       column 5: symbolic string-based base quality summary

############
#  SET UP  #
############

#__STANDARD PERL PACKAGES
   use strict;
   use warnings;
   use Getopt::Long;
   use FileHandle;

#__BASOVAC PACKAGE
   use Statistics::BasoVac;

###################   this section contains the only user-configurable code
#  CONFIGURATION  #   in the program
###################

#__DEFAULTS
   my $defaults = {

   #__NORMAL CONTAMINATION BY TUMOR
      'normal_contam_by_tumor' => 0,

   #__TUMOR CONTAMINATION BY NORMAL
      'tumor_contam_by_normal' => 0,

   #__NORMAL VARIANT FREQUENCY
      'normal_variant_freq' => 0.5,

   #__TUMOR VARIANT FREEQUENCY
      'tumor_variant_freq' => 0.5,

   #__TUMOR BACKGROUND MUTATION RATE
      'tumor_background_mutation_rate' => 0.000001,
   };

####################
#  PRE-PROCESSING  #
####################

#__NEED 64-BIT MACHINE TO RUN SAMTOOLS (PARSE UNAME -A CUZ -P RETURNS "UNKNOWN")
   my $hardware = `uname -a`;
   unless ($hardware =~ /x86\_64/) {
      warn "must run this code on a 64-bit machine --- for example at the " .
       "genome center type 'open64' to get an interactive session on one " .
       "of the 64-bit blade machines\n";
      exit;
   }

#__NON-CONFIGURABLE DEFAULTS FOR COMMAND-LINE ARGS - SEE HELP OUTPUT BELOW
   my ($help, $dir) = (0, ".");
   my ($in_file, $prob) = ({}, {});
   my ($ofile);

#__PARSE COMMAND LINE
   my $status = &GetOptions (
      "dir=s" => \$dir,
      "help" => \$help,
      "nfile=s" => \$in_file->{'nfile'},
      "ofile=s" => \$ofile,
      "rfile=s" => \$in_file->{'rfile'},
      "tfile=s" => \$in_file->{'tfile'},
      "normal_contam_by_tumor=f" => \$prob->{'normal_contam_by_tumor'},
      "tumor_contam_by_normal=f" => \$prob->{'tumor_contam_by_normal'},
      "normal_variant_freq=f" => \$prob->{'normal_variant_freq'},
      "tumor_variant_freq=f" => \$prob->{'tumor_variant_freq'},
      "tumor_background_mutation_rate=f" => \$prob->{'tumor_background_mutation_rate'},
   );
   die "could not parse command line" unless $status;

#__PRINT HELP AND EXIT IF DIRECTED TO
   if ($help) {
      print "=" x 60, "\n  $0\n", "=" x 60, "\n\n";
      print "Synopsis: Bayesian Somatic Variants in Cancer\n";
      print "\n";
      print "File- and Directory-Related Command Line Arguments:\n";
      print "  --dir[string]      output goes to this directory ('$dir')\n";
      print "  --rfile[string]    input human reference file (FASTA file)\n";
      print "  --nfile[string]    input of normal sample (BAM file)\n";
      print "  --ofile[string]    output of this program\n";
      print "  --tfile[string]    input of tumor sample (BAM file)\n";
      print "\n";
      print "Parameter-Related Command Line Arguments:\n";
      print "  --normal_contam_by_tumor[f]  prob of tumor read in normal\n";
      print "  --tumor_contam_by_normal[f]  prob of normal read in tumor\n";
      print "  --normal_variant_freq[f]     normal sample variant frequency\n";
      print "  --tumor_variant_freq[f]      tumor sample variant frequency\n";
      print "  --tumor_background_mutation_rate[f]  prob of base mutation\n";
      print "\n";
      print "Miscellaneous Command Line Arguments:\n";
      print "  --help             print this message and quit\n";
      print "\n";
      print "Notes:\n";
      print "  There is little checking done on the input files, other than\n";
      print "that they exist. Samtools processes these input files into the\n";
      print "parseable ascii 'pileup' format and this information is then\n";
      print "passed to the BasoVac algorithm.\n";
      exit;
   }

#__COMMAND LINE CHECKING: FILES (SPECIFICATION AND EXISTENCE CHECKS)
   foreach my $file (qw/nfile rfile tfile/) {
      die "need to specify a '$file' file" unless defined $in_file->{$file};
      die "cant find '$file' file" unless -e $in_file->{$file};
   }

#__COMMAND LINE CHECKING: PROBABILITY PARAMETERS
   foreach my $param (qw/normal_contam_by_tumor tumor_contam_by_normal normal_variant_freq tumor_variant_freq tumor_background_mutation_rate/) {
######die "need to specify $param" unless defined $prob->{$param};

   #__SET A DEFAULT IF ARG NOT SPECIFIED
      $prob->{$param} = $defaults->{$param} unless defined $prob->{$param};

   #__SANITY CHECK
      die "$param must be a probability" unless
            $prob->{$param} >= 0 && $prob->{$param} <= 1;
   }

#__QUERY FOR MACHINE NAME
   my $machine;
   if (defined $ENV{'HOSTNAME'}) {
      $machine = $ENV{'HOSTNAME'};
   } else {
      $machine = `uname -n`;
      chomp $machine;
   }

#__USE COMMAND-LINE FILE NAME (WITH PATH) FOR THE OUTPUT
   if (defined $ofile && $ofile) {
      $ofile = join ("/", $dir, $ofile);

#__OR MAKE UP ARBITRARY NAME
   } else {
      $ofile = $dir . "/";
      $ofile .= join ("_", "output", "basovac", $$);
      $ofile .= ".txt";
   }

#__MAKE SURE DIRECTORY EXISTS
   if ($dir ne ".") {
      unless (-e $dir && -d $dir) {
         mkdir $dir || die "could not create output directory '$dir'";
      }
   }

#__OPEN OUTPUT FILE
   unlink $ofile if -e $ofile;
   open (OUT, ">$ofile") || die "could not open output file '$ofile'";

#__OUTPUT FILE HEADER
   print OUT "#  BAYESIAN SOMATIC VARIANTS IN CANCER\n#\n";
   print OUT "#  ", "=" x 60, "\n#\n";
   print OUT "#  BasoVac calculation run through script $0\n#\n";
   print OUT "#  run on machine: $machine\n#\n";
   print OUT "#  files\n#\n";
   print OUT "#     reference sequence file:  $in_file->{'rfile'}\n#\n";
   print OUT "#     normal sample sequence file:  $in_file->{'nfile'}\n#\n";
   print OUT "#     tumor sample sequence file:  $in_file->{'tfile'}\n#\n";
   print OUT "#  global parameters\n#\n";
   print OUT "#     normal_contam_by_tumor:         $prob->{'normal_contam_by_tumor'}\n#\n";
   print OUT "#     tumor_contam_by_normal:         $prob->{'tumor_contam_by_normal'}\n#\n";
   print OUT "#     normal_variant_freq:            $prob->{'normal_variant_freq'}\n#\n";
   print OUT "#     tumor_variant_freq:             $prob->{'tumor_variant_freq'}\n#\n";
   print OUT "#     tumor_background_mutation_rate: $prob->{'tumor_background_mutation_rate'}\n#\n";
#  print OUT "#\n";
#  print OUT "#  column descriptors (useful for parsing and collating)\n#\n";
#  print OUT "#     col 1: chromosome number\n#\n";
#  print OUT "#     col 2: chromosome position\n#\n";
   print OUT "#  ", "=" x 60, "\n#\n";

#####################
#  MAIN PROCESSING  #
#####################

#__RUN SAMTOOLS ON NORMAL SAMPLE TO GET "PILEUP" FILE AND FORWARD TO PIPE
   my $fhn = new FileHandle;
   $fhn->open ("samtools pileup -f $in_file->{'rfile'} $in_file->{'nfile'} |")
      || die "could not run samtools on normal BAM file '$in_file->{'nfile'}'";

#__RUN SAMTOOLS ON TUMOR SAMPLE TO GET "PILEUP" FILE AND FORWARD TO PIPE
   my $fht = new FileHandle;
   $fht->open ("samtools pileup -f $in_file->{'rfile'} $in_file->{'tfile'} |")
      || die "could not run samtools on tumor BAM file '$in_file->{'tfile'}'";

#__GET FIRST LINE OF EACH FILE AND ESTABLISH MARKER STATUS
   my ($marker, $curr_pos, $n_ref, $t_ref) = query_both_next ($fhn, $fht);

#__ANALYZE ALL POSITIONS SHARED IN COMMON
   while (1) {

   #__ANALYZE THIS POSITION IF WE'RE SET
      if ($marker eq "set") {

      #__GET PARAMS FROM THE NORMAL SAMPLE
         my ($n_tot_reads, $n_supp_reads, $n_qvals) = parse_pileup ($n_ref);

      #__GET PARAMS FROM THE TUMOR SAMPLE
         my ($t_tot_reads, $t_supp_reads, $t_qvals) = parse_pileup ($t_ref);

      #__NEW BASOVAC OBJECT FOR THIS POSITION
         my $basovac_obj = Statistics::BasoVac->new ({

         #__READ COUNTS
            normal_total_reads     => $n_tot_reads,
            normal_support_reads   => $n_supp_reads,
            tumor_total_reads      => $t_tot_reads,
            tumor_support_reads    => $t_supp_reads,

         #__MUTUAL CONTAMINATION
            normal_contam_by_tumor => $prob->{'normal_contam_by_tumor'},
            tumor_contam_by_normal => $prob->{'tumor_contam_by_normal'},

         #__VARIANT FRQUENCIES
            normal_variant_freq    => $prob->{'normal_variant_freq'},
            tumor_variant_freq     => $prob->{'tumor_variant_freq'},

         #__READ QUALITY VALUES (LIST REFERENCES)
            normal_read_qualities  => $n_qvals,
            tumor_read_qualities   => $t_qvals,

         #__TUMOR BACKGROUND MUTATION RATE
            tumor_background_mutation_rate => $prob->{'tumor_background_mutation_rate'},
         });

      #__GET POSTERIOR PROBABILITY OF ALL THE SCENARIOS
#        print "PLACEHOLDER: calling basovac for position '$curr_pos'\n";
#        my $pval_som = $basovac_obj->prob_somatic_variant;
         my $pval_het = $basovac_obj->prob_heterozygous_variant;
         my $pval_hom = $basovac_obj->prob_homozygous_variant;
         my $pval_loh = $basovac_obj->prob_loh;
         my $pval_nne = $basovac_obj->prob_nne;

      #__OUTPUT
         print OUT "POSITION $curr_pos\n";
         print OUT "  params\n";
         print OUT "    normal_total_reads:    $n_tot_reads\n";
         print OUT "    normal_support_reads:  $n_supp_reads\n";
         print OUT "    tumor_total_reads:     $t_tot_reads\n";
         print OUT "    tumor_support_reads:   $t_supp_reads\n";
         print OUT "    normal_read_qualities: " . join (' ', @{$n_qvals}) . "\n";
         print OUT "    tumor_read_qualities:  " . join (' ', @{$t_qvals}) . "\n";
         print OUT "  resulting P-values\n";
         print OUT "    prob heterozygous variant: $pval_het\n";
         print OUT "    prob homozygous variant:   $pval_hom\n";
         print OUT "    prob loss-of-hetero:       $pval_loh\n";
         print OUT "    prob non-notable event:    $pval_nne\n";
         print OUT "\n";

      #__GET NEXT LINE OF EACH FILE AND ESTABLISH THE MARKER STATUS
         ($marker, $curr_pos, $n_ref, $t_ref) = query_both_next ($fhn, $fht);

   #__ELSE INCREMENT TUMOR UP TO NORMAL IF LATTER CURRENTLY HOLDS THE MARKER
      } elsif ($marker eq "n") {
         ($marker, $curr_pos, $t_ref) =
               query_nonmarker_until ($fht, $marker, $curr_pos);

   #__ELSE INCREMENT NORMAL UP TO TUMOR IF LATTER CURRENTLY HOLDS THE MARKER
      } elsif ($marker eq "t") {
         ($marker, $curr_pos, $n_ref) =
               query_nonmarker_until ($fhn, $marker, $curr_pos);

   #__ELSE WE'RE DONE
      } elsif ($marker eq "done") {
         last;
      }
   }

##################
#  POST PROCESS  #
##################

#__WRAP-UP OUTPUT AND CLOSE FILE
   close (OUT);
   $fhn->close;
   $fht->close;

################################################################################
#                                                                              #
#                            S U B R O U T I N E S                             #
#                                                                              #
################################################################################

#  ============
#  PARSE PILEUP   parse symbolic read & qual strings in samtools pileup format
#  ============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  There is not very much documentation on this format, only the short
#  description at "http://samtools.sourceforge.net/pileup.shtml". Dave Larson
#  also explained this in detail and furnished some of the regular expressions
#  used here.
#
#  some typical lines in pileup format
#
#  2 25310746 C 28 t..$.,.Tt..T,,Tt..tttt,t,t.T^F, G##@F8#HHF/A;F<@/FBF8F=F#FHB
#  2 25310797 G 14 .,..,,,,,,,,..                  FD##?HH8GEFBFH
#  2 25310803 A 8  ,.,..,..                        C3H>)A:H
#  2 25310778 T 12 ,,$,,..,..,..                   CBCCFCHE@>BH
#  2 25310676 C 2  .^].                            =@
#
#  This routine takes a reference to a list of the already-split fields of a
#  single "pileup" line.
#
#  HELP FROM DAVE LARSON
#
# This should convert the bases into Phred scale probabilities that they are
# incorrect. Forgive my abuse of map.
# 
# my @prob_wrong = map { ord($_)-33 } split("","<<<+;<<<<<<<<<<<=<;<;7<&");
# 
# This should strip out any asterixes (deleted bases), end of read
# characters or beginning of read characters. Print was there to verify it
# matched the number of reads in the samtools pileup.
# 
# my $string = q{,.$.....,,.,.,...,,,.,..^+.};
# $string =~ s/([*\$])|([\^].{1})//g;
# print $string,"\t",length($s),"\n";

sub parse_pileup {
   my ($fields) = @_;

#__TOTAL NUMBER OF COVERING READS ASSIGNED DIRECTLY
   my $total_reads = $fields->[3];

#__MUST PARSE SYMBOLS IN READ STRING TO GET NUMBER READS SUPPORTING REF ALLELE
   my $read_string = $fields->[4];

#__FILTER '*', '$', & "CARAT FOLLOWED BY ANY CHARACTER" SYMBOLS FROM READ STRING
   $read_string =~ s/([*\$])|([\^].{1})//g;

#__FILTER ALL INSERTION AND DELETION SUB-STRINGS
#
#  Here, we have to first determine the exact size of the insertion or deletion
#  by reading the field size and then we remove exactly that number of bases
#  following (and along with) the numeric specifier. This process will leave
#  any variant bases intact in that probably rare case of a variant immediately
#  following an insertion/deletion. We do each instance one-by-one because
#  there's no guarantee that multiple instances would be the same size. The
#  basic method is to find an instance, then create a new regexp representing
#  that *specific* instance, i.e. having an exact quantifier "{number}" read
#  from the field size replacing the one-or-more quanitifier "+", so that the
#  exact number of bases are deleted.

   while (1) {
      if ($read_string =~ /([-+])([0-9]+)[ACGTNacgtn]+/) {
         my ($sign, $num) = ($1, $2);
         my $quantified_regexp = quotemeta ($sign);
         $quantified_regexp .= $num . '[ACGTNacgtn]{' . $num . '}';
         $read_string =~ s/$quantified_regexp//;
      } else {
         last;
      }
   }

#__SANITY CHECK: STRING LENGTH MUST BE EQUAL TO NUMBER OF COVERING BASES
   die "string length of '$read_string' not equal to total reads '$total_reads'"
         unless $total_reads == length ($read_string);

#__COUNT VARIANTS IN STRING (SEE EG http://www.perlmonks.org/?node_id=527973)
   my $var_count = () = $read_string =~ /[ACGTNacgtn]/g;

#__NUMBER OF READS SUPPORTING REFERENCE
   my $ref_reads = $total_reads - $var_count;

#__DECODE THE LOCAL QUALITIES FOR EACH COVERING BASE (BASED ON DAVE LARSON'S
#  IMPLEMENTATION OF THE RULE THAT THE SYMBOL'S ASCII CODE MINUS 33 IS THE
#  QUALITY VALUE --- SEE E.G. http://samtools.sourceforge.net/pileup.shtml,
#  WHICH SEEMS TO BE THE ONLY REAL "PILEUP" DOCUMENTATION)

   my $qvals = [map { ord($_)-33 } split("", $fields->[5])];

#__RETURN RESULTS
   return ($total_reads, $ref_reads, $qvals);
}

#  ======================
#  QUERY NON-MARKER UNTIL   read non-marker file until caught-up with other file
#  ======================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub query_nonmarker_until {
   my ($fh, $marker, $curr_pos) = @_;
   my @line;

#__RETURN APPROPRIATE FLAGS IF FILE IS ALREADY EXHAUSTED
   return ("done", 0, []) if $fh->eof;

#__READ LINES OF THE NON-MARKER FILE
   while (<$fh>) {

   #__PARSE LINE
      chomp;
      @line = split '\t';

   #__IF POSITIONS ARE SYNCHRONIZED THEN MARK AS SUCH AND RETURN
      if ($line[1] == $curr_pos) {
         $marker = "set";
         last;

   #__IF POSITION PASSES CURRENT MARK THEN SWITCH MARKERS AND RETURN
      } elsif ($line[1] > $curr_pos) {
         $curr_pos = $line[1];
         if ($marker eq "n") {
            $marker = "t";
         } else {
            $marker = "n";
         }
         last;
      }
   }
   return ($marker, $curr_pos, \@line);
}

#  ===============
#  QUERY BOTH NEXT   read 1 line of each file and set marker flag appropriately
#  ===============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub query_both_next {
   my ($fhn, $fht) = @_;

#__RETURN APPROPRIATE FLAGS IF EITHER FILE IS ALREADY EXHAUSTED
   return ("done", 0, [], []) if $fhn->eof;
   return ("done", 0, [], []) if $fht->eof;

#__GET A LINE (BASE POSITION) OF EACH FILE
   my $norm_line = <$fhn>;
   my $tum_line = <$fht>;
   chomp ($norm_line, $tum_line);

#__NAIVE-PARSE THE LINE
   my @normal_line = split ('\t', $norm_line);
   my @tumor_line = split ('\t', $tum_line);

#__ESTABLISH CURRENT POSITION AS LOWER NUMBER PLAUSIBLY COMMON TO BOTH SETS
   my ($curr_pos, $marker) = ($normal_line[1], "n");
   if ($tumor_line[1] > $curr_pos) {
      $curr_pos = $tumor_line[1];
      $marker = "t";
   } elsif ($tumor_line[1] == $curr_pos) {
      $marker = "set";
   }

#__RETURN MARKER, PROVISIONAL POSITION, AND POINTERS TO THE PARSED LINES
   return ($marker, $curr_pos, \@normal_line, \@tumor_line);
}

