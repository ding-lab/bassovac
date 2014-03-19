#!/gsc/bin/perl
#
#   THIS SCRIPT IS A WRAPPER AROUND THE Bassovac.pm SOMATIC VARIANT CALLER
#
#   TO DO A FAST TEST --- I.E. NO PILEUP RUN BECAUSE WE HARD-CODE ALL OF THAT
#
#   INOF
#
#
#   ./dnmt3a-fasttest.pl --normal_purity 0.95 --tumor_purity 0.97 --ofile r882_dnmt3a_aml31_np_95_tp_97.dat
#

############
#  SET UP  #
############

#__STANDARD PERL PACKAGES
   use strict;
   use warnings;
   use Getopt::Long;
   use FileHandle;

#__BASOVAC PACKAGE
   use Statistics::Bassovac;

###################   this section contains the only user-configurable code
#  CONFIGURATION  #   in the program
###################

#__DEFAULTS
   my $defaults = {

   #__NORMAL PURITY
      'normal_purity' => 1,

   #__TUMOR PURITY
      'tumor_purity' => 1,

   #__NORMAL VARIANT FREQUENCY
      'normal_variant_freq' => 0.5,

   #__TUMOR VARIANT FREEQUENCY
      'tumor_variant_freq' => 0.5,

   #__TUMOR BACKGROUND MUTATION RATE FOR AML
      'tumor_background_mutation_rate' => 0.0000018,
   };

####################
#  PRE-PROCESSING  #
####################

#__NON-CONFIGURABLE DEFAULTS FOR COMMAND-LINE ARGS - SEE HELP OUTPUT BELOW
   my ($help, $dir) = (0, ".");
   my ($in_file, $prob) = ({}, {});
   my ($ofile);

#__PARSE COMMAND LINE
   my $status = &GetOptions (
      "dir=s" => \$dir,
      "help" => \$help,
      "ofile=s" => \$ofile,
      "normal_purity=f" => \$prob->{'normal_purity'},
      "tumor_purity=f" => \$prob->{'tumor_purity'},
      "normal_variant_freq=f" => \$prob->{'normal_variant_freq'},
      "tumor_variant_freq=f" => \$prob->{'tumor_variant_freq'},
      "tumor_background_mutation_rate=f" => \$prob->{'tumor_background_mutation_rate'},
   );
   die "could not parse command line" unless $status;

#__PRINT HELP AND EXIT IF DIRECTED TO
   if ($help) {
      print "=" x 60, "\n  $0\n", "=" x 60, "\n\n";
      print "Synopsis: Bayesian Somatic Variants in Cancer\n";
      print "SPECIAL VERSION TO EXAMINE R882 MUTATION IN DNMT3A\n";
      print "\n";
      print "THIS IS THE *FAST* VERSION THAT IS HARD-CODED --- NO PILE-UP\n";
      print "\n";
      print "File- and Directory-Related Command Line Arguments:\n";
      print "  --dir[string]      output goes to this directory ('$dir')\n";
      print "  --ofile[string]    output of this program\n";
      print "\n";
      print "Parameter-Related Command Line Arguments:\n";
      print "  --normal_purity[f]  prob of 'normal' read bing from normal\n";
      print "  --tumor_purity[f]  prob of 'tumor' read being from tumor\n";
      print "  --normal_variant_freq[f]     normal sample variant frequency\n";
      print "  --tumor_variant_freq[f]      tumor sample variant frequency\n";
      print "  --tumor_background_mutation_rate[f]  prob of base mutation\n";
      print "\n";
      print "Miscellaneous Command Line Arguments:\n";
      print "  --help             print this message and quit\n";
      print "\n";
      exit;
   }

#__COMMAND LINE CHECKING: PROBABILITY PARAMETERS
   foreach my $param (qw/normal_purity tumor_purity normal_variant_freq tumor_variant_freq tumor_background_mutation_rate/) {

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
   print OUT "#  Bassovac calculation run through script $0\n#\n";
   print OUT "#  run on machine: $machine\n#\n";
   print OUT "#  files\n#\n";
   print OUT "#  global parameters\n#\n";
   print OUT "#     normal_purity:         $prob->{'normal_purity'}\n#\n";
   print OUT "#     tumor_purity:         $prob->{'tumor_purity'}\n#\n";
   print OUT "#     normal_variant_freq:            $prob->{'normal_variant_freq'}\n#\n";
   print OUT "#     tumor_variant_freq:             $prob->{'tumor_variant_freq'}\n#\n";
   print OUT "#     tumor_background_mutation_rate: $prob->{'tumor_background_mutation_rate'}\n#\n";
   print OUT "#  ", "=" x 60, "\n#\n";

#####################
#  MAIN PROCESSING  #
#####################

#__ANALYZE ALL POSITIONS SHARED IN COMMON
   my $norm_var_freq = 0.5;

#__ACTUAL DNMT3A IN AML31
#  my $n_qvals = [qw/33 26 18 34 34 32 34 38 34 27 39 31 30 34 39/];
#  my $t_qvals = [qw/38 2 2 31 37 23 2 39 39 37 14 32 26 37 27 31 14 37 33 37 23 37 28 37 2 37 39 33/];
#  my $norm_tot_reads = 15;
#  my $norm_sup_reads = 14;
#  my $tum_tot_reads  = 28;
#  my $tum_sup_reads  = 15;
#  my $bins_normal = 3;
#  my $bins_tumor = 3;

#__SIMULATED CASE -- RIGHT NOW TOO MANY READS (I THINK) CUZ IT DOESN'T RETURN
#  my $n_qvals = [qw/33 26 18 34 34 32 34 38 34 27 39 31 30 34 39 23 34 2 2 2 2 2 2 2/];
#  my $t_qvals = [qw/38 2 2 31 37 23 2 39 39 37 14 32 26 37 27 31 14 37 33 37 23 37 28 37 2 37 39 33 35 35 34 23 2 2 2 2 2 2 15 15 13 17 28 2 2 2/];
#  my $norm_tot_reads = 15;
#  my $norm_sup_reads = 14;
#  my $tum_tot_reads  = 28;
#  my $tum_sup_reads  = 24;
#  my $bins_normal = 2;
#  my $bins_tumor = 2;

#__TEST: WE THOUGHT 2 VS 3 BIN SHOULD BE VERY DIFFERENT BUT ACTUALLY NOT SO MUCH
#  my $n_qvals = [qw/2 2 3 2 2 12 10 10 11 13 13 10 38 38 40 36 32 35 38/];
#  my $t_qvals = [qw/2 2 2 4 2 2 17 16 19 20 17 18 39 39 37 37 38 39 36 34 32/];
   my $n_qvals = [qw/2 2 3 2 2 12 10 10 11 13 13 10 38 38 40 36 32 35 38/];
   my $t_qvals = [];
   for (my $i = 0; $i <= 30; $i++) {
      push @{$t_qvals}, 38;
   }
   for (my $i = 0; $i <= 100; $i++) {
      push @{$t_qvals}, 4;
   }
   my $norm_tot_reads = scalar @{$n_qvals};
   my $norm_sup_reads = scalar @{$n_qvals} - 1;
   my $tum_tot_reads  = scalar @{$t_qvals};
   my $tum_sup_reads  = int (scalar @{$t_qvals} / 2) + 1;
   my $bins_normal = 2;
   my $bins_tumor = 2;

#__RUN BASSOVAC
   foreach my $count (0..100) {
      my $tum_var_freq = $count / 100;

   #__NEW BASOVAC OBJECT FOR THIS POSITION
      my $basovac_obj = Statistics::Bassovac->new ({

      #__READ COUNTS
         normal_total_reads     => $norm_tot_reads,
         normal_support_reads   => $norm_sup_reads,
         tumor_total_reads      => $tum_tot_reads,
         tumor_support_reads    => $tum_sup_reads,

      #__PURITIES
         normal_purity => $prob->{'normal_purity'},
         tumor_purity => $prob->{'tumor_purity'},

      #__VARIANT FRQUENCIES
         normal_variant_freq    => $norm_var_freq,
         tumor_variant_freq     => $tum_var_freq,

      #__NUMBER OF QUALITY BINS
         normal_num_q_bins => $bins_normal,
         tumor_num_q_bins => $bins_tumor,

      #__READ QUALITY VALUES (LIST REFERENCES)
         normal_read_qualities  => $n_qvals,
         tumor_read_qualities   => $t_qvals,

      #__TUMOR BACKGROUND MUTATION RATE
         tumor_background_mutation_rate => $prob->{'tumor_background_mutation_rate'},
      });
# die "SCRIPTL: END OF TEST";

   #__GET POSTERIOR PROBABILITY OF ALL THE SCENARIOS
      my $pval_het = $basovac_obj->prob_heterozygous_variant;
      my $pval_hom = $basovac_obj->prob_homozygous_variant;
      my $pval_loh = $basovac_obj->prob_loh;
      my $pval_nne = $basovac_obj->prob_nne;

   #__OUTPUT INFO HEADER ON FIRST PASS
      if ($count == 0) {
         print OUT "# DNMT3A MUTATION IN FAST SCRIPT\n";
         print OUT "#   params\n";
         print OUT "#     normal_total_reads:    $norm_tot_reads\n";
         print OUT "#     normal_support_reads:  $norm_sup_reads\n";
         print OUT "#     tumor_total_reads:     $tum_tot_reads\n";
         print OUT "#     tumor_support_reads:   $tum_sup_reads\n";
         print OUT "#     normal_read_qualities: " . join (' ', @{$n_qvals}) . "\n";
         print OUT "#     tumor_read_qualities:  " . join (' ', @{$t_qvals}) . "\n";
         print OUT "#     NORMAL VARIANT FREQ:  $norm_var_freq\n";
         print OUT "#     column 1: tumor variant frequency\n";
         print OUT "#     column 2: probability of heterozygous variant\n";
      }

   #__OUTPUT HET VARIANT PROBABILITY FOR THIS TUMOR VARIANT FREQ
      print OUT "$tum_var_freq\t$pval_het\n";
   }

##################
#  POST PROCESS  #
##################

#__WRAP-UP OUTPUT AND CLOSE FILE
   close (OUT);

################################################################################
#                                                                              #
#                            S U B R O U T I N E S                             #
#                                                                              #
################################################################################

