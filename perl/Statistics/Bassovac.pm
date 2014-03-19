package Statistics::Bassovac;

#  SUPERCEDES BasoVac.pm
#
#  1. This package is framed in terms of "purity" rather than "contamination"
#
#  2. Has convolution approximation (2 & 3 bin special cases) for variable read
#     qualities
#
#  3. Has Ramanujan approximation of exact binomial (see below)
#
#  NUMERICAL ARTEFACTS:
#
#  The Poisson p_mass_vector did not work correctly, as implemented (found on
#  Malachi Griffith's PNC4 KRAS mutation) because of a bona-fide round-off error
#  that significantly changes the final answer, so we simply now comment that
#  out. It has been replaced with a Ramanujan approximation of the exact
#  binomial result, which seems to work for the PNC4 KRAS mutation parameters.
#
#  REWORK / CORRECTION / TUMOR SUBCLONE MASS FRACTION
#
#  In November 2012, we overhauled the theory to make proper distinction
#  between purity and tumor subclone mass fraction (see "variant validation"
#  notes, pp 71-94). With respect to programming implementation, these
#  modifications are described in "travis_modification.pdf" (because Travis
#  Abbott will modifiy the production C++ of this code). We made appropriate
#  modifications here for the Perl package.

################################################################################
##                                                                            ##
##                         C O N F I G U R A T I O N                          ##
##                                                                            ##
##    this section contains the only user-configurable code in the package    ##
##                                                                            ##
################################################################################

#__SPECIFY DEFAULTS FOR ALL THOSE VARIABLES THAT ARE NOT DIRECTLY OBSERVED
   use constant DEFAULTS => {

   #__THESE APPLY TO BOTH SAMPLE TYPES (TUMOR AND NORMAL)
      purity => 1,
      num_q_bins => 1,

   #__THESE APPLY ONLY TO THE NORMAL SAMPLE TYPE
      homozygous_variant_rate => 0.0005, # from SoapSNP (li:2009:b:seq:snps)
      heterozygous_variant_rate => 0.001, # from SoapSNP (li:2009:b:seq:snps)

   #__THESE APPLY ONLY TO THE TUMOR SAMPLE TYPE
      background_mutation_rate => 0.000002, # common value
      subclone_mass_fraction => 1, # default is a single dominant clone
   };

#__VERSION
   our $VERSION = 1.2;

################################################################################
##                                                                            ##
##                                S E T - U P                                 ##
##                                                                            ##
################################################################################

#__STANDARD PERL PACKAGES
   use strict;
   use Carp;

#__CONSTANT OF PI
   use constant PI => 4*atan2 1, 1;

################################################################################
##                                                                            ##
##         I N T R O D U C T O R Y   P O D   D O C U M E N T A T I O N        ##
##                                                                            ##
################################################################################

=head1 NAME

Statistics::Bassovac - BAyesian Scheme for SOmatic VAriants in Cancer

=head1 SYNOPSIS

	use Statistics::Bassovac;

	$basovac_obj = Statistics::Bassovac->new ($ref_to_hash_of_arguments);
	$pval = $basovac_obj->prob_heterozygous_variant;
	$pval = $basovac_obj->prob_homozygous_variant;
	$pval = $basovac_obj->prob_somatic_variant;
	$pval = $basovac_obj->prob_loh;
	$pval = $basovac_obj->prob_nne;

=head1 DESCRIPTION

Cancer genomics now routinely involves sequencing tumor
and normal sample pairs, from which the immediate desired
results are somatic events like mutation and loss of
heterozygosity.
Quite a few algorithms have been proposed for normal-only (germline)
mutation detection, but the tumor-normal pair presents an additional problem
of mutual suffusion of the samples, i.e. impurities within the normal sample
from tumor material and I<vice versa>, and tumor heterogeneity, i.e.
"clonality".
This algorithm is designed to account for these computational confounders in
a systematic way, weaving them formally into the well-known method of Bayesian
inversion.

=head2 Nature of the Algorithm

The algorithm is based on the basic idea of dependent analysis, i.e. I<Bayesian
inference>, in which the probability of a specific genotype, given the observed
data, is computed by inverting the probablity of the observed data, given the
genotype.
For example, the scenario for a homozygous variant is:

	                        P(N_11 & T00) * P(DATA | N_11 & T00)
	P (N_11 & T00 | DATA) = ------------------------------------
	                                      P(DATA)

where C<N_11> is homozygous reference in the
normal and C<T00> is homozygous variant in the
tumor.
The primary novelty of the approach is in the proper
determination of the Bernoulli probabilities of each
individual observed read in the tumor and normal read
sets.
Because this is a somatic mutation caller and not a
genotyper, other nucleotide-specific factors, such as
transition-transversion bias, do not have to be
considered.

=head2 Assumptions in the Test

The current implementation makes several assumptions:

=over 4

=item * the tumor and normal data are independent of one another, although
        each depend upon both tumor and normal genotypes via the artefact of
        mutual impurities

=item * suitable estimates of sample purities are available from some
        independent source like pathology reports

=item * mutations in the two alleles of a tumor are independent of one
        another and a function of an estimated background rate - there may
        be biological mechanisms that render invalid this assumption of
        independence

=item * candidate tumor mutation site is diploid, i.e. copy number higher than
	2 is not considered. (Actually, the Bassovac theory already considers
	this factor in a general way. It simply has not been implemented in this
	version of the code.)

=back

=head1 VERSION

	1.0 original implementation

	1.1 extension to use convolution (summation form) for handling
	    distribution of read qualities

	1.2 extension to use tumor subclone mass fraction in a rigorous manner

=head1 LEGALESE

=head2 Author and Copyright

Michael C. Wendl

S<mwendl@wustl.edu>

Copyright (C) 2010, 2011, 2012 Washington University in Saint Louis

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

=head2 Ported Versions

This is the original "research and development" code implementation
for the underlying Bassovac algorithm, but it can also be
used for actual computations if read tallies are not too
large.
Bassovac has been ported to C<C++>, as well and has been tested
against this code by Travis Abbott and that implementation will
be much more useful for larger calculations, e.g. on a whole
genome.

=head1 METHODS

The available methods are as follows.

=cut

################################################################################
##                                                                            ##
##                      P R O G R A M M E R   N O T E S                       ##
##                                                                            ##
################################################################################
#
#  THE OBJECT'S STRUCTURE LOOKS LIKE THIS
#
#  $obj = {
#
#  #__QUANTITIES RELATED TO THE COMPUTED PROBABILITIES
#     probs => {
#
#     #__ALL DATA IE P(D)
#        entire_data => 0.8282,
#
#     #__GENOTYPES AND DATA IE P(D AND S_{Njk} AND S{Tmn}) NESTED AS J,K,M,N
#        dat_gentypes => [
#           [
#              [ [0.03, 0.004], [0.08, 0.01] ],
#              [ [0.067, 0.0002], [0.11, 0.03] ]
#           ],
#           [
#              [ [0.078, 0.034], [0.102, 0.04] ],
#              [ [0.045, 0.011], [0.098, 0.089] ]
#           ],
#        ],
#     },
#
#  #__QUANTITIES RELATED TO THE NORMAL SAMPLE
#     N => {
#
#     #__ARGUMENTS (HAVE COUNTERPARTS FOR TUMOR)
#        purity => 0.8, # e.g. 80% pure normal and 20% impurity from tumor
#        total_reads => 30,
#        support_reads => 25,
#        read_qualities  => [31, 42, 44, 33, ...],
#        num_q_bins  => 2,
#
#     #__ARGUMENTS: NORMAL-ONLY
#        homozygous_variant_rate => 0.0005,
#        heterozygous_variant_rate => 0.001,
#
#     #__DERIVED QUANTITIES (HAVE COUNTERPARTS FOR TUMOR)
#        bin_read_err_pvals => [0.0008313, 0.0001709],
#        bin_read_counts => [12, 18],
#        bin_bernoulli_pvals => [0.003, 0.0005],
#     },
#
#  #__QUANTITIES RELATED TO THE TUMOR SAMPLE
#     T => {
#
#     #__ARGUMENTS (HAVE COUNTERPARTS FOR NORMAL)
#        purity => 0.85, # e.g. 85% pure tumor and 15% impurity from normal
#        total_reads => 65,
#        support_reads => 12,
#        read_qualities   => [34, 36, 40, 29, ...],
#        num_q_bins  => 3,
#
#     #__ARGUMENTS: TUMOR-ONLY
#        background_mutation_rate => 0.000002,
#        subclone_mass_fraction => 0.1,
#
#     #__DERIVED QUANTITIES (HAVE COUNTERPARTS FOR NORMAL)
#        bin_read_err_pvals => [0.6309573, 0.0398107, 0.0050119],
#        bin_read_counts => [13, 20, 32],
#        bin_bernoulli_pvals => [0.00332, 0.000434, 0.000056],
#     },
#  };
#
################################################################################
##                                                                            ##
##                         P U B L I C   M E T H O D S                        ##
##                                                                            ##
################################################################################

#  ===
#  NEW   create a new object
#  ===   ~~~~~~~~~~~~~~~~~~~

=head2 new

This is the usual object constructor, which takes a hash
of several required arguments and many additional optional
ones.
The latter each default to a nominal value if not
specified.

	my $basovac_obj = Statistics::Bassovac->new ({

	#__TOTAL READ COUNTS (MANDATORY)
	   normal_total_reads     => 32,
	   tumor_total_reads      => 46,

	#__NUMBERS OF READS SUPPORTING *REFERENCE* ALLELE (MANDATORY)
	   normal_support_reads   => 27,
	   tumor_support_reads    => 13,

	#__READ QUALITY VALUES (MANDATORY)
	   normal_read_qualities  => [31, 42, 44, 33, ...],
	   tumor_read_qualities   => [34, 36, 40, 29, ...],

	#__PRITY ESTIMATES (OPTIONAL)
	   normal_purity => 0.90,
	   tumor_purity => 0.55,

	#__NUMBER OF BINS TO DISCRETIZE QUALITY VALUES (OPTIONAL)
	   normal_num_q_bins => 2,
	   tumor_num_q_bins => 3,

	#__GERMLINE VARIANT RATE (OPTIONAL)
	   normal_homozygous_variant_rate => 0.0005,
	   normal_heterozygous_variant_rate => 0.001,

	#__TUMOR BACKGROUND MUTATION RATE (OPTIONAL)
	   tumor_background_mutation_rate => 0.000002,

	#__TUMOR SUBCLONE MASS FRACTION (OPTIONAL)
	   tumor_subclone_mass_fraction => 0.04,
	});

The method will silently ignore any other arguments not listed
above.

=cut

sub new {
   my $class = shift;

#__A HASH OF THE RELEVANT ARGUMENTS
   my ($input) = @_;

#__NAMES OF THE INDIVIDUAL ARGUMENTS WE EXPECT
   my @counts = qw/
      normal_total_reads      tumor_total_reads
      normal_support_reads    tumor_support_reads
   /;
   my @floats = qw/
      normal_purity  tumor_purity
      normal_homozygous_variant_rate
      normal_heterozygous_variant_rate
      tumor_background_mutation_rate
      tumor_subclone_mass_fraction
   /;
   my @integer_lists = qw/
      normal_read_qualities  tumor_read_qualities 
   /;
   my @bin_tags = qw/
      normal_num_q_bins   tumor_num_q_bins
   /;

#__INITIALIZE OBJECT
   my $obj = {};

#__PROCESS BIN SPECS (NATURAL NUMBERS; OPTIONAL)
   foreach my $bin_tag (@bin_tags) {

   #__DISTILL THE CONSTITUENT PARTS OF THE ARG NAME
      my ($type, $key) = &text_to_indeces ($bin_tag);

   #__TRY TO PROCESS ARG IF IT IS SPECIFIED
      if (defined $input->{$bin_tag}) {

      #__ARG MUST BE A NATURAL NUMBER
         if ($input->{$bin_tag} =~ /^\d+$/ && $input->{$bin_tag} >= 1) {
            $obj->{$type}->{$key} = $input->{$bin_tag};

      #__OTHERWISE THERE'S SOME EXTERNAL PROBLEM WITH THE CALLER
         } else {
            croak "arg '$bin_tag' must be a natural number >= 1";
         }

   #__OTHERWISE ASSIGN THE DEFAULT
      } else {
         $obj->{$type}->{$key} = DEFAULTS->{$key};
      }
   }

#__PROCESS QUALITY VALUES (LIST; MANDATORY)
   foreach my $list (@integer_lists) {

   #__DISTILL THE CONSTITUENT PARTS OF THE ARG NAME
      my ($type, $key) = &text_to_indeces ($list);

   #__TRY TO PROCESS ARG IF IT IS SPECIFIED
      if (defined $input->{$list}) {

      #__THIS ARG MUST POINT TO A LIST
         croak "arg '$list' not a list" unless ref($input->{$list}) eq 'ARRAY';

      #__STORE THE LIST OF QUALITIES
         $obj->{$type}->{$key} = $input->{$list};

      #__PROCESS INTO BIN ERROR PROBABILITIES FOR ACTUAL CALCULATION LATER
         ($obj->{$type}->{"bin_read_err_pvals"},
          $obj->{$type}->{"bin_read_counts"}) = &bin_read_counts_and_error_probs (
            $input->{$list}, $obj->{$type}->{"num_q_bins"}
         );

   #__OTHERWISE THIS IS A FATAL ERROR
      } else {
         croak "missing arg '$list'";
      }
   }

#__PROCESS LATENT AND ESTIMATED VARIABLES (FLOATS; OPTIONAL)
   foreach my $float (@floats) {

   #__DISTILL THE CONSTITUENT PARTS OF THE ARG NAME
      my ($type, $key) = &text_to_indeces ($float);

   #__TRY TO PROCESS ARG IF IT IS SPECIFIED
      if (defined $input->{$float}) {

      #__EACH FLOAT IS FUNDAMENTALLY A PROBABLITY
         if ($input->{$float} =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ &&
                   $input->{$float} >= 0 && $input->{$float} <= 1) {
            $obj->{$type}->{$key} = $input->{$float};

      #__OTHERWISE THERE'S SOME EXTERNAL PROBLEM WITH THE CALLER
         } else {
            croak "arg '$float' must be a probability and '$input->{$float}' " .
                  "does not appear to conform --- please contact the authors " .
                  "if you feel this message is in error";
         }

   #__OTHERWISE ASSIGN THE DEFAULT
      } else {
         $obj->{$type}->{$key} = DEFAULTS->{$key};
      }
   }

#__PROCESS OBSERVED READ COUNTS (NATURAL NUMBERS; REQUIRED)
   foreach my $count (@counts) {

   #__DISTILL THE CONSTITUENT PARTS OF THE ARG NAME
      my ($type, $key) = &text_to_indeces ($count);

   #__TRY TO PROCESS ARG IF IT IS SPECIFIED
      if (defined $input->{$count}) {

      #__ARG MUST BE A NATURAL NUMBER
         if ($input->{$count} =~ /^\d+$/ && $input->{$count} >= 0) {
            $obj->{$type}->{$key} = $input->{$count};

      #__OTHERWISE THERE'S SOME EXTERNAL PROBLEM WITH THE CALLER
         } else {
            croak "arg '$count' must be a natural number >= 0";
         }

   #__OTHERWISE THIS IS A FATAL ERROR
      } else {
         croak "missing arg '$count'";
      }
   }

#__SANITY CHECK: NUMBER OF QUALITY VALUES EQUALS NUMBER OF TOTAL READS

#  TBD

#__SANITY CHECK: NUMBER OF SUPPORTING READS CANNOT EXCEED NUMBER OF TOTAL READS
   foreach my $type (qw/T N/) {
      croak "sample type $type: supporting reads cannot exceed total reads"
         if $obj->{$type}->{'support_reads'} > $obj->{$type}->{'total_reads'};
   }

#__BLESS INTO CLASS AND RETURN OBJECT
   bless $obj, $class;
   return $obj;
}

#  =======================
#  PROB HOMOZYGOUS VARIANT   Bayes' probability of a homozygous variant
#  =======================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 prob_homozygous_variant

Gives the probability of a homozygous variant, i.e. the configuration
where the normal is homozygous reference and the tumor is homozygous
non-reference.

=cut

sub prob_homozygous_variant {
   my $obj = shift;
   my $prob_data = $obj->probability_of_data;
   return $obj->{'probs'}->{'dat_gentypes'}->[1]->[1]->[0]->[0] / $prob_data;
}

#  =========================
#  PROB HETEROZYGOUS VARIANT   Bayes' probability of a heterozygous variant
#  =========================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 prob_heterozygous_variant

Gives the probability of a heterozygous variant, i.e. the configuration where
the normal is homozygous reference and one or the other allele of the tumor is
non-reference.

=cut

sub prob_heterozygous_variant {
   my $obj = shift;
   my $prob_data = $obj->probability_of_data;
   my $phv = $obj->{'probs'}->{'dat_gentypes'}->[1]->[1]->[0]->[1] / $prob_data;
   $phv  +=  $obj->{'probs'}->{'dat_gentypes'}->[1]->[1]->[1]->[0] / $prob_data;
# DEBUG
# print "   P(het) = $phv\n";
# DEBUG
   return $phv;
}

#  ====================
#  PROB SOMATIC VARIANT   Bayes' probability of a het or homoezygous variant
#  ====================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 prob_somatic_variant

This gives the probablity of either a homozygous or heterozygous
variant, i.e. it is the sum of the probabilities of the two
scenarios.

=cut

sub prob_somatic_variant {
   my $obj = shift;
   my $somatic_prob = $obj->prob_homozygous_variant;
   $somatic_prob += $obj->prob_heterozygous_variant;
   return $somatic_prob;
}

#  ========
#  PROB LOH   Bayes' probability of loss-of-heterozygosity variant
#  ========   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 prob_loh

Gives the probability of a loss-of-heterozygosity event, i.e. where the
normal sample is heterozygous, but the tumor is homozygous for the
reference.

=cut

sub prob_loh {
   my $obj = shift;
   my $prob_data = $obj->probability_of_data;
   my $loh = $obj->{'probs'}->{'dat_gentypes'}->[1]->[0]->[1]->[1] / $prob_data;
   $loh  +=  $obj->{'probs'}->{'dat_gentypes'}->[0]->[1]->[1]->[1] / $prob_data;
   return $loh;
}

#  ========
#  PROB NNE   Bayes' probability of a Non-Notable Event
#  ========   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Algorithm note: we do not simply just return 1 - (sum of notable events).
#  Explicit calculation of all configurations allows the user to sum the
#  probs from all cases to check if they sum to 1

=head2 prob_nne

Gives the probability of a non-notable event, which is any scenario
that is I<not> a type of somatic mutation or loss-of-heterozygosity
event.

=cut

sub prob_nne {
   my $obj = shift;
   my $prob_data = $obj->probability_of_data;
   my $nne = $obj->{'probs'}->{'dat_gentypes'}->[1]->[1]->[1]->[1] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[0]->[1]->[0]->[1] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[0]->[1]->[1]->[0] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[0]->[1]->[0]->[0] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[1]->[0]->[0]->[1] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[1]->[0]->[1]->[0] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[1]->[0]->[0]->[0] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[0]->[0]->[1]->[1] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[0]->[0]->[0]->[1] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[0]->[0]->[1]->[0] / $prob_data;
   $nne  +=  $obj->{'probs'}->{'dat_gentypes'}->[0]->[0]->[0]->[0] / $prob_data;
   return $nne;
}

################################################################################
##                                                                            ##
##                       P R I V A T E   M E T H O D S                        ##
##                                                                            ##
################################################################################

###########################################
#  MODEL-SPECIFIC MATHEMATICAL FUNCTIONS  #
###########################################

#  ===================
#  PROBABILITY OF DATA   calculates prob of each case and overall prob of data
#  ===================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  calculates the entire probability of observing the data from considering
#  all possible 16 composite binary genotypes
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way

sub probability_of_data {
   my $obj = shift;

#__JUST RETURN THE VALUE IF WE'VE ALREADY COMPUTED IT
   return $obj->{'probs'}->{'entire_data'}
      if defined $obj->{'probs'}->{'entire_data'};

#__ALL POSSIBLE COMPOSITE BINARY GENOTYPES: NORMAL SAMPLE INDEX 1
   my $prob_of_data = 0;
   foreach my $j (0, 1) {

   #__NORMAL SAMPLE INDEX 2
      foreach my $k (0, 1) {

      #__TUMOR SAMPLE INDEX 1
         foreach my $m (0, 1) {

         #__TUMOR SAMPLE INDEX 2
            foreach my $n (0, 1) {
# DEBUG_06DEC
# print "  (J,K,M,N) = ($j, $k, $m, $n)\n";
# DEBUG_06DEC

            #__COMPOSITE PRIOR OF NORMAL AND TUMOR GENOTYPES OF THIS (J,K,M,N)
               my $prior_probs = $obj->prior_prob_genotypes ($j, $k, $m, $n);
# DEBUG: OK TO HERE FOR MALACHI'S PNC4 KRAS MUTATION

            #__PROBABILITY OF DATA ON NORMAL GIVEN THE COMPOSITE GENOTYPES
               my $prob_data_norm =
                  $obj->prob_observed_data_given_genotypes('N', $j, $k, $m, $n);
# DEBUG_06DEC
# print "      prob normal data given genotype: $prob_data_norm\n";
# DEBUG_06DEC

            #__PROBABILITY OF DATA ON TUMOR GIVEN THE COMPOSITE GENOTYPES
            #  (INDEX SWITCH: SEE "VARIANT VALIDATION" NOTES PP 23)
               my $prob_data_tum =
                  $obj->prob_observed_data_given_genotypes('T', $m, $n, $j, $k);
# DEBUG_06DEC
# print "      prob tumor  data given genotype: $prob_data_tum\n";
# DEBUG_06DEC

            #__JOINT PROB OF DATA & GENOTYPES: P(D & S_{Njk} & S_{Tmn})
               my $prob_joint = $prior_probs * $prob_data_norm * $prob_data_tum;

# DEBUG_06DEC
# print "   (J,K,M,N) = ($j, $k, $m, $n): joint P = $prob_joint\n";
# DEBUG_06DEC
            #__SAVE THIS TO DATA STRUCTURE
               $obj->{'probs'}->{'dat_gentypes'}->[$j]->[$k]->[$m]->[$n]
                  = $prob_joint;

            #__INCREMENT P(D)
               $prob_of_data += $prob_joint;
            }
         }
      }
   }
# DEBUG2
#  die "end of checking";
# DEBUG2

#__SANITY CHECKS: MUST BE A P-VAL BUT ALSO CANNOT BE ZERO
   croak "probability of data must be > 0" unless $prob_of_data > 0;
   croak "probability of data cannot exceed 1" unless $prob_of_data <= 1;
# DEBUG_06DEC
# print "   total P(D) = $prob_of_data\n";
# DEBUG_06DEC

#__ALL OK: SAVE THE VALUE AND RETURN IT
   $obj->{'probs'}->{'entire_data'} = $prob_of_data;
   return $prob_of_data;
}

#  ====================
#  PRIOR_PROB_GENOTYPES   prior probability of genotypes of the paired samples
#  ====================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  calculates the prior probability of the intersection of normal and tumor, ie
#              _
#  P (S_{Njk} | | S_{Tmn}) = P (S_{Njk}) * P(S_{Tmn} | S_{Njk}),
#
#  given specific values of (j, k, m, n)
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way

sub prior_prob_genotypes {
   my $obj = shift;

#__ARGS: (J,K) ARE INDECES OF NORMAL GENOTYPE AND (M,N) ARE INDECES OF THE TUMOR
   my ($j, $k, $m, $n) = @_;

#__NORMAL PRIOR IS EITHER HOMOZYGOUS
   my $composite_prior;
   if ($j == $k) {

   #__REFERENCE: HETERO RATE ALREADY INCLUDES BOTH 0,1 & 1,0 SO NO FACTOR OF 2
      if ($j) {
         $composite_prior = 1 - $obj->{'N'}->{'homozygous_variant_rate'} -
                                  $obj->{'N'}->{'heterozygous_variant_rate'};

   #__VARIANT: SIMPLY THE HOMOZYGOUS VARIANT RATE
      } else {
         $composite_prior = $obj->{'N'}->{'homozygous_variant_rate'};
      }

#__OR HETEROZYGOUS: THIS IS A SPECIFIC CASE (0,1 OR 1,0) SO MUST DIVIDE BY 2
   } else {
      $composite_prior = $obj->{'N'}->{'heterozygous_variant_rate'} / 2;
   }

#__TUMOR PRIOR IS CONDITIONED UPON NORMAL CONFIGURATION
   $composite_prior *= $obj->piecewise_psi ($j, $m);
   $composite_prior *= $obj->piecewise_psi ($k, $n);
# DEBUG: OK TO HERE FOR MALACHI'S PNC4 KRAS MUTATION
# print "   (J,K,M,N) = ($j, $k, $m, $n): composite prior prob = $composite_prior\n";
# DEBUG

#__RETURN COMPOSITE PRIOR OF THE GIVEN NORMAL AND GIVEN TUMOR BINARY GENOTYPES
   return $composite_prior;
}

#  ==================================
#  PROB_OBSERVED_DATA_GIVEN_GENOTYPES   probability of data for sample type i
#  ==================================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  This is the probability of the observed data for sample type i, given
#  the specific genotypes for both samples as specified in the argument list.
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way

sub prob_observed_data_given_genotypes {
   my $obj = shift;
   my ($i, $j, $k, $m, $n) = @_;

#__SORT OUT THE TYPES
   my $l = "N";
   $l = "T" if $i eq "N";

#__BERNOULLI PROBABILITY OF A READ REPORTING THE REFERENCE ALLELE
# DEBUG2
# print "BERNOULLI NUM $i ";
# DEBUG2
   $obj->{$i}->{"bin_bernoulli_pvals"} =
         $obj->bernoulli_read_prob_ref_allele ($i, $j, $k, $l, $m, $n);

#__PREP
   my $pval;
   my $num_bins_i = scalar @{$obj->{$i}->{"bin_bernoulli_pvals"}};

#__CALCULATE P-VALUE WITHOUT CONVOLUTION IF ONLY 1-BIN
   if ($num_bins_i == 1) {

   #__CALCULATE PROB MASS VECTOR FOR JUST THIS BIN: NOT MOST EFFICIENT APPROACH
      my $p_bernoul = $obj->{$i}->{"bin_bernoulli_pvals"}->[0];
      my $n_reads = $obj->{$i}->{"total_reads"};
      my $k_reads = $obj->{$i}->{"support_reads"};
# DEBUG_06DEC
# print "        bernoulli prob of ref allele: $p_bernoul\n";
# print "        total number of reads for sample: $n_reads\n";
# print "        number of reads supporting ref: $k_reads\n";
# DEBUG_06DEC
      my $p_masses = &p_mass_vector ($p_bernoul, $n_reads, $k_reads);

   #__SINCE N >= K FOR THIS CASE THE P-VAL WE WANT IS SIMPLY THE K-TH INDEX
      $pval = $p_masses->[$k_reads];
# DEBUG_06DEC
# print "        P(data|genotype) = $pval\n";
# DEBUG_06DEC

#__ELSE DO 2-BIN CONVOLUTION
   } elsif ($num_bins_i == 2) {
      $pval = $obj->probability_2_bin_summation ($i);

#__ELSE DO 3-BIN CONVOLUTION
   } elsif ($num_bins_i == 3) {
      $pval = $obj->probability_3_bin_summation ($i);

#__ELSE CROAK SINCE WE'RE NOT SET-UP FOR ANYTHING ELSE
   } else {
      croak "we can only handle 1, 2, or 3 bins at the moment --- please " .
            "change argument(s)";
   }

#__FINAL DISPOSITION
   if ($pval >= 0 && $pval <= 1) {
      return $pval;
   } else {
      croak "cannot compute probability of observed data, given genotypes " .
            "for configuration (i, j, k, l, m, n) = ($i, $j, $k, $l, $m, $n)";
   }
}

#  ==============================
#  BERNOULLI_READ_PROB_REF_ALLELE   Bernoulli prob that read reports ref allele
#  ==============================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way
#
#  Fundamentally, the read error probability is a property inherent to the
#  read itself, irrespective of what sample it actually comes from. Were we
#  to be doing the exact, read-specific (brute-force) expansion, there would
#  be no ambiguity whatsoever. However, since we are doing the convolution
#  approximation, a particular read's error value is really a composite of
#  those in the same bin. Here, we use the error p-value of the bin associated
#  with sample i for both sample i itself (purity) *and* sample l (1 - purity),
#  according to the fact that shrinking bins to single reads would converge
#  to the correct result for the single read itself *and* because there's no
#  guarantee that there would be a suitable bin value for sample l. This is
#  opposite of what we did in the earliest version. Dave Larson and I both
#  felt this was the right choice.

sub bernoulli_read_prob_ref_allele {
   my $obj = shift;
   my ($i, $j, $k, $l, $m, $n) = @_;

##_GET VALUES OF THE SWITCHING FUNCTION
## my $phi_ijk = $obj->piecewise_phi ($i, $j, $k);
## my $phi_lmn = $obj->piecewise_phi ($l, $m, $n);

#__GET VALUES OF SWITCHING FUNCTION (NEW CORRECTED VERSION WITH TUMOR SUBCLONE
#  MASS FRACTION) --- IMPLEMENTATION FOLLOWS (COMPLETELY GENERAL) THEOREM 2 IN
#  SUPPLEMENTAL INFORMATION OF THE PAPER RATHER THAN THE SLIGHTLY MORE CONCRETE
#  RESULT IN "VARIANT VALIDATION" NOTES ON PP 94
   my $phi_1_ijk = $obj->piecewise_phi_1 ($i, $j, $k);
   my $phi_2_ijkl = $obj->piecewise_phi_2 ($i, $l, $m, $n);
# DEBUG2
#  print "  (i,j,k,l,m,n) = ($i, $j, $k, $l, $m, $n): numerator = $phi_1_ijk + $phi_2_ijkl\n";
# DEBUG2

#__CALCULATE THE BERNOULLI PROBABILITY FOR EACH CONVOLUTION BIN
   my $bernoulli_pvals = [];
   foreach my $read_err_pval (@{$obj->{$i}->{"bin_read_err_pvals"}}) {

   #__CALCULATE BERNOULLI P-VAL OF A READ OF THIS TYPE REPORTING THE REF ALLELE
#     my $pval = $obj->{$i}->{'purity'} *
#        ((1 - $phi_ijk)*(1 - $read_err_pval) + $read_err_pval * $phi_ijk / 3) +
#        (1 - $obj->{$i}->{'purity'}) *
#        ((1 - $phi_lmn)*(1 - $read_err_pval) + $read_err_pval * $phi_lmn / 3);
####    THIS IS THE BERNOULLI PVAL OF THE *VARIANT* ALLELE
####  my $pval = (1 - 4*$read_err_pval/3) * ($phi_1_ijk + $phi_2_ijkl) / 2
####             + $read_err_pval;
####    THIS IS THE BERNOULLI PVAL OF THE *VARIANT* ALLELE

      my $pval = 1 - (1 - 4*$read_err_pval/3) * ($phi_1_ijk + $phi_2_ijkl) / 2
                 - $read_err_pval;

   #__SANITY CHECK
      croak "bernoulli value '$pval' not a p-val" unless
            $pval >= 0 && $pval <= 1;
#########
print "S($i, $j, $k) /\\ S($l, $m, $n): err = $read_err_pval,  P = $pval\n";
# <STDIN>;
#########

   #__SAVE TO LIST
      push @{$bernoulli_pvals}, $pval;
   }

#__RETURN VALUE
   return $bernoulli_pvals;
}

#################################################
#  MATHEMATICAL METHODS FOR PROBABILITY MASSES  #
#################################################

#  ===========================
#  PROBABILITY 2 BIN SUMMATION   summation implementation for 2 convolution bins
#  ===========================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  SEE BASCON NOTES PP 6
#
#  This routine does the naive summation approach for 2-bin convolution. This
#  should be changed to use the Fourier approach (DFT) when you have time.

sub probability_2_bin_summation {
   my $obj = shift;
   my ($type) = @_;

#__BIN BERNOULLI SUCCESS PROBABILITIES
   my ($p1, $p2) = @{$obj->{$type}->{"bin_bernoulli_pvals"}};

#__BIN READ COUNTS
   my ($n1, $n2) = @{$obj->{$type}->{"bin_read_counts"}};

#__NUMBER OF READS SUPPORTING REFERENCE
   my $k = $obj->{$type}->{'support_reads'};

#__GET BOTH BIN PROBABILITY MASS VECTORS
   my $pvals_bin_1 = &p_mass_vector ($p1, $n1, $k);
   my $pvals_bin_2 = &p_mass_vector ($p2, $n2, $k);

#__EXECUTE THE CONVOLUTION
   my $pval = 0;
   for (my $i = 0; $i <= $k; $i++) {
      $pval += $pvals_bin_1->[$i] * $pvals_bin_2->[$k - $i];
   }
   return $pval;
}

#  ===========================
#  PROBABILITY 3 BIN SUMMATION   summation implementation for 3 convolution bins
#  ===========================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  SEE BASCON NOTES PP 6
#
#  This routine does the naive summation approach for 3-bin convolution. This
#  should be changed to use the Fourier approach (DFT) when you have time.

sub probability_3_bin_summation {
   my $obj = shift;
   my ($type) = @_;

#__BIN BERNOULLI SUCCESS PROBABILITIES
   my ($p1, $p2, $p3) = @{$obj->{$type}->{"bin_bernoulli_pvals"}};

#__BIN READ COUNTS
   my ($n1, $n2, $n3) = @{$obj->{$type}->{"bin_read_counts"}};

#__NUMBER OF READS SUPPORTING REFERENCE
   my $k = $obj->{$type}->{'support_reads'};

#__GET ALL BIN PROBABILITY MASS VECTORS
   my $pvals_bin_1 = &p_mass_vector ($p1, $n1, $k);
   my $pvals_bin_2 = &p_mass_vector ($p2, $n2, $k);
   my $pvals_bin_3 = &p_mass_vector ($p3, $n3, $k);

#__EXECUTE THE CONVOLUTION
   my $pval = 0;
   for (my $i = 0; $i <= $k; $i++) {
      my $sub_sum_pval = 0;
      for (my $j = 0; $j <= $k - $i; $j++) {
         $sub_sum_pval += $pvals_bin_1->[$j] * $pvals_bin_2->[$k - $i - $j];
      }
      $pval += $sub_sum_pval * $pvals_bin_3->[$i];
   }
   return $pval;
}

#  =============
#  P MASS VECTOR   compute probability mass vector for convolution calculation
#  =============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Computes probability masses for a bin from 0 to k using hierarchy of methods,
#  (binomial, Poisson, asymptotic Poisson), where there is appropriate 0-padding
#  if k exceeds the bin size, whereby the contribution of that term within the
#  convolution will correctly be 0. This routine is not completely bullet-proof,
#  for example, we compute binomial coefficients using row-wise left-to-right
#  recursion in Pascal's triangle, so once this overflows there's no going back,
#  even though terms sufficiently further to the right are certain to be
#  computable. However, in many cases, one of the other two terms in the
#  product, "p" or "q", will underflow prior to this point (we have observed
#  this), so the issue may not be too serious.
#
#  Specifically, we calculate binomial coefficients by marching across the
#  required row in Pascal's Triangle, using the following recursion
#
#               n!                            n!
#  C(n,k) = ---------        C(n,k-1) = ---------------
#           (n-k)! k!                   (n-k+1)! (k-1)!
#
#
#    C(n,k)     n! (n-k+1)! (k-1)!    (n-k+1) (n-k)! (k-1)!     n-k+1
#  ---------- = ------------------- = ---------------------- = -------
#   C(n,k-1)       (n-k)! k! n!           (n-k)! k (k-1)!         k
#
#
#  therefore:  C(n,k) = (n-k+1) C(n,k-1) / k    where C(n,0) = 1
#
#  Limited on standard hardware to about C(1028,514) = 7.1560510548779e+307
#  which is within the 1.797693134862316e+308 upper limit i.e. the value of
#  (2-2^{-52})*2^{1023} --- Perl returns 'inf' if this is exceeded
#
#
#  \BEGIN{DEPRECATED}
#
#      THE ROUTINE NO LONGER DOES THIS PART BECAUSE THERE ARE CASES THAT
#     SILENTLY RETURN A COMPLETELY INCORRECT ANSWER BECAUSE OF ROUND-OFF
#     ERROR IN POISSON APPROX. FOR THE MOMENT, WE HAVE REPAIRED THIS PROBLEM
#     BY DOING A RAMANUJAN APPROXIMATION (BELOW) ON THE EXACT BINOMIAL, BUT
#     EVEN THIS HAS NOT BEEN THOROUGHLY TESTED AND MAY STILL FAIL NUMERICALLY
#     FOR CERTAIN COMBINATIONS OF PARAMETERS.
#
#     \BEGIN{NEW_RAMANUJAN}
#
#                          n!       k  /     \ n-k
#         P(n,k,p)  =  --------- * p  | 1 - p |
#                      (n-k)! k!       \     /
#
#                             _                              _
#                            |      n!       k  /     \ n-k   |
#                        ln  |  --------- * p  | 1 - p |      |
#                            |_ (n-k)! k!       \     /      _|
#                   =  e
#
#                                     k            n-k
#                        ln(n!) + ln(p  ) + ln(1-p)   - ln(k!) - ln[(n-k)!]
#                   =  e
#
#                     EVAL THESE USING RAMANUJAN APPROX.
#                        ______________/\_____________
#                       /                             \
#                        ln(n!) - ln(k!) - ln[(n-k)!] + k*ln(p) + (n-k)ln(1-p)
#                   =  e
#
#     \END{NEW_RAMANUJAN}
#
#  In case of overflow, as can happen when the total number of reads > 1028,
#  we can approximate as
#                                  k                                    k
#              n!                 n                           k    (p*n)
#  C(n,k) = --------- \approx ----------   whereby    C(n,k)*p  = ---------
#           k! (n-k)!             k!                                 k!
#
#  which leads directly to the Poisson distribution once asymtotic approximation
#  is also applied to the remaining part of the binomial equation. If "p" is
#  sufficiently small, then p*n may be small enough to raise to the k power
#  without overflowing. WE HAVE NOT MATHEMATICALLY EXAMINED THE VALUES FOR
#  WHICH THIS WILL RESOLVE A BINOMIAL OVERFLOW.
#
#  Failing that, this implementation also fields an exponential-log form of the
#  Poisson eqn, coupled with an asymptotic approximation of the factorial term,
#  which uses Ramanujan's approximation, to compute the probability. See, for
#  example, "plus-minus test" notes pp 29 for an explanation.
#
#  If Ramanujan's approximation proves yet insufficient for certain cases, we
#  may have to implement higher-order asymptotics, e.g. as on pp 111 in
#  Knuth "The Art of Computer Programming, Fundamental Algorithms" Vol 1, or
#  go to extended precision arithmetic, for which there would be much-increased
#  CPU cost.
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way
#
#  \END{DEPRECATED}

sub p_mass_vector {
   my ($p, $n, $k) = @_;

#__INITIALIZE VECTOR
   my $p_masses = [];
   for (my $i = 0; $i <= $k; $i++) {
      push @{$p_masses}, 0;
   }

#__PROBABILITY ELEMENTS RELATED TO BINOMIAL MODEL
   my ($binom_coeff, $p_term, $q_term) = (1, 1, (1 - $p)**$n);
# DEBUG
#   print "        p_term: $p_term\n";
#   print "        q_term: $q_term\n";
# DEBUG

#__PROBABILITY ELEMENTS RELATED TO POISSON MODEL
   my ($mu, $mu_sup_i, $factorial) = ($n*$p, 1, 1);
   my $e_minus_mu = exp (-$mu);

#__PROBABILITY ELEMENTS RELATED TO ASYMPTOTIC POISSON MODEL
   my $log_mu = log ($mu);

#__INITIAL MASS ELEMENT (I=0): BINOMIAL FIRST THEN POISSON ELSE LEAVE AS ZERO
   if ($q_term > 0) {
      $p_masses->[0] = $binom_coeff * $p_term * $q_term;
   } elsif ($e_minus_mu > 0) {
      $p_masses->[0] = $e_minus_mu;
   }
# DEBUG
# print "          i = 0 (first element): $p_masses->[0]\n";
# DEBUG

#__FILL MASS VECTOR WITH REMAINING ELEMENTS
   for (my $i = 1; $i <= $k; $i++) {
      my $pmass;

   #__BINOMIAL COMPONENTS
      $binom_coeff *= ($n - $i + 1) / $i;
      $p_term *= $p;
      $q_term /= 1 - $p;

   #__POISSON COMPONENTS
      $factorial *= $i;
      $mu_sup_i *= $mu;

   #__MASS TERM: TRY EXACT BINOMIAL FIRST
      if ($binom_coeff ne 'inf' && $p_term > 0 && $q_term > 0) {
         $pmass = $binom_coeff * $p_term * $q_term;

# DEBUG
# print "          i = $i: exact binomial\n";
# print "          i = $i: binom_coeff is $binom_coeff   p_term is $p_term   q_term is $q_term\n";
# DEBUG
# old   #__MASS TERM: ELSE TRY POISSON APPROXIMATION
# old      } elsif ($factorial ne 'inf' && $mu_sup_i ne 'inf' &&
# old               $mu_sup_i > 0 && $e_minus_mu > 0) {
# old# DEBUG
# old# print "          i = $i: poisson approx\n";
# old# print "              mu = $mu\n";
# old# print "              e^(-mu) = $e_minus_mu\n";
# old# print "              k = $k\n";
# old# print "              mu^k = $mu_sup_i\n";
# old# print "              k! = $factorial\n";
# old# DEBUG
# old         croak "Poisson approximation cannot be trusted as of PNC4 KRAS case\n no P-value returned\n data are too deep and/or read qualities too high\n";
# old         $pmass = $mu_sup_i * $e_minus_mu / $factorial;
# old
# old   #__MASS TERM: ELSE TRY ASYMPTOTIC POISSON AND LEAVE AS 0 IF STILL UNDERFLOW
# old      } else {
# old         croak "Ramanujan approximation of Poisson approximation cannot be trusted as of PNC4 KRAS case\n no P-value returned\n data are too deep and/or read qualities too high\n";
# old# DEBUG
# old# print "          i = $i: ramanujan approx\n";
# old# DEBUG
# old
# old      #__RAMANUJAN APPROXIMATION FOR LOG(I!)
# old         my $ln_i_factorial = $i * (log($i) - 1)
# old                              + log($i * (1 + 4*$i*(1 + 2*$i))) / 6
# old                              + log (PI) / 2;
# old
# old      #__COMPUTE PVAL FROM EXPONENTIATED FORM
# old         $pmass = exp ($i * $log_mu - $mu - $ln_i_factorial);
# ELSE TRY RAMANUJAN APPROX OF EXACT BINOMIAL
      } else {
# DEBUG
# print "          i = $i: ramanujan approx to binomial\n";
# DEBUG
         my $arg = ramanujan_log_factorial ($n);
         $arg += $k * log ($p);
         $arg += ($n - $k) * log (1 - $p);
         $arg -= ramanujan_log_factorial ($k);
         $arg -= ramanujan_log_factorial ($n - $k);
         $pmass = exp ($arg);
      }

   #__CHECK LEGITMACY OF PROBABILITY MASS AND SAVE IT
      if ($pmass >= 0 && $pmass <= 1) {
         $p_masses->[$i] = $pmass;

   #__OR CROAK
      } else {
         croak "cannot compute prob mass for (p, n, k, i) = ($p, $n, $k, $i)";
      }

   #__LEAVE IF STACK (THIS ROW OF PASCAL'S TRIANGLE) IS COMPLETE
      last if $i == $n;
   }

#__RETURN RESULT
   return $p_masses;
}

sub ramanujan_log_factorial {
   my ($n) = @_;
   return 0 if $n == 0;
   return $n * log($n) - $n + log ($n * (1 + 4 * $n * (1 + 2 * $n)))/6 +
           log (PI) / 2;
}

######################################
#  LOW-LEVEL MATHEMATICAL FUNCTIONS  #
######################################

#  =============
#  PIECEWISE PSI   computes the special switching function called "psi"
#  =============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way

sub piecewise_psi {
   my $obj = shift;
   my ($j, $m) = @_;
   if ($j == $m) {
      if ($j) {
         return 1 - $obj->{'T'}->{'background_mutation_rate'};
      } else {
         return 1 - $obj->{'T'}->{'background_mutation_rate'}/3;
      }
   } else {
      if ($j) {
         return $obj->{'T'}->{'background_mutation_rate'};
      } else {
         return $obj->{'T'}->{'background_mutation_rate'}/3;
      }
   }
}

#  ===============
#  PIECEWISE PHI 2   computes the special switching function called "phi_2"
#  ===============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way
#
#  Implementation follows (completely general) Theorem 2 in supplemental
#  information of the paper rather than the slightly more concrete result in
#  "variant validation" notes on PP 94, but was checked case-by-case against
#  an example of p_T = 0.4, p_N = 0.3, mass fraction = 0.9 against the long-hand
#  results on pp 92-93 of those notes. SEE piecewise_phi_1 FOR EXAMPLE OUTPUT
#  OF A PARTICULAR TEST.

sub piecewise_phi_2 {
   my $obj = shift;
   my ($i, $l, $m, $n) = @_;
   my $phi_4 = 1;
   $phi_4 = $obj->{'T'}->{'subclone_mass_fraction'} if $l eq 'T';
   if ($m == $n) {
      if ($m) {
         return 0;
      } else {
         return 2 * (1 - $obj->{$i}->{'purity'}) * $phi_4;
      }
   } else {
      return (1 - $obj->{$i}->{'purity'}) * $phi_4;
   }
}

#  ===============
#  PIECEWISE PHI 1   computes the special switching function called "phi_1"
#  ===============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way
#
#  Implementation follows (completely general) Theorem 2 in supplemental
#  information of the paper rather than the slightly more concrete result in
#  "variant validation" notes on PP 94, but was checked case-by-case against
#  an example of p_T = 0.4, p_N = 0.3, mass fraction = 0.9 against the long-hand
#  results on pp 92-93 of those notes
#
#  TEST CASE WITH p_N = 0.3, p_T = 0.4, subclone mass fraction = 0.9
#
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 0, 0, T, 0, 0): numerator = 0.6 + 1.26
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 0, 0, N, 0, 0): numerator = 0.72 + 1.2
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 0, 0, T, 0, 1): numerator = 0.6 + 0.63
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 0, 1, N, 0, 0): numerator = 0.36 + 1.2
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 0, 0, T, 1, 0): numerator = 0.6 + 0.63
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 1, 0, N, 0, 0): numerator = 0.36 + 1.2
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 0, 0, T, 1, 1): numerator = 0.6 + 0
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 1, 1, N, 0, 0): numerator = 0 + 1.2
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 0, 1, T, 0, 0): numerator = 0.3 + 1.26
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 0, 0, N, 0, 1): numerator = 0.72 + 0.6
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 0, 1, T, 0, 1): numerator = 0.3 + 0.63
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 0, 1, N, 0, 1): numerator = 0.36 + 0.6
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 0, 1, T, 1, 0): numerator = 0.3 + 0.63
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 1, 0, N, 0, 1): numerator = 0.36 + 0.6
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 0, 1, T, 1, 1): numerator = 0.3 + 0
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 1, 1, N, 0, 1): numerator = 0 + 0.6
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 1, 0, T, 0, 0): numerator = 0.3 + 1.26
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 0, 0, N, 1, 0): numerator = 0.72 + 0.6
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 1, 0, T, 0, 1): numerator = 0.3 + 0.63
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 0, 1, N, 1, 0): numerator = 0.36 + 0.6
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 1, 0, T, 1, 0): numerator = 0.3 + 0.63
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 1, 0, N, 1, 0): numerator = 0.36 + 0.6
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 1, 0, T, 1, 1): numerator = 0.3 + 0
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 1, 1, N, 1, 0): numerator = 0 + 0.6
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 1, 1, T, 0, 0): numerator = 0 + 1.26
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 0, 0, N, 1, 1): numerator = 0.72 + 0
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 1, 1, T, 0, 1): numerator = 0 + 0.63
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 0, 1, N, 1, 1): numerator = 0.36 + 0
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 1, 1, T, 1, 0): numerator = 0 + 0.63
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 1, 0, N, 1, 1): numerator = 0.36 + 0
#  BERNOULLI NUM N   (i,j,k,l,m,n) = (N, 1, 1, T, 1, 1): numerator = 0 + 0
#  BERNOULLI NUM T   (i,j,k,l,m,n) = (T, 1, 1, N, 1, 1): numerator = 0 + 0

sub piecewise_phi_1 {
   my $obj = shift;
   my ($i, $j, $k) = @_;
   my $phi_4 = 1;
   $phi_4 = $obj->{'T'}->{'subclone_mass_fraction'} if $i eq 'T';
   if ($j == $k) {
      if ($j) {
         return 0;
      } else {
         return 2 * $obj->{$i}->{'purity'} * $phi_4;
      }
   } else {
      return $obj->{$i}->{'purity'} * $phi_4;
   }
}

#  =====================
#  READ_ERR_HARMONIC_AVG   error p-val harmonic avg starting from quality vals
#  =====================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  True "read qualities" are simply a semantic alternative to error p-values.
#  In an "expected value" context, a p-value is a measure of the *rate* of
#  error, e.g. 1 error expected in 10,000 bases. Given that this is a rate-based
#  parameter, the proper averaging technique is the harmonic method (see e.g.
#  moroney:1956:3:math:history)

sub read_err_harmonic_avg {
   my (@qvalues) = @_;
   my $num_vals = scalar @qvalues;
   my $total_denom = 0;
   foreach my $qval (@qvalues) {
      croak "quality value '$qval' must be an integer" unless $qval =~ /^\d+$/;
      croak "quality value '$qval' must be > 0" unless $qval > 0;
      $total_denom += 10**($qval/10);
   }
   return $num_vals/$total_denom;
}

#  =========
#  NUMERICAL   modfier to sort numerically
#  =========   ~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub _numerical_ {$a<=>$b}

#  ===============================
#  BIN READ COUNTS AND ERROR PROBS   perform tight binning process
#  ===============================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  This routine performs binning based on sorting the quality values (lowest
#  to highest) and then placing bin boundaries where differences between
#  consecutive quality values are greatest. (It actually orders the differences
#  too, from highest to lowest, so that, by proceeding in that order, the
#  fixed number of bins that are created capture the largest differences.)
#
#  Here's an example from the DNMT3A mutation in AML31 tumor
#
#  QVALS: 2 2 2 2  14 14 23 23 26 27 28 31 31 32 33 33 37 37 37 37 37 37 37 37..
#  DIFFS:  0 0 0 12  0  9  0  3  1  1  3  0  1  1  0  4  0  0  0  0  0  0  0....
#
#  note that the diffs are "in between" the quality values.
#
#    DIFF = 12
#      positions: 3, 
#    DIFF = 9
#      positions: 5, 
#    DIFF = 4
#      positions: 15, 
#    DIFF = 3
#      positions: 7, 10, 
#    DIFF = 1
#      positions: 8, 9, 12, 13, 23, 24, 
#    DIFF = 0
#      positions: 0, 1, 2, 4, 6, 11, 14, 16, 17, 18, 19, 20, 21, 22, 25, 26, 
#    MAXIMUM POSSIBLE BINS = 12
#
#    Assuming we requested 5 bins, then the routine would have made these
#    calculations
#
#    BIN DEMARCATIONS
#      3, 5, 7, 15, 27, 
#
#    BIN 1: 2 2 2 2                              err pval = 0.6309573 (4 reads)
#    BIN 2: 14 14                                err pval = 0.0398107 (2 reads)
#    BIN 3: 23 23                                err pval = 0.0050119 (2 reads)
#    BIN 4: 26 27 28 31 31 32 33 33              err pval = 0.0008313 (8 reads)
#    BIN 5: 37 37 37 37 37 37 37 37 38 39 39 39  err pval = 0.0001709 (12 reads)
#
#    Routine returns a list of the reads per bin and a list of the corresponding#    bin error probabilities the above example would return
#
#    ([0.6309573, 0.0398107, 0.0050119, 0.0008313, 0.0001709], [4, 2, 2, 8, 12])

sub bin_read_counts_and_error_probs {
   $DB::single = 1;
   my ($qvals, $num_bins) = @_;

#__LOCAL DEBUG FLAG (NOT PERMANENT)
#  my $debug = 0;

#__SORT THE LIST OF QUALITY VALUES FROM LOWEST TO HIGHEST
   @{$qvals} = sort _numerical_ @{$qvals};

#------------ DEBUG
#  if ($debug) {
#     print "\n\nQ VALS ARE: " . join (' ', @{$qvals}) . "\n";
#     print "NUMBER OF REQUESTED BINS = $num_bins\n\n";
#  }
#------------ DEBUG

#__CATALOG THE DIFFERENCES BETWEEN CONSECUTIVE QVALS AND THE PLACES THEY OCCUR
#  NOTING THAT RESULTING POSITIONS ARE INHERENTLY ORDERED FROM LOWEST TO HIGHEST
   my $diffs = {};
   for (my $i = 0; $i < $#{$qvals}; $i++) {
      my $diff = $qvals->[$i + 1] - $qvals->[$i];
      push @{$diffs->{$diff}}, $i;
   }

#__LET "0 DIFF" REPRESENT A SINGLE BIN (IN THAT THERE'S ONLY 1 ELEMENT IN ITS
#  LIST), E.G. IN CASE ALL QVALS IN THE LIST ARE IDENTICAL, ALSO, WE DO THIS IN
#  ALL CASES, EVEN WHEN NO $diffs->{'0'} IS DEFINED, SO THAT IT PROPERLY ADDS
#  AN EXTRA TALLY TO ACCOUNT FOR THE FACT THAT THESE ARE DEFINED "BETWEEN QVALS"
#  SO WE'LL INHERENTLY BE 1 SHORT ON THE MAX
#  $diffs->{'0'} = [0] if defined $diffs->{'0'};
   $diffs->{'0'} = [0];

#__FIND THE MAXIMUM NUMBER OF BINS THESE Q-VALUES *COULD* BE DIVIDED INTO
#  GIVEN THAT THERE MUST BE A DIFFERENCE OF AT LEAST "1" TO PERFORM A DIVISION
   my $max_bins = 0;
   foreach my $diff (keys %{$diffs}) {
      $max_bins += scalar @{$diffs->{$diff}};
   }

#__ACTUAL NUMBER OF BINS WE WILL SPLIT INTO IS THE LOWER OF THE TWO
   $num_bins = $max_bins if $max_bins < $num_bins;

#------------ DEBUG
#  if ($debug) {
#     foreach my $diff (reverse sort _numerical_ keys %{$diffs}) {
#        print "DIFF = $diff\n  positions: ";
#        foreach my $pos (@{$diffs->{$diff}}) {
#           print "$pos, ";
#        }
#        print "\n";
#     }
#     print "MAXIMUM POSSIBLE BINS = $max_bins\n";
#     print "NUMBER OF BINS WE WILL DIVIDE INTO = $num_bins\n";
#  }
#------------ DEBUG

#__IF WE'VE JUST ONE BIN THEN COMPUTE OVERALL ERROR PROBABILITY AND RETURN
#  my $bins = [];
   my ($bin_read_err_pvals, $bin_read_counts) = ([], []);
   if ($num_bins == 1) {
      my $err_prob = &read_err_harmonic_avg (@{$qvals});
#     push @{$bins}, [$err_prob, scalar @{$qvals}];
      push @{$bin_read_err_pvals}, $err_prob;
      push @{$bin_read_counts}, scalar @{$qvals};

#__ELSE FOLLOW-THROUGH WITH THE BINNING PROCESS
   } else {

   #__DISTILL LIST OF BIN DEMARCATIONS (WHICH ARE JUST LIST POSITIONS) ORDERED
   #  FROM MOST EXTREME DIFFERENCES TO LEAST
       my $bin_demarcs = [];
       foreach my $diff (reverse sort _numerical_ keys %{$diffs}) {
          foreach my $pos (@{$diffs->{$diff}}) {

          #__THIS IS THE DIVDER FOR THE NEXT BINNING DIVISION WE'LL MAKE
             push @{$bin_demarcs}, $pos;

          #__STOP WHEN NUMBER OF DIVIDERS IS 1 LESS THAN NUMBER OF NEEDED BINS
          #  (I.E. THERE'S INHERENTLY 1 FEWER DIVIDER THAN THE NUMBER OF BINS)
             goto DONE_WITH_DEMARCATIONS if scalar @{$bin_demarcs} == $num_bins - 1;
          }
       }
       DONE_WITH_DEMARCATIONS:

   #__NOW THAT WE HAVE THE DEMARCATIONS WE CAN ORDER THESE FROM LOWEST TO
   #  HIGHEST SO AS TO EASILY SLICE-UP THE LIST OF QUALITY VALUES INTO BINS
      @{$bin_demarcs} = sort _numerical_ @{$bin_demarcs};

   #__PUSH THE LAST INDEX ON SINCE THAT IS ALWAYS THE END BOUNDARY
      push @{$bin_demarcs}, $#{$qvals};

#------------ DEBUG
#  if ($debug) {
#     print "BIN DEMARCATIONS\n  ";
#     foreach my $pos (@{$bin_demarcs}) {
#        print "$pos, ";
#     }
#     print "\n";
#  }
#------------ DEBUG

   #__CREATE THE BINS
      my $lower_demarc = 0;
      foreach my $upper_demarc (@{$bin_demarcs}) {

      #__THIS SLICE IS THE CURRENT BIN
         my @bin = @{$qvals}[$lower_demarc .. $upper_demarc];

      #__ERROR PROBABILITY AND READ COUNT FOR THIS BIN
         my $err_prob = &read_err_harmonic_avg (@bin);
##       push @{$bins}, [$err_prob, scalar @bin];
         push @{$bin_read_err_pvals}, $err_prob;
         push @{$bin_read_counts}, scalar @bin;

#------------ DEBUG
#        if ($debug) {
#           print " THIS BIN: " . join (' ', @bin) . "\n";
#        }
#------------ DEBUG

      #__MOVE TO NEXT BIN
         $lower_demarc = $upper_demarc + 1;
      }
   }

#------------ DEBUG
#  if ($debug) {
#     print "error pvals and read counts\n";
#     foreach my $err_prob (@{$bin_read_err_pvals}){
#        print "   $err_prob\n";
#     }
#     foreach my $read_count (@{$bin_read_counts}){
#        print "   $read_count\n";
#     }
# #   foreach my $bin (@{$bins}){
# #      my ($err_prob, $num_reads) = @{$bin};
# #      print "   $err_prob ($num_reads reads)\n";
# #   }
#     print "\n";
#  }
#------------ DEBUG

#__RETURN RESULT
#  return $bins;
   return ($bin_read_err_pvals, $bin_read_counts);
}

######################################
#  PARSING AND TEXT-RELATED METHODS  #
######################################

#  ===============
#  TEXT_TO_INDECES   resolve components of verbose arg name for internal storage
#  ===============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  This is a private method, so we trust all arguments and do not check or
#  validate them here in any way

sub text_to_indeces {
   my ($text) = @_;
   my ($type, $key);

#__GET THE TYPE OF SAMPLE
   if ($text =~ /^normal/) {
      $type = "N";
   } else {
      $type = "T";
   }

#__GET THE NAME OF THE KEY UNDER WHICH TO STORE THIS INFORMATION
#
#  PROGRAMMING NOTE: NOW THAT ALL INPUT TAGS CORRESPOND EXACTLY TO THEIR DATA
#  STRUCTURE TAGS (BECAUSE WE NO LONGER HAVE THE QUALITY<->READ ERROR MAPPING)
#  YOU CAN SIMPLIFY THIS WHEN YOU HAVE A CHANCE

   if ($text =~ /purity/) {
      $key = "purity";
   } elsif ($text =~ /support_reads/) {
      $key = "support_reads";
   } elsif ($text =~ /num_q_bins/) {
      $key = "num_q_bins";
   } elsif ($text =~ /total_reads/) {
      $key = "total_reads";
   } elsif ($text =~ /read_qualities/) {
      $key = "read_qualities";
   } elsif ($text =~ /homozygous_variant_rate/) {
      $key = "homozygous_variant_rate";
   } elsif ($text =~ /heterozygous_variant_rate/) {
      $key = "heterozygous_variant_rate";
   } elsif ($text =~ /background_mutation_rate/) {
      $key = "background_mutation_rate";
   } elsif ($text =~ /subclone_mass_fraction/) {
      $key = "subclone_mass_fraction";
   }

#__RETURN RESULTS
   return ($type, $key);
}

################################################################################
##                                                                            ##
##             T R A I L I N G   P O D   D O C U M E N T A T I O N            ##
##                                                                            ##
################################################################################

################################################################################
##                                                                            ##
##                                -  E N D  -                                 ##
##                                                                            ##
################################################################################

1;
