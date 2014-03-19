#include "bvprob/Bassovac.hpp"
#include "bvprob/PBin.hpp"
#include "utility/Lut.hpp"
#include "bvprob/ResultFormatter.hpp"
#include "bvprob/Sample.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <iomanip>

using namespace std;

const double normHetVarRate = 0.001;
const double normHomVarRate = 0.002;
const double tumBgMutRate   = 0.003;

TEST(TestExpectedResult, againstPerl) {
    Lut::init();
    unsigned normalTotalReads = 40;
    unsigned normalSupportingReads = 39;
    unsigned tumorTotalReads = 93;
    unsigned tumorSupportingReads = 89;

    double normalPurity = 0.96;
    double tumorPurity = 0.8;
    double tumorMassFraction = 0.1;

    uint8_t normalReadQualities[] = {
        20, 20, 20, 22, 28, 23, 23, 23, 23, 25,
        25, 30, 30, 30, 30, 30, 30, 31, 32, 32,
        35, 35, 35, 35, 35, 35, 35, 35, 35, 35,
        38, 38, 38, 38, 38, 38, 38, 38, 38, 38,
    };
    int numNormalReadQualities = sizeof(normalReadQualities)/sizeof(normalReadQualities[0]);
    uint8_t tumorReadQualities[] = {
        20, 20, 20, 22, 28, 23, 23, 23, 23, 25,
        25, 30, 30, 30, 30, 30, 30, 31, 32, 32,
        25, 30, 30, 30, 30, 30, 30, 31, 32, 32,
        35, 35, 35, 35, 35, 35, 35, 35, 35, 35,
        35, 35, 35, 35, 35, 35, 35, 35, 35, 35,
        35, 35, 35, 35, 35, 35, 35, 35, 35, 35,
        38, 38, 38, 38, 38, 38, 38, 38, 38, 38,
        38, 38, 38, 38, 38, 38, 38, 38, 38, 38,
        38, 38, 38, 38, 38, 38, 38, 38, 38, 38,
        38, 38, 38
    };
    int numTumorReadQualities = sizeof(tumorReadQualities)/sizeof(tumorReadQualities[0]);

    sort(normalReadQualities, normalReadQualities + numNormalReadQualities);
    sort(tumorReadQualities, tumorReadQualities + numTumorReadQualities);

    Sample normal;
    normal.setValues(
        normalTotalReads,
        normalSupportingReads,
        0.0, // variant freq
        normalPurity,
        tumorMassFraction*(1.0-normalPurity),
        normalReadQualities,
        numNormalReadQualities,
        2);

    Sample tumor;
    tumor.setValues(
        tumorTotalReads,
        tumorSupportingReads,
        0.0, // variant freq
        tumorMassFraction*tumorPurity,
        1.0-tumorPurity,
        tumorReadQualities,
        numTumorReadQualities,
        2);
        
    double tumorBgMutationRate = 0.0000018;
    double normalHetVariantRate = 0.001;
    double normalHomVariantRate = 0.0005;
    Bassovac bv(normal, tumor, normalHetVariantRate, normalHomVariantRate, tumorBgMutationRate);
    cout << "      hom var: " << bv.homozygousVariantProbability() << "\n";
    cout << "      het var: " << bv.heterozygousVariantProbability() << "\n";
    cout << "      som var: " << bv.somaticVariantProbability() << "\n";
    cout << "          loh: " << bv.lossOfHeterozygosityProbability() << "\n";
    cout << "uninteresting: " << bv.nonNotableEventProbability() << "\n";
}

#if 0
TEST(TestExpectedResult, knownResult) {
    Lut::init();
    // Expected results as computed by perl version for this test:
    double expected[] = {
        3.59999336446043e-06, 9.40439137242974e-06, 2.2932130663145e-05, 5.26083185723004e-05,
        0.000114245490332158, 0.000236019982890073, 0.000465734766314079, 0.000880758048035665,
        0.00160067398875547, 0.00280201565033682, 0.00473340333367698, 0.00772798951632526,
        0.0122084814617112, 0.0186786207468407, 0.0276946666796615, 0.0398122863261926,
        0.055509300600671, 0.0750929547543096, 0.0986097062373996, 0.125781648779977,
        0.155991794709351, 0.188328594368296, 0.221681959332572, 0.2548668780169,
        0.286744171113004, 0.316313166910393, 0.342763862278863, 0.365489570436077,
        0.384070026388573, 0.398238037748338, 0.407841494154702, 0.412809258256655,
        0.413126031226559, 0.408818634788001, 0.399954406500588, 0.386651247052986,
        0.369097813381229, 0.347580986381765, 0.32251582736129, 0.294470929250176,
        0.264180081908436, 0.232530789929323, 0.200522868815962, 0.169196834009944,
        0.139541003329485, 0.112394885741348, 0.0883702798790751, 0.0678078687275174,
        0.0507769581312774, 0.037113776594739, 0.0264846069597731, 0.0184569115447004,
        0.0125641714204782, 0.00835580195585947, 0.00542937542384105, 0.00344667747918263,
        0.00213734990395305, 0.00129439292640894, 0.000765285426780347, 0.000441532673042524,
        0.000248466158748162, 0.000136298166122502, 7.28378159423546e-05, 3.78937108387204e-05,
        1.91775646377716e-05, 9.43369181968974e-06, 4.50660492818565e-06, 2.08874120793572e-06,
        9.3830228749506e-07, 4.08079455997443e-07, 1.71622147938886e-07, 6.97059330712936e-08,
        2.73041869433566e-08, 1.02990316590497e-08, 3.73473174725675e-09, 1.29970200347457e-09,
        4.33216281720641e-10, 1.38012582284323e-10, 4.19249657753089e-11, 1.21130336081988e-11,
        3.31917335407245e-12, 8.59888098533705e-13, 2.09881481828359e-13, 4.80765986142037e-14,
        1.02900968400884e-14, 2.04778080844137e-15, 3.76776459108271e-16, 6.36830170935639e-17,
        9.81462074806305e-18, 1.36731475722517e-18, 1.70439840541262e-19, 1.87794943173205e-20,
        1.80213030252369e-21, 1.47895047100436e-22, 1.01434831104643e-23, 5.64309001620183e-25,
        2.44633779684936e-26, 7.80994541841935e-28, 1.68668707082255e-29, 2.14071631951768e-31,
    };


    vector<double> normalQualities = {
        2, 2, 3, 2, 2, 12, 10, 10, 11, 13, 13, 10, 38, 38, 40, 36, 32, 35, 38
    };

    vector<double> tumorQualities(31, 38.0);
    for (uint32_t i = 0; i < 101; ++i)
        tumorQualities.push_back(4);
    
    const uint32_t normalTotal = normalQualities.size();;
    const uint32_t normalSupporting = normalTotal - 1;
    const double normalVariantFrequency = 0.5;
    const double normalPurity = 1.0;
    const uint32_t normalBins = 2;

    const uint32_t tumorTotal = tumorQualities.size();
    const uint32_t tumorSupporting = tumorTotal / 2 + 1;
    const double tumorPurity = 1.0;
    const double tumorMassFraction = 1.0;
    const uint32_t tumorBins = 2;

    const double normalHetVariantRate = 0.001;
    const double normalHomVariantRate = 0.0005;
    const double tumorBgMutationRate = 0.0000018;

    sort(normalQualities.begin(), normalQualities.end());
    sort(tumorQualities.begin(), tumorQualities.end());

    for (uint32_t i = 0; i < 100; ++i) {
        Sample normal;
        Sample tumor;

        const double tumorVariantFrequency = i/100.0;
        normal.setValues(normalTotal, normalSupporting, normalVariantFrequency, normalPurity,
            tumorMassFraction*(1 - normalPurity),
            &normalQualities[0], normalQualities.size(), normalBins);
        tumor.setValues(tumorTotal, tumorSupporting, tumorVariantFrequency, tumorMassFraction*tumorPurity,
            1 - tumorPurity,
            &tumorQualities[0], tumorQualities.size(), tumorBins);
        Bassovac bv(normal, tumor, normalHetVariantRate, normalHomVariantRate, tumorBgMutationRate);
        double result = bv.somaticVariantProbability();
        EXPECT_NEAR(expected[i], result, 1e-13)
             << "Failed with tumorVariantFrequency=" << tumorVariantFrequency;
    }
}

TEST(TestExpectedResult, funkyData) {
    Lut::init();
    // Expected results as computed by perl version for this test:
    double expected[] = {
        3.99999187710743e-06, 3.7282989035804e-06, 3.47257984402911e-06, 3.2320442702547e-06,
        3.00593335747991e-06, 2.79351892717615e-06, 2.59410250939906e-06, 2.40701442443324e-06,
        2.23161288354686e-06, 2.06728310865665e-06, 1.91343647070363e-06, 1.76950964653984e-06,
        1.63496379412632e-06, 1.50928374584207e-06, 1.39197721970425e-06, 1.28257404829923e-06,
        1.18062542522427e-06, 1.08570316883948e-06, 9.97399003129666e-07, 9.15323855475483e-07,
        8.39107171133377e-07, 7.68396244223652e-07, 7.02855565025977e-07, 6.42166183381568e-07,
        5.86025088001245e-07, 5.34144601478521e-07, 4.86251790806809e-07, 4.42087893199829e-07,
        4.01407757014232e-07, 3.6397929757343e-07, 3.29582967691594e-07, 2.98011242696755e-07,
        2.69068119751895e-07, 2.42568631272927e-07, 2.18338372242395e-07, 1.96213041217746e-07,
        1.76037994832973e-07, 1.5766781559244e-07, 1.40965892755656e-07, 1.25804016111767e-07,
        1.12061982442524e-07, 9.96272144724553e-08, 8.83943921049846e-08, 7.82650957431997e-08,
        6.91474614939919e-08, 6.09558480542613e-08, 5.3610515077883e-08, 4.70373128221216e-08,
        4.1167382872173e-08, 3.59368697425113e-08, 3.12866431537091e-08, 2.71620307833986e-08,
        2.35125612900353e-08, 2.02917174081233e-08, 1.74566989135577e-08, 1.49681952577376e-08,
        1.27901676690999e-08, 1.08896405207228e-08, 9.23650176264544e-09, 7.80331221754871e-09,
        6.56512353844128e-09, 5.49930462699265e-09, 4.58537631115494e-09, 3.80485408071365e-09,
        3.14109867940699e-09, 2.57917435225272e-09, 2.10571454672067e-09, 1.7087948663888e-09,
        1.3778130757202e-09, 1.10337595459797e-09, 8.77192801254618e-10, 6.91975382232583e-10,
        5.41344128011913e-10, 4.19740372941288e-10, 3.22344438108234e-10, 2.44999355784324e-10,
        1.84140034081072e-10, 1.36727660452123e-10, 1.00189142677305e-10, 7.23613859640417e-11,
        5.14402048015976e-11, 3.59336682036005e-11, 2.46196769742483e-11, 1.65075716336061e-11,
        1.08035696373793e-11, 6.87983052653798e-12, 4.24694764216539e-12, 2.52966504089409e-12,
        1.44561824629307e-12, 7.86897471564688e-13, 4.0423194891028e-13, 1.93594000921118e-13,
        8.50214893515244e-14, 3.34572766301882e-14, 1.14043354656703e-14, 3.19514190190191e-15,
        6.73987179453859e-16, 9.08437583067776e-17, 5.42085785154113e-18, 4.48698391733568e-20
    };    

    string qnormal = "99;BBBBBBBBBBBB";
    string qtumor = "BBBBBBB";

    vector<double> normalQualities;
    vector<double> tumorQualities;

    for (uint32_t i = 0; i < qnormal.size(); ++i)
        normalQualities.push_back(qnormal[i]-33);

    for (uint32_t i = 0; i < qtumor.size(); ++i)
        tumorQualities.push_back(qtumor[i]-33);

    const uint32_t normalTotal = normalQualities.size();;
    const uint32_t normalSupporting = normalTotal;
    const double normalVariantFrequency = 0.5;
    const double normalPurity = 1.0;
    const uint32_t normalBins = 2;

    const uint32_t tumorTotal = tumorQualities.size();
    const uint32_t tumorSupporting = tumorTotal;
    const double tumorPurity = 1.0;
    const double tumorMassFraction = 1.0;
    const uint32_t tumorBins = 2;

    const double normalHetVariantRate = 0.001;
    const double normalHomVariantRate = 0.0005;
    const double tumorBgMutationRate = 2e-6;

    sort(normalQualities.begin(), normalQualities.end());
    sort(tumorQualities.begin(), tumorQualities.end());


    ResultFormatter f(&cout, false, 16);
    for (uint32_t i = 0; i < 100; ++i) {
        Sample normal;
        Sample tumor;
        const double tumorVariantFrequency = i/100.0;
        normal.setValues(normalTotal, normalSupporting, normalVariantFrequency, normalPurity,
            tumorMassFraction*(1-normalPurity),
            &normalQualities[0], normalQualities.size(), normalBins);
        tumor.setValues(tumorTotal, tumorSupporting, tumorVariantFrequency, tumorMassFraction*tumorPurity,
            1-tumorPurity,
            &tumorQualities[0], tumorQualities.size(), tumorBins);
        Bassovac bv(normal, tumor, normalHetVariantRate, normalHomVariantRate, tumorBgMutationRate);
        double result = bv.somaticVariantProbability();
        EXPECT_NEAR(expected[i], result, 1e-13)
             << "Failed with tumorVariantFrequency=" << tumorVariantFrequency;
    }
}
#endif
