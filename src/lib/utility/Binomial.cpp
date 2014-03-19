#include "Binomial.hpp"

#include "Lut.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace Binomial {
    using namespace std;

    double pdf(double p, int n, int k) {
        double logbinc;
        if (n == k || k == 0)
            logbinc = 0;
        else
            logbinc = Lut::lgamma(n+1) - Lut::lgamma(k+1) - Lut::lgamma(n-k+1);

        return exp(logbinc + log(p)*k + log1p(-p)*(n-k));
    }

    // Returns P(Z=k) where
    //  Z = X + Y
    //  X ~ Binomial(p1, n1)
    //  Y ~ Binomial(p2, n2)
    //  X indep. Y
    //
    // See the test case in test/lib/bvprob/TestBinomial.cpp for a worked
    // example.
    double pdfConvolve2(double p1, double p2, int n1, int n2, int k) {
        double lp1 = log(p1);
        double lq1 = log1p(-p1);
        double lp2 = log(p2);
        double lq2 = log1p(-p2);
        int begin = std::max(0, k-n2);
        int limit = std::min(k, n1);

        // bcTop is log(n1! * n2!)
        double bcTop = Lut::lgamma(n1+1) + Lut::lgamma(n2+1);

        // Log gamma lookup arrays
        double const* lgamma_up1 = Lut::lgamma_arr(begin + 1);
        double const* lgamma_down1 = Lut::lgamma_arr(n1 - begin + 1);
        double const* lgamma_down2 = Lut::lgamma_arr(k - begin + 1);
        double const* lgamma_up2 = Lut::lgamma_arr(n2 - (k - begin) + 1);

        int x = 0;
        double rv(0.0);
        for (int i = begin; i <= limit;
              ++i, ++x
            , ++lgamma_up1, --lgamma_down1
            , --lgamma_down2, ++lgamma_up2
            )
        {
            int i2 = k - i;

            assert(i + i2 == k);
            assert(i <= n1);
            assert(i2 <= n2);

            // bc = log(n1 choose i * n2 choose i2)
            double bc = bcTop
                - *lgamma_up1     // lgamma(i + 1)
                - *lgamma_down1   // lgamma(n1 - i+1)
                - *lgamma_down2   // lgamma(i2 + 1)
                - *lgamma_up2     // lgamma(n2 - i2 + 1)
                ;

            // log(dbinom(i, n1, p1) * dbinom(i2, n2, p2))
            double logValue = bc +
                i*lp1 + (n1-i)*lq1 +
                i2*lp2 + (n2-i2)*lq2;

            rv += std::exp(logValue);
        }

        return rv;
    }
}
