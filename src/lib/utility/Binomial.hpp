#pragma once

namespace Binomial {
    double pdf(double p, int n, int k);
    double pdfConvolve2(double p1, double p2, int n1, int n2, int k);
}
