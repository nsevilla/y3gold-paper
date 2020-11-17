### needs python 2.7
from scipy.integrate import quad
import math
import numpy as np

class efficiencyError:
    """
    Calculates the error on an 'efficiency'
    Check out http://home.fnal.gov/~paterno/images/effic.pdf
    and http://nbviewer.ipython.org/github/nsevilla/exploratory_notebook/blob/master/compute_sgsep.ipynb
    This class basically copied from there.
    call like:
    k = 1568
    N = 1575
    confidence_level = 0.678
    print 'N=', N, 'k=', k
    print efficiencyError(N, k, confidence_level).calculate()

    k = number of hits
    N = number of events
    confidence_level = required confidence interval [0,1]

    class author: Ben Hoyle
    Date: 30 Oct 2015
    """

    def __init__(self, N, k, confidence_level):
        self.N = N
        self.k = k
        self.confidence_level = confidence_level

    def _pdf(self, eps, N, k):
        lnorm = math.lgamma(N + 2) - math.lgamma(k + 1) - math.lgamma(N - k + 1)
        val = lnorm + k * np.log(eps) + (N - k) * np.log(1 - eps)
        return np.exp(val)

    def calculate(self):
        N_, k_, confidence_level_ = self.N, self.k, self.confidence_level
        #print float(k_),float(N_),float(k_) / float(N_)
        eps = float(k_) / float(N_)
        #integral for error estimate
        interval_list = []
        # not much difference using numpy and this way I can use append
        step = 0.001 #0.0001
        xa = np.arange(eps-0.05, np.minimum(eps+0.05,0.999999), step)
        counter = 0
        for alpha in xa:
            xb = np.arange(alpha, np.minimum(eps+0.05,0.999999), step)
            for beta in xb:
                I = quad(self._pdf, alpha, beta, args=(N_, k_))
                counter = counter + 1
                if I[0] > confidence_level_:
                    interval_list.append([beta - alpha, alpha, beta])
                    break

        # should match confidence level
        minimum_integral = np.amin(interval_list, axis=0)[0]
        minimum_index = np.argmin(interval_list, axis=0)[0]

        return eps, interval_list[minimum_index][1], interval_list[minimum_index][2]


def test_eff():
    """
    unit test for case give in notes
    """
    k = 1568
    N = 1575
    confidence_level = 0.678
    print('N=', N, 'k=', k)
    res = efficiencyError(N, k, confidence_level).calculate()
    np.testing.assert_equal( (np.abs(res[0] - 0.9955555555555555) < 0.001 ) * (np.abs(res[1] - 0.0013444444444392634) < 0.001) * (np.abs( res[2] - 0.0020555555555603622) < 0.001), True)

#test_eff()
