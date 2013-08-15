
#!/usr/bin/env python
from __future__ import division
from scipy import integrate, optimize
from scipy.special import beta
from simpson_quant import simpson_quant_hp, MaxSumReached
from utils import plot_fn
import mpmath
import math


# Global
VERBOSE = False


def invert_ppf_if_needed(orig_ppf):
    """ This decorator clean up the logic of evaluating the ppf using the inverse ratio x2/x1 when we would
    estimate that it is > 1. This speeds up some computations using the simpson_quant_hp. """
    def wrapper(self, q, **kw_args):
        if self.inverted:
            inv_ppf = orig_ppf(self.inverted, 1-q, **kw_args)
            if inv_ppf == 'NA':
                return inv_ppf
            return 1 / inv_ppf
        else:
            return orig_ppf(self, q, **kw_args)
    return wrapper


class BetaRat(object):
    """ BetaRat class - representation of the Beta Ratio distribution. """
    def __init__(self, a1, a2, b1, b2, no_inverting=False, prior=(1.0, 1.0)):
        """ Arguments a1, a2, b1, b2 are the beta parameters for distributions X1 and X2. Now can specify a
        prior which will be tacked on to (a1, b1) and (a2, b2). """
        # See invert_pdf_if_needed for explanation of this madness
        self.a1, self.a2, self.b1, self.b2 = a1 + prior[0], a2 + prior[0], b1 + prior[1], b2 + prior[1]
        if a1*b2 > a2*b1 and not no_inverting:
            self.inverted = self.invert()
        else:
            self.inverted = False

        self.A = beta(self.a1, self.b1) * beta(self.a2, self.b2)
        self.Blt = beta(self.a1 + self.a2, self.b2) 
        self.Bgt = beta(self.a1 + self.a2, self.b1)

    def __repr__(self):
        rep = "BetaRat({}, {}, {}, {})".format(self.a1, self.a2, self.b1, self.b2)
        if self.inverted:
            rep += "<inv>"
        return rep

    def invert(self):
        """ Return inverted version of this beta ratio... """
        return BetaRat(self.a2, self.a1, self.b2, self.b1, prior=(0,0))

    def h2f1_l(self, w, hypf=mpmath.hyp2f1):
        """ Application of the hypergeometric function for computing the pdf varies depending on whether w < 1
        or not. Left hand side of the function. """
        return hypf(self.a1 + self.a2, 1 - self.b1, self.a1 + self.a2 + self.b2, w)

    def h2f1_r(self, w, hypf=mpmath.hyp2f1):
        """ Right hand side of the function. """
        return hypf(self.a1 + self.a2, 1 - self.b2, self.a1 + self.a2 + self.b1, 1.0/w)

    def pdf(self, w, hypf=mpmath.hyp2f1):
        """
        Probability Density Function.
        hypf argument points to the function to be used for computation of the H2F1 function. This defaults to
        `mpmath.hyp2f1`. There is also an implementation which relies on the recursively factored form of the
        function, accessible as `betarat.rf_hyp2f1`. This function can be faster in cases of moderate failure
        counts, but not generally, and only works when the b_{1,2} are integers. As such, this function only
        works with priors Beta(a', b') where a' is an integer.
        """
        if w == 0:
            return 0
        elif w <= 1:
            return (self.Blt *
                (w ** (self.a1 -1)) * 
                self.h2f1_l(w, hypf=hypf) /
                self.A)

        else:
            return (self.Bgt *
                (w ** -(1+self.a2)) *
                self.h2f1_r(w, hypf=hypf) /
                self.A)

    def pdfs(self, ws):
        """ PDF of list - enables easier application of scipy.integrate for CDF. """
        return [self.pdf(w) for w in ws]

    def cdf(self, w):
        """ Cumulative Density Function. """
        return integrate.quadrature(self.pdfs, 0, w)[0]

    def map(self):
        """
        Maximum A Posteriori: compute the mode of the posterior distribution. This can be more robust than the
        median (ppf(0.5)) as an estimate of the most likely ratio in cases where one of the counts of
        negatives is zero, as this situation can lead to a long tail and very high median values. In
        situations where this is not a problem, the median and MAP tend to agree fairly well.
        """
        return optimize.brent(lambda w: - self.pdf(w))

    def lt_pdf(self, t):
        "Log transformed posterior density function"
        return math.exp(t) * self.pdf(math.exp(t))

    def lt_map(self):
        """Exponential of the MAP of the log transformed PDF: This estimator more equally treats relative
        probability ratios less more equally than those greater than oen. However, this comes at the expense
        of a higher mean square error as an estimator."""
        return math.exp(optimize.brent(lambda t: - self.lt_pdf(t)))
    
    @invert_ppf_if_needed
    def simpson_ppf(self, q, max_sum=10, **kw_args):
        """
        Quantile function (AKA Percentile Point Function).
        This implementation uses the simpson_quantile method in betarat.simpson_quant.
        Keyword args are the same as for simpson_quant_hp.
        """
        try:
            return simpson_quant_hp(self.pdf, q, max_sum=max_sum, **kw_args)
        except MaxSumReached:
            return "NA"

    def optim_ppf(self, q, maxiter=75, **kw_args):
        """
        Quantile function (AKA Percentile Point Function).
        This implementation uses scipy.optimize to solve for CDF(x) = Q. While in some cases, the
        simpson_quantile implementation ends up being faster, there are also cases where it is WAY slower, and
        in general it seems to be less accurate as well. As such, we generally advocate using this
        implementation. Keyword args are passed along to optimize.brenth.
        """
        def find_outer(b=1):
            if self.cdf(b) > q:
                return b
            else:
                return find_outer(b*2)
        return optimize.brenth(lambda x: self.cdf(x) - q, 0, find_outer(), maxiter=maxiter, **kw_args)

    def ppf(self, q, method="optim", **kw_args):
        """
        Quantile function (AKA Percentile Point Function).
        This convenience function wraps either optim_ppf or simpson ppf, as specied by `method` argument.
        Any keyword args are passed along to the appropriate function.
        """
        if method == 'optim':
            fn = self.optim_ppf
        elif method == 'simpson':
            fn = self.simpson_ppf
        else:
            raise ValueError, "Unrecognized method. Valid method options are optim and simpson."
        return fn(q, **kw_args)

    def plot_pdf(self, a=0, b=5, points=100):
        """ Convenience function for plotting a BetaRat using pylab over the range (a, b) with points
        points. As one might expect, this requires the pylab library, and it is suggested that the ipython
        library be used for this. For example, running `ipython --pylab` and importing this module should
        allow you to use this function. """
        plot_fn(self.pdf, a, b, points=points, label=self.__repr__())

    def plot_ltpdf(self, a=-5, b=5, points=100):
        """ Convenience function for plotting a BetaRat using pylab over the range (a, b) with points
        points. As one might expect, this requires the pylab library, and it is suggested that the ipython
        library be used for this. For example, running `ipython --pylab` and importing this module should
        allow you to use this function. """
        plot_fn(self.lt_pdf, a, b, points=points, label=self.__repr__())

