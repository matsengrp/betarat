#!/usr/bin/env python

from __future__ import division
from scipy import integrate, optimize
from scipy.special import beta
from simpson_quant import simpson_quant_hp, MaxSumReached
import mpmath
import pylab


# Global
VERBOSE = False



def plot_fn(f, a, b, points=100, label=None):
    """ Nice utility function for plotting functions using pylab."""
    inc = float(b-a) / points
    xs = pylab.arange(a, b, inc)
    ys = [f(x) for x in xs]
    if not label:
        label = str(f)
    pylab.plot(xs, ys, label=label)
    pylab.legend()


class BetaRat(object):
    """ Bata Ratio class - instantiate for an object on which to compute pdf etc """
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

    def h2f1_l(self, w, hypf=mpmath.hyp2f1):
        """ Left hand side of the function. """
        return hypf(self.a1 + self.a2, 1 - self.b1, self.a1 + self.a2 + self.b2, w)

    def h2f1_r(self, w, hypf=mpmath.hyp2f1):
        """ Right hand side of the function. """
        return hypf(self.a1 + self.a2, 1 - self.b2, self.a1 + self.a2 + self.b1, 1.0/w)

    def invert(self):
        """ Return inverted version of this beta ratio... """
        return BetaRat(self.a2, self.a1, self.b2, self.b1, prior=(0,0))

    def invert_ppf_if_needed(orig_ppf):
        """ This evaluates the integral using the inverse ratio X2/x1 if we would naively estimate that ratio
        to be > 1. """
        def wrapper(self, q, **kw_args):
            if self.inverted:
                inv_ppf = orig_ppf(self.inverted, 1-q, **kw_args)
                if inv_ppf == 'NA':
                    return inv_ppf
                return 1 / inv_ppf
            else:
                return orig_ppf(self, q, **kw_args)
        return wrapper

    def __repr__(self):
        rep = "BetaRat({}, {}, {}, {})".format(self.a1, self.a2, self.b1, self.b2)
        if self.inverted:
            rep += "<inv>"
        return rep

    def pdfs(self, ws):
        """ PDF of list - basically so we can apply scipy.integrate for CDF. """
        return [self.pdf(w) for w in ws]

    def cdf(self, w):
        """ Cumulative density function. """
        return integrate.quadrature(self.pdfs, 0, w)[0]

    def pdf(self, w, hypf=mpmath.hyp2f1):
        """ Probability density function. """
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

    def map(self):
        """ Maximum a posteriori: compute the mode of the posterior distribution. """
        return optimize.brent(lambda w: - self.pdf(w))
    
    @invert_ppf_if_needed
    def ppf(self, q, max_sum=10, **kw_args):
        """ Quantile function. Keyword args are the same as for simpson_quant_hp. """
        try:
            return simpson_quant_hp(self.pdf, q, max_sum=max_sum, **kw_args)
        except MaxSumReached:
            return "NA"

    def stupid_ppf(self, q, **kw_args):
        """ Quantile function. Keyword args are the same as for simpson_quant_hp. """
        def find_outer(b=1):
            if self.cdf(b) > q:
                return b
            else:
                return find_outer(b*2)
        return optimize.brenth(lambda x: self.cdf(x) - q, 0, find_outer())


    def plot_pdf(self, a=0, b=5, points=100):
        """ Convenience function for plotting a BetaRatio using pylab over the range (a, b) with points points"""
        plot_fn(self.pdf, a, b, points=points, label=self.__repr__())


