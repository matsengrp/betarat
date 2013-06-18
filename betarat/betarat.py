#!/usr/bin/env python

from __future__ import division
from scipy import integrate, optimize
from scipy.special import beta
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



class EvalSet(object):
    """ This class keeps track of evaluations of a function under question and enables us to fetch the x
    values from the evaluations arose """
    def __init__(self, f, a, h):
        """ f = function to evaluate; a = initial f input; h = incremement between succesive value of x """
        self.f = f
        self.a, self.h = a, h
        self.evaluations = []

    def __getitem__(self, index):
        """ Fetch correponding index item, or create it if doesn't exist (just appends, so don't want to use
        this if the index is not len(self.evaluations) """
        try:
            return self.evaluations[index]
        except IndexError:
            y = self.f(self.a + self.h*index)
            self.evaluations.append(y)
            return y

    def x_at(self, index):
        """ Computes the value of x which produced the f evaluation at a given index in evaluations list """
        return self.a + index * self.h


class MissingMass(Exception):
    pass

class MaxSumReached(Exception):
    pass


def determine_eval_set(i, n):
    """ Awesome Aaron G. Magic Sauce - operation that find the least significant digit in the
    bitwise and operation, as this will be the index of the eval set that we will need when n
    is the depth of our eval_set_list (init dept==0)."""
    i %= 2 ** n
    for e in xrange(n, -1, -1):
        if i & 1:
            return e
        i >>= 1
    return 0

def inner_index(i, n, j):
    """ The 'inner_index', returned from this function for figuring out from the outer index (i), the depth
    (n) and the eval_set_index (j), is the current index of the value in the corresponding eval_set's
    evaluation list. The conditional on the exponent is to account for the distance between the first level
    being the same as the second """
    lower_exp = (n - j) if j == 0 else (n - j) + 1
    return int((i - 2 ** (n - j)) / 2 ** lower_exp)

def navigate(eval_sets, max_iter):
    """ This is how we iterate through the eval sets and their evaluations """
    depth = len(eval_sets) - 1
    for outer_index in xrange(1, max_iter):
        es_index = determine_eval_set(outer_index, depth)
        yield outer_index, eval_sets[es_index], inner_index(outer_index, depth, es_index)


def simpson_quant_hp(f, q, a=0, tolerance=5e-4, h_init=0.005, max_sum=None):
    """ Simpson Half Plane Quantiles - uses a half infinite interval transform to give finite bounds for
    integration of pdf (all that is needed for the beta_rat distribution since D = [0, \inf)). This function
    could be modified to use full-infinite integral transform for general use case scenario. 
    
    max_sum: If the summation exceeds this value, the assumption will be that something has gone terribly
    wrong. """

    # XXX - Should still really put something in as we were doing that flips the function around the x-axis. Wouldn't
    # be hard and would likely save us a lot of hassle
    
    g = lambda t: a + t / (1 - t)
    g_ = lambda t: 1 / (1 - t) ** 2
    f_ = lambda t: f(g(t)) * g_(t)

    def next_level(eval_sets, h):
        # This handles the next eval_set down
        stop_at = q * (3/h)
        cur_sum = f_(0)
        # Depth is 0-based, so we set this before we've added the latest EvalSet
        depth = len(eval_sets)
        if VERBOSE:
            print "Depth level:", depth,
        es_h = h if depth == 0 else h * 2.0
        es_a = h / 2 ** (depth)
        eval_sets.append(EvalSet(f_, es_a, es_h))
        # Don't want to have a ZeroDivisionError
        max_iter = int(1.0/h)
        for outer_index, eval_set, inner_index in navigate(eval_sets, max_iter):
            val = eval_set[inner_index]
            odd = outer_index % 2
            if odd:
                cur_sum += val * 4
            else:
                if cur_sum + val > stop_at:
                    return cur_sum + val, eval_set.x_at(inner_index)
                cur_sum += val * 2
        # The asumption here is that if we haven't found the mass we are looking for at this point, then we
        # are probably missing mass that is all lumped up by itself somewhere, and need to greatly decrease
        # our incremement size
        raise MissingMass

    def still_converging(sums):
        def not_within_tol(s1, s2):
            return abs(1 - s1 / s2) > tolerance

        if len(sums) < 3:
            return True
        return not_within_tol(sums[-1], sums[-2]) and not_within_tol(sums[-2], sums[-3])

    eval_sets = []
    x = None
    sums = []

    h = h_init
    while still_converging(sums):
        try:
            cur_sum, x = next_level(eval_sets, h)
            cur_sum *= (h/3)
        except MissingMass:
            print "Missing mass - decreasing h_init"
            return simpson_quant_hp(f, q, a=a, max_sum=max_sum, tolerance=tolerance, h_init=h_init*(0.1))
        if VERBOSE:
            print "   Val:", cur_sum,
            print "   x:", x,
            print "   g(x):", g(x)
        if max_sum:
            if cur_sum > max_sum:
                raise MaxSumReached
        sums.append(cur_sum)
        h /= 2.0

    # This is required to transform the limit of integration back into the original unbounded coordinates
    return g(x)



