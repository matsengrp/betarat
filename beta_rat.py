#!/usr/bin/env python

from __future__ import division
import scipy
import scipy.special


VERBOSE = False


__version__ = '0.0.1'

B = scipy.special.beta
h2f1 = scipy.special.hyp2f1
gamma = scipy.special.gamma


def bigA(a1, a2, b1, b2):
    return 

class BetaRat(object):
    """ Bata Ratio class - instantiate for an object on which to compute pdf etc """
    def __init__(self, a1, a2, b1, b2):
        """ Arguments a1, a2, b1, b2 are the beta parameters for distributions X1 and X2 """
        self.a1, self.a2, self.b1, self.b2 = a1, a2, b1, b2
        self.A = B(a1, b1) * B(a2, b2)
        self.Blt = B(a1 + a2, b2) 
        self.Bgt = B(a1 + a2, b1) 

    def pdf(self, w):
        """ Probability density function at value w """
        if w <= 1:
            return self.Blt * (w ** (self.a1 -1)) * h2f1(self.a1 + self.a2, 1 - self.b1, self.a1 + self.a2 +
                    self.b2, w) / self.A
        else:
            return self.Bgt * (w ** -(1+self.a2)) * h2f1(self.a1 + self.a2, 1 - self.b2, self.a1 + self.a2 +
                    self.b1, 1/w) / self.A

    def ppf(self, q, **kw_args):
        return simpson_quant_hp(self.pdf, q, **kw_args)



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



def determine_eval_set(i, n):
    """ Awesome Aaron Magic Sauce - operation that find the least significant digit in the
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


def simpson_quant_hp(f, q, a=0, tolerance=5e-4, h_init=0.005):
    """ Simpson Half Plane Quantiles - uses a half infinite interval transform to give finite bounds for
    integration of pdf (all that is needed for the beta_rat distribution since D = [0, \inf)). This function
    could be modified to use full-infinite integral transform for general use case scenario. """

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
            print "Depth level:", depth
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
        def within(s1, s2):
            return abs(1 - s1 / s2) > tolerance

        if len(sums) < 3:
            return True
        return within(sums[-1], sums[-2]) and within(sums[-2], sums[-3])

    def final_iterpolate():
        pass

    eval_sets = []
    x = None
    sums = []

    h = h_init
    while still_converging(sums):
        try:
            cur_sum, x = next_level(eval_sets, h)
        except MissingMass:
            print "Missing mass - decreasing h_init"
            return simpson_quant_hp(h, q, a=a, tolerance=tolerance, h_init=h_init*(0.1))
        sums.append(cur_sum * (h/3))
        h /= 2.0
    
    # This is required to transform the limit of integration back into the original unbounded coordinates
    return g(x)



def normal_test():
    global VERBOSE
    VERBOSE = True

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('q', default=0.25, type=float)
    parser.add_argument('-t', default=1e-5, type=float)
    args = parser.parse_args()
    from scipy.stats import norm
    from time import time
    t1 = time()
    solution = simpson_quant_hp(norm.pdf, args.q, tolerance=args.t)
    t2 = time()
    print "Actual solition: ", norm.ppf(0.5 + args.q)
    print "My solution:     ", solution
    print "\nTime: ", t2 - t1, "s\n"



def cli():
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""
    Beta Rat (v {})!!! Awesome contingency table stats!

    For table

                X1   X2
         succ   a    b
         fail   c    d

    Evaluates quantile q of X1/X2.""".format(__version__)
            )

    parser.add_argument('a', type=int, help='# sucesses in X1')
    parser.add_argument('b', type=int, help='# sucesses in X2')
    parser.add_argument('c', type=int, help='# failures in X1')
    parser.add_argument('d', type=int, help='# failures in X2')

    parser.add_argument('q', type=float, help='quantile', default=0.05)

    parser.add_argument('--prior-succ', type=float, help='prior on beta(a*, b*)', default=1)
    parser.add_argument('--prior-fail', type=float, help='prior on beta(a*, b*)', default=1)

    #parser.add_argument('--prob-diff', type=float, default=0.0, help='Test that P1 - P2 > prob-diff')

    args = parser.parse_args()

    a, b, c, d = [getattr(args, x) for x in ['a', 'b', 'c', 'd']]

    beta_rat = BetaRat(a + args.prior_succ, b + args.prior_succ, c + args.prior_fail, d + args.prior_fail)
    result = beta_rat.ppf(args.q)

    print "\nP ( P1 / P2 < q ) =", result, '\n'




if __name__ == '__main__':
    normal_test()




# Here there be dragons
def h3F2(a1, a2, c, b1, d, z):
    g_part = gamma(d) / (gamma(c) * gamma(d-c))
    def integrand(t):
        t_part = t ** (c-1) * (1 - t) ** (d-c-1)
        return t_part * h2f1(a1, a2, b1, t*z)
    return g_part * scipy.integrate.quad(integrand, 0, 1)

def cdf(a1, a2, b1, b2, w):
    pass


