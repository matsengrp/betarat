#!/usr/bin/env python

from __future__ import division
import scipy
import scipy.special


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
        if q <= 0.5:
            return simpson_quant(self.pdf, 0, q, **kw_args)
        else:
            # In this case, we use a transform to alter the problem so that we can compute the area under the
            # curve from the left. This helps us avoid issues with missing the mass really early on when
            # curves are supremely steep
            left_fun = lambda t: self.pdf((t-1)/t) /t ** 2
            a = simpson_quant(left_fun, 0, 1-q, remove_first_point=True, **kw_args)
            return (a - 1) / a



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


class MaxIterationsExceeded(Exception):
    pass

class MaxDepthExceeded(Exception):
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

def navigate(eval_sets):
    """ This is how we iterate through the eval sets and their evaluations """
    depth = len(eval_sets) - 1
    outer_index = 1
    while True:
        es_index = determine_eval_set(outer_index, depth)
        yield outer_index, eval_sets[es_index], inner_index(outer_index, depth, es_index)
        outer_index += 1


def simpson_quant(f, a, q, tolerance=1e-5, h_init=0.001, max_init_iter=1e5, remove_first_point=False):

    def next_level(eval_sets, h):
        stop_at = q * (3/h)
        if remove_first_point:
            cur_sum = 0
        else:
            cur_sum = f(a)
        es_h = h if h == h_init else h * 2.0
        es_a = h / 2 ** (len(eval_sets) - 1)
        eval_sets.append(EvalSet(f, es_a, es_h))
        x = None
        for outer_index, eval_set, inner_index in navigate(eval_sets):
            val = eval_set[inner_index]
            outer_index %= 2
            if outer_index:
                cur_sum += val * 4
            else:
                if cur_sum + val > stop_at or (outer_index > max_init_iter):
                    x = eval_set.x_at(inner_index)
                    break
                cur_sum += val * 2
            outer_index += 1
        return cur_sum, x

    def still_converging(cur_sum, last_sum):
        if not last_sum:
            return True
        return abs(1 - cur_sum / last_sum) > tolerance

    eval_sets = []
    last_sum = None
    cur_sum = None
    x = None

    h = h_init
    level = 0
    while still_converging(cur_sum, last_sum):
        level += 1
        last_sum = cur_sum
        cur_sum, x = next_level(eval_sets, h)
        cur_sum *= (h/3)
        h /= 2.0

    return x



def simple_simpson_quant(f, a, q, h=10e-6):
    S = f(a) + 4*f(a+h)
    stop_at = q * (3/h)

    i = 2
    while S + f(a + i*h) < stop_at:
        # even
        S += 2 * f(a + i*h)
        i += 1
        # odd
        S += 4 * f(a + i*h)
        i += 1

    return a + i*h


def normal_test():
    from scipy.stats import norm
    print "Actual solition: ", norm.ppf(0.75)
    print "My solution:     ", simpson_quant(norm.pdf, 0.0, 0.25)



def cli():
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""
    Beta ratio (v {}) for contingency table

                X1   X2
         succ   a    b
         fail   c    d

    Evaluates probability that P1 / P2 < q.""".format(__version__)
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
    cli()




# Here there be dragons
def h3F2(a1, a2, c, b1, d, z):
    g_part = gamma(d) / (gamma(c) * gamma(d-c))
    def integrand(t):
        t_part = t ** (c-1) * (1 - t) ** (d-c-1)
        return t_part * h2f1(a1, a2, b1, t*z)
    return g_part * scipy.integrate.quad(integrand, 0, 1)

def cdf(a1, a2, b1, b2, w):
    pass


