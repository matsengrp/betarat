"""
This module implements an adapted version of Simpson integration appropriate for computing quantiles. It
functions by integrating along the domain until it reaches the point x where the desired mass Q has been
obtained. It repeats this process for smaller and smaller step sizes until x has converged within some
tolerance.

Currently, this implemementation is limited in that the primary function it provides, simpson_quant_hp, only
works on a half plain. The underlying machinery could be made to work in general using a different transform,
but only the half plain is needed here.
"""

import betarat

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
    """MissingMass - if you expect that it should be possible to integrate up to a certain value (for
    example, if you have a probability distribution, it should be possible to construct an integral with any
    value 0 < x < 1), but after reaching the end of the domain this mass has not been obtained, this is an
    indication that your step size is too large, and that it needs to be greatly decreased."""
    pass

class MaxSumReached(Exception):
    """MaxSumReached - With probability distritions, your integral should never greatly exceed 1.0. Observing
    this sort of behaviour is a sign that you might be experiencing numerical difficulties in your
    computations (or that something else is wrong). This exception can be handled as you wish in such cases.
    """
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
    integration of pdf (all that is needed for the beta_rat distribution since D = [0, \inf)).
    
    max_sum: If the summation exceeds this value, the assumption will be that something has gone wrong (useful
    if you know, for example, that your integration should never exceed 1.0 - allows you to catch numerical
    instabilities, as found in compoutation of the hypergeometric).
    tolerance: Convergence criteria
    h_init: initial incremement for Simpson integration
    """

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
        if betarat.VERBOSE:
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
        if betarat.VERBOSE:
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



