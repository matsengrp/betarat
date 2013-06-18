
cdef extern from "rf_hyp2f1_core.h":
    double rf_hyp2f1_core(double a, int b, double c, double w, int mpf_prec)

def rf_hyp2f1(double a, int b, double c, double w, int mpf_prec=300):
    """ Computes the hypergeometric function 2F1(a, b; c; z) using a finitely terminating, recursively
    factored form of the function, as seen in Concrete Mathematics (Graham, Knuth, Patashnik - section 5.5). """
    assert b <= 0 and isinstance(b, (int, long))
    return rf_hyp2f1_core(a, b, c, w, mpf_prec)
