
from betarat import BetaRat, VERBOSE, defaults
from betarat.version import __version__
import argparse
import time 


# Table count string for CLI usage message
table_string = """
                X1   X2
         succ   a    b
         fail   c    d
         """

def setup_common_args(subparser):
    # Hmm... just discovered the parents attribute. Might want to try using this
    subparser.add_argument('a', type=int, help='# sucesses in X1')
    subparser.add_argument('b', type=int, help='# sucesses in X2')
    subparser.add_argument('c', type=int, help='# failures in X1')
    subparser.add_argument('d', type=int, help='# failures in X2')

    def cs_arg(arg_string):
        return arg_string.split(',')

    prior_group = subparser.add_mutually_exclusive_group()
    prior_group.add_argument('--prior', type=float, help='Prior on Beta distributions. [default: %(default)s]', nargs=2, default=(1.0,1.0))
    prior_group.add_argument('--jeff', help='Use Jeffreys Prior on Beta distributions.', action='store_const',
            dest='prior', const=(0.5, 0.5))
    prior_group.add_argument('--rare',
            help='Use prior Beta(0.5, 1.0): corresponds to a belief that failures are more likely than successes.',
            action='store_const', dest='prior', const=(0.5, 1.0))

    subparser.add_argument('-m', '--quadr-maxiter', type=int, default=defaults['quadr_maxiter'],
            help="Value of maxiter passed to the scipy.integrate.quadrature function for CDF computation.")
    subparser.add_argument('-M', '--optim-maxiter', type=int, default=defaults['optim_maxiter'],
            help="Value of maxiter passed to scipy.optimize.brenth for MAP and PPF computation.")
    subparser.add_argument('--no-inverting', action='store_true', default=False,
            help="""Unless this flag is specified, betarat may compute the desired metrics by transforming the
            values computed from the inverse BetaRatio distribution.""" )
    subparser.add_argument('-v', '--verbose', action='store_true', default=False)


def setup_cli_br(args):
    br_params = [getattr(args, x) for x in ['a', 'b', 'c', 'd']]
    beta_rat = BetaRat(*br_params, no_inverting=args.no_inverting, prior=args.prior)
    if VERBOSE:
        print "Inverted?: ", beta_rat.inverted
    return beta_rat


def setup_ppf_args(subparsers):
    ppf_args = subparsers.add_parser('ppf',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Compute ppf(q) of beta ratio (X1/X2) distribution for table\n{}'.format(table_string))

    setup_common_args(ppf_args)
    ppf_args.add_argument('--h-init', type=float, help='initial step size in (0,1). [default: %(default)s]', default=0.005)
    ppf_args.add_argument('--simpson', action="store_true", default=False,
            help="""Use simpson_ppf method instead of solving directly for cdf = q. This method is a modified
            form of the simpson method of integration, adapted for computation of percentiles. While this
            is often faster than solving for cdf = q directly, it can also be significantly slower in cases
            of very peaked distributions, and is generally much less accurate. As such, we advise
            against this option.""")
    ppf_args.add_argument('q', type=float, help='quantile [default: %(default)s]', nargs="?", default=0.05)

    def func(args):
        br = setup_cli_br(args)
        if args.simpson:
            result = br.ppf(args.q, method="simpson", h_init=args.h_init)
        else:
            result = br.ppf(args.q, quadr_maxiter=args.quadr_maxiter, optim_maxiter=args.optim_maxiter)

        print "\nPPF({}) = {}\n".format(args.q, result)
        
    ppf_args.set_defaults(func=func)


def setup_cdf_args(subparsers):
    cdf_args = subparsers.add_parser('cdf',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Compute cdf(w) of beta ratio (X1/X2) distribution for table\n{}'.format(table_string))
    setup_common_args(cdf_args)
    cdf_args.add_argument('w', type=float, help='integrate ppf from 0 to w [default: %(default)s]', nargs="?",
            default=1.0)

    def func(args):
        br = setup_cli_br(args)
        result = br.cdf(args.w, quadr_maxiter=args.quadr_maxiter)

        print "\nCDF({}) = {}\n".format(args.w, result)

    cdf_args.set_defaults(func=func)

def setup_map_args(subparsers):
    map_args = subparsers.add_parser('map',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Compute MAP of beta ratio (X1/X2) distribution for table\n{}'.format(table_string))
    setup_common_args(map_args)

    def func(args):
        br = setup_cli_br(args)
        result = br.map(quadr_maxiter=args.quadr_maxiter, optim_maxiter=args.quadr_maxiter)

        print "\nMAP = {}\n".format(result)

    map_args.set_defaults(func=func)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""
    Beta Rat (v {})!!! Awesome Bayesian contingency table stats!

    For table
    {}
    compute either PPF, CDF or MAP of Beta ratio distribution (X1/X2).""".format(__version__, table_string))

    subparsers = parser.add_subparsers(title='subcommands', help='additional help')

    setup_ppf_args(subparsers)
    setup_cdf_args(subparsers)
    setup_map_args(subparsers)

    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose

    if VERBOSE:
        t0 = time.time()

    args.func(args)

    if VERBOSE:
        t1 = time.time()
        print "Time taken:", t1-t0, "sec\n"


if __name__ == '__main__':
    main()


