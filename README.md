# BetaRat

The _BetaRat_ distribution is the ratio of two _Beta_ distributed random variables.
Since Bayesian analysis of a binomially distributed random variable naturally leads to the _Beta_ distribution, the _BetaRat_ distribution can be used as a way of comparing two binomially distributed random variables.
This can be particularly useful for considering analysis of contingency tables from a Bayesian perspective.
Indeed, the "gold standard" of contingency table tests for small to moderate trial counts is the Fisher exact test.
This test has been shown to be inappropriate in cases where the margins of the contingency table are not fixed (that is, in almost any situation other than [the lady tasting tea](http://en.wikipedia.org/wiki/Lady_tasting_tea), for which it was designed).
In such cases, it has been shown that the test is overly conservative, producing more false negatives than one would expect from a given significance level.
For a full treatment of this see [Sekhon](http://polmeth.wustl.edu/media/Paper/SekhonFisherTest.pdf).

The _BetaRat_ distribution has appropriate error rates, behaves more sensibly in cases where there are few or no successes in one of the categories, and is more readily interpretable as the relative ratio of the probabilities of the two independent binomial distributions.
This distribution can be used as a drop in replacement for the Fisher exact test by considering the _CDF(1.0)_ as one would a P-value.
This can be used in deciding whether there is evidence of a statistically higher underlying probability of success in one category versus another.
We recommend using the [Maximum a Posteriori](http://en.wikipedia.org/wiki/Maximum_a_posteriori_estimation) as the 'best' estimate of the relative probability ratio.


## CLI

Once installed, this package can be used via a command line interface with the `betarat` command.
The `betarat` command takes two subcommands, `cdf`, `pdf` and `ppf`.
Help for each of these commands can be obtained by entering `betarat [cmd] -h` at the command line.

Some example usage...

    betarat cdf 5 7 20 19
    => CDF(1.0) = 0.71149240327133

    betarat cdf 5 7 20 19 2.0
    => CDF(2.0) = 0.977779952999632

    betarat map 5 7 20 19
    => MAP = 0.62937009157218


## Library use

This package can be used within other python scripts as a library by importing from the `betarat` package.

    from betarat import BetaRat

    br = BetaRat(4, 5, 1, 9)

    print "PDF(0.8) is", br.pdf(0.8)
    print "CDF(1.0) is", br.cdf(1.0)
    print "Maximum a Posteriori (MAP) is", br.map()
    print "0.05 quantile, PPF(0.05), is", br.ppf(0.05)


## A note about quantiles

For computation of quantiles (`BetaRat.pdf`), two methods are available.
The suggested and default method is use of `scipy.optimize` to directly solve for the CDF(x) = q.
An alternative method employs a modified implementation of Simpson's method of numerical integration which "stops' after a certain.
This implementation can be chosen by passing `method='simpson'` to `BetaRat.pdf`.


## Installation

Currently, BetaRat has only been tested on Linux systems, but OSX installation should be fairly easy.
Below, we'll guide you step by step through installation on Ubuntu Linux, and provide some general guidance for dependencies on other 'nix systems.
If you do set this up on another system and would like to share any useful tips, please feel free to add a page to [the wiki](https://github.com/fhcrc/betarat/wiki).
If you have trouble, please submit a [github issue](https://github.com/fhcrc/betarat/issues), and we'll try to help you out.

### Pre-requisites

First, you'll need a version of python (this program has only been tested on python 2.7) installed.
Many Linux distros come with python 2.7 preinstalled, as does OSX.

You will also need the following python libraries:

* numpy
* scipy
* cython
* mpmath
* setuptools (and/or pip)

You will also need the `libgmp` library to be accessible on your system.

If you are running Ubuntu, you should be able to handle all of these dependencies with the following

    sudo apt-get install libgmp-dev python-pip python-scipy python-numpy python-mpmath cython

If you are running OSX or another flavor of Linux, we recommend first installing [pip](https://pypi.python.org/pypi/pip), and using that to install the remaining python packages.
These other packages may require additional libraries (such as libblas and liblapack) which you will have to install manually.
The python libraries' websites should provide all necessary information on these lines.

One last note: mpmath is currently incompatible with gmpy2, but should work fine with gmpy(1).
If you need gmpy2 for some other use on your system, you may have to set up a [virtualenv](https://pypi.python.org/pypi/virtualenv) for one project or the other.

### With that out of the way...

Now we can install BetaRat!
From a terminal, enter these lines:

    # I like to download source code to a src directory in my home folder
    mkdir -p ~/src
    cd ~/src

    # download and unzip the betarat source code
    wget https://github.com/fhcrc/betarat/archive/master.zip -O betarat.zip
    unzip betarat.zip
    cd betarat-master

    # install (sudo may not be necessary, depending on how you set up your python environment)
    sudo python setup.py install

And there you have it!
You can try running `betarat -h` from the command line to test your installation.

