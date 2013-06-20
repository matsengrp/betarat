

def plot_fn(f, a, b, points=100, label=None):
    """ Nice utility function for plotting functions using pylab."""
    # Don't want to force users to have pylab, since many may not actually need this function
    try:
        import pylab
    except ImportError:
        print "Could not import pylab. Make sure that the library is installed and try again."
        return
    inc = float(b-a) / points
    xs = pylab.arange(a, b, inc)
    ys = [f(x) for x in xs]
    if not label:
        label = str(f)
    pylab.plot(xs, ys, label=label)
    pylab.legend()


