# ------------------------------------------
#
# Matthieu Defrance - ULB
# commandify
#
# ------------------------------------------

def commandify(f, allowPositionnalArguments = True):
    # import inspect
    import argparse

    # get function signature
    fargs     = f.__code__.co_varnames[:f.__code__.co_argcount]
    fname     = f.__code__.co_name
    fname     = f.__name__
    fhelp     = f.__doc__

    fdefaults = f.__defaults__ or {}
    fdefaults = dict(zip(fargs[-len(fdefaults):], fdefaults))

    # add arguments
    parser = argparse.ArgumentParser(description = '%s' % (fhelp))
    for farg in fargs:
        if fdefaults.has_key(farg) or allowPositionnalArguments == False:
            parser.add_argument('--%s' % farg, help = '%s' % farg)
        else:
            parser.add_argument('%s' % farg, help = '%s' % farg)

    # parse args and call function
    args   = parser.parse_args()
    kwargs = dict(args._get_kwargs())
    f(**kwargs)
