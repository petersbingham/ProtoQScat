import inspect

def getArgDesc(func, args):
    a = inspect.getargspec(func)
    if a.defaults is not None:
        d = dict(zip(a.args[-len(a.defaults):],a.defaults))
    else:
        d = {}
    d.update(args)
    return str(d).replace(':', '').replace('\'', '')