import inspect
from tabulate import tabulate

def getArgDesc(func, args):
    a = inspect.getargspec(func)
    if a.defaults is not None:
        d = dict(zip(a.args[-len(a.defaults):],a.defaults))
    else:
        d = {}
    d.update(args)
    return str(d).replace(':', '').replace('\'', '').replace('{', '(').replace('}', ')')

def getFormattedHTMLTable(contents, floatSigFigs, headers):
    outStr = tabulate(contents, headers=headers, numalign="decimal", floatfmt='.'+str(floatSigFigs)+'f', tablefmt="html")
    return outStr.replace("<table>", "<table cellspacing=\"10\" style=\"float: left\">")