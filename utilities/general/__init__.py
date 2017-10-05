import inspect
from tabulate import tabulate

def getArgDesc(func, args):
    a = inspect.getargspec(func)
    if a.defaults is not None:
        d = dict(zip(a.args[-len(a.defaults):],a.defaults))
    else:
        d = {}
    d.update(args)
    argStr = "("
    first = True
    for arg in a.args:
        if arg in d:
            if not first:
                argStr += ", "
            else:
                first = False
            argStr += arg + " " + str(d[arg])
    argStr += ")"
    return argStr
    #return str(d).replace(':', '').replace('\'', '').replace('{', '(').replace('}', ')')

NOTABULATEFORMAT = "noformat"  #tabulate automatically formats if it look like a string
def getFormattedHTMLTable(contents, headers, floatFmtFigs=None, numalign=None, stralign=None, border=False):
    if floatFmtFigs is None:
        floatfmt = "g"   #This is the default value in the tabulate parameter list.
    else:
        floatfmt = '.'+str(floatFmtFigs)+'f'
    outStr = tabulate(contents, headers=headers, numalign=numalign, stralign=stralign, floatfmt=floatfmt, tablefmt="html").replace(NOTABULATEFORMAT,"")
    if border:
        replaceStr = "<table style=\"border-collapse:collapse;\" border=\"1\" cellpadding=\"5\" cellspacing=\"0\">"
        outStr = outStr.replace("<table>", replaceStr)
        replaceStr = "<th style=\"background: lightgrey;\" "
        return outStr.replace("<th", replaceStr)
    else:
        replaceStr = "<table cellspacing=\"10\" style=\"float: left\">"
        return outStr.replace("<table>", replaceStr)