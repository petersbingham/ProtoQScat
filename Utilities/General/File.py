from sys import platform as _platform

def sep():
    if _platform == "win32":
        return "\\"
    else:
        return "/"

def fixPath(path):
    if _platform == "win32":
        return u"\\\\?\\"+unicode(path)
    else:
        return path