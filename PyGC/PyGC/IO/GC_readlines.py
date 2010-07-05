def GCreadlines(fobject, comment='#'):
    """Strip files from comments.

    fobject: file object to be read
    comment: string indicating comment.
    """

    lines = []
    for l in fobject.readlines():
        # ignore comments and empty lines
        l = l.strip()
        pos = l.find(comment)
        if pos>-1:
            l = l[:pos]
        if len(l)==0:
            continue
        lines.append(l)
    return lines
