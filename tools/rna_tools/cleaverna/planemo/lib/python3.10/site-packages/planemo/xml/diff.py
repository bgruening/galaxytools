def diff(x1, x2, reporter=None):
    """Return 0 if and only if the XML has the same content."""
    compare = xml_compare(x1, x2, reporter)
    return_val = 0 if compare else 1
    return return_val


# From
# bitbucket.org/ianb/formencode/src/tip/formencode/doctest_xml_compare.py
# with (PSF license)
def xml_compare(x1, x2, reporter=None):
    if reporter is None:

        def reporter(x):
            return None

    if x1.tag != x2.tag:
        reporter(f"Tags do not match: {x1.tag} and {x2.tag}\n")
        return False
    for name, value in x1.attrib.items():
        if x2.attrib.get(name) != value:
            reporter("Attributes do not match: %s=%r, %s=%r\n" % (name, value, name, x2.attrib.get(name)))
            return False
    for name in x2.attrib.keys():
        if name not in x1.attrib:
            reporter("x2 has an attribute x1 is missing: %s\n" % name)
            return False
    if not text_compare(x1.text, x2.text):
        reporter(f"text: {x1.text!r} != {x2.text!r}\n")
        return False
    if not text_compare(x1.tail, x2.tail):
        reporter(f"tail: {x1.tail!r} != {x2.tail!r}\n")
        return False
    return _compare_children(x1, x2, reporter)


def _compare_children(x1, x2, reporter):
    cl1 = list(x1)
    cl2 = list(x2)
    if len(cl1) != len(cl2):
        reporter("children length differs, %i != %i\n" % (len(cl1), len(cl2)))
        return False
    i = 0
    for c1, c2 in zip(cl1, cl2):
        i += 1
        if not xml_compare(c1, c2, reporter=reporter):
            reporter("children %i do not match: %s\n" % (i, c1.tag))
            return False
    return True


def text_compare(t1, t2):
    if not t1 and not t2:
        return True
    return (t1 or "").strip() == (t2 or "").strip()
