def escape_regex(regex):
    """
    Escapes all common metacharacters of a given regular expression with a '\'.
    Parameters
    ----------
    regex: str
        Regular expression that should be excaped.
    Returns
    -------
        Escaped regular expression.
    """
    modyfied_regex = ""
    for a in regex:
        if a in r"^[]{}().$*\+|?<>=":
            modyfied_regex += fr"\{a}"
        else:
            modyfied_regex += a
    return modyfied_regex