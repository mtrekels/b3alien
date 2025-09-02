from IPython import get_ipython


def detect_runtime():
    """
    Detects the runtime environment where the code is executed.

    Returns:
        str: One of "Jupyter Notebook", "IPython Terminal", or "Standard Python Script"
    """
    try:
        from IPython import get_ipython
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return "Jupyter Notebook"
        elif shell == 'TerminalInteractiveShell':
            return "IPython Terminal"
        else:
            return "Other IPython"
    except (ImportError, AttributeError, NameError):
        return "Standard Python Script"

def in_jupyter():
    """
        Returns True if running inside a Jupyter Notebook or Lab.
    """
    return detect_runtime() == "Jupyter Notebook"

def in_ipython():
    """
        Returns True if running inside any IPython shell (not standard Python).
    """
    return detect_runtime() != "Standard Python Script"

def in_script():
    """
        Returns True if running in a standard Python script (non-interactive).
    """
    return detect_runtime() == "Standard Python Script"