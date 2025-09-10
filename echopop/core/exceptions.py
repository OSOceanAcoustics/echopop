class EchopopValidationError(Exception):
    """
    Custom exception class for Echopop validation errors.

    This exception is raised when validation of input data, parameters, or
    configurations fails within the Echopop package. It wraps the original
    exception to provide consistent error handling.

    Parameters
    ----------
    exception : Exception, optional
        The original exception that caused the validation error, by default None

    Attributes
    ----------
    exception : Exception or None
        The wrapped original exception

    Examples
    --------
    >>> try:
    ...     # Some validation code
    ...     raise ValueError("Invalid parameter")
    ... except ValueError as e:
    ...     raise EchopopValidationError(e)
    """

    def __init__(self, exception: Exception = None):
        self.exception = exception
        super().__init__(str(exception) if exception else "")
