class EmptyFileError(OSError):
    """
    This exception is raised when reading an input file with no content.
    """

    def __init__(self, path: str = "", message: str = "File is empty.") -> None:
        self.path = path
        self.message = message

    def __str__(self):
        if self.path:
            self.message = f"File '{self.path}' is empty."

        return self.message
