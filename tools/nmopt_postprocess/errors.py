"""Errors shared by post-processing stages."""


class PostprocessError(RuntimeError):
    """An input artifact cannot be processed by a post-processor."""
