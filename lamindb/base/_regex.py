"""Python regex semantics shared by feature definitions and curator backends."""

import re

from lamindb.errors import ValidationError


def compile_regex(pattern: str) -> re.Pattern:
    if not isinstance(pattern, str):
        raise TypeError("Feature.regex must be a string or None")
    try:
        return re.compile(pattern)
    except re.error as error:
        raise ValidationError(f"Invalid feature regex {pattern!r}: {error}") from error


def matches_regex(value: object, pattern: re.Pattern) -> bool:
    """Match a complete non-null string without coercing its value."""
    return isinstance(value, str) and pattern.fullmatch(value) is not None


def regex_error(feature) -> str:
    return f"Feature '{feature.name}' must fully match regex {feature.regex!r}"
