from __future__ import annotations

import re

REDACTED_SECRET_VALUE = "***REDACTED***"  # noqa: S105
SENSITIVE_PARAM_KEY_PATTERN = re.compile(
    r"(^|[_\-.])(api[_-]?key|access[_-]?key|secret|token|password|passwd|private[_-]?key|client[_-]?secret)($|[_\-.])"
)

# Match only quoted literals in assignments, e.g.:
# - my_secret = "value"
# - my.secret: "value"
# - mySecret := "value"
# We intentionally do not match unquoted RHS values to avoid false positives like
# type annotations (`api_key: str`) or variable forwarding (`api_key=api_key`).
_KEY_VALUE_ASSIGNMENT_PATTERN = re.compile(
    r"(?P<prefix>(?P<key>[A-Za-z_][A-Za-z0-9_.\-]*)\s*(?P<op>:=|=|:)\s*)"
    r"(?P<value>(?P<quote>['\"`])(?P<quoted>.*?)(?P=quote))"
)

# Match: os.environ["API_KEY"] = "value"
_ENV_ASSIGNMENT_PATTERN = re.compile(
    r"(?P<prefix>os\.environ\[\s*(?P<kquote>['\"])(?P<key>[^'\"]+)(?P=kquote)\s*\]\s*=\s*)"
    r"(?P<value>(?P<quote>['\"`])(?P<quoted>.*?)(?P=quote))"
)

# Match: {"client_secret": "value"}
_QUOTED_KEY_ASSIGNMENT_PATTERN = re.compile(
    r"(?P<prefix>(?P<kquote>['\"])(?P<key>[^'\"]+)(?P=kquote)\s*:\s*)"
    r"(?P<value>(?P<quote>['\"`])(?P<quoted>.*?)(?P=quote))"
)

# We intentionally treat env lookups as safe/re-runnable references, not embedded secrets.
# Examples that should remain unchanged:
# - api_key = os.getenv("OPENAI_API_KEY")
# - api_key = getenv("OPENAI_API_KEY")
# - api_key = os.environ["OPENAI_API_KEY"]
# - api_key = os.environ.get("OPENAI_API_KEY")
_ENV_REFERENCE_VALUE_PATTERN = re.compile(
    r"^(os\.getenv\(.+\)|getenv\(.+\)|os\.environ\[[^\]]+\]|os\.environ\.get\(.+\))$"
)


def normalize_sensitive_key_name(key: str) -> str:
    normalized_key = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", key)
    normalized_key = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", normalized_key).lower()
    return normalized_key


def is_sensitive_param_key(key: str) -> bool:
    return bool(SENSITIVE_PARAM_KEY_PATTERN.search(normalize_sensitive_key_name(key)))


def _redact_assignment_match(match: re.Match[str]) -> str:
    key = match.group("key")
    if not is_sensitive_param_key(key):
        return match.group(0)
    # Redact only hardcoded values, not environment-based references.
    # This preserves reproducibility for source code that reads secrets from env vars.
    raw_value = match.group("value")
    if _ENV_REFERENCE_VALUE_PATTERN.match(raw_value):
        return match.group(0)
    quote = match.group("quote")
    redacted_value = (
        f"{quote}{REDACTED_SECRET_VALUE}{quote}"
        if quote is not None
        else REDACTED_SECRET_VALUE
    )
    return f"{match.group('prefix')}{redacted_value}"


def redact_secrets_in_source_code(source_code: str) -> tuple[str, int]:
    redaction_count = 0

    def replace_with_count(match: re.Match[str]) -> str:
        nonlocal redaction_count
        replaced = _redact_assignment_match(match)
        if replaced != match.group(0):
            redaction_count += 1
        return replaced

    redacted = _ENV_ASSIGNMENT_PATTERN.sub(replace_with_count, source_code)
    redacted = _KEY_VALUE_ASSIGNMENT_PATTERN.sub(replace_with_count, redacted)
    redacted = _QUOTED_KEY_ASSIGNMENT_PATTERN.sub(replace_with_count, redacted)
    return redacted, redaction_count
