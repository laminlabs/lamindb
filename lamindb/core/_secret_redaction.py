from __future__ import annotations

import ast
import io
import re
import tokenize

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
_LITERAL_START = r"(?=(?i:br|rb|fr|rf|tr|rt|[rubft])?['\"`])"
_KEY_VALUE_ASSIGNMENT_PATTERN = re.compile(
    r"(?P<prefix>(?P<key>[A-Za-z_][A-Za-z0-9_.\-]*)\s*(?P<op>:=|=|:)\s*)"
    + _LITERAL_START
)

# Match: os.environ["API_KEY"] = "value"
_ENV_ASSIGNMENT_PATTERN = re.compile(
    r"(?P<prefix>os\.environ\[\s*(?P<kquote>['\"])(?P<key>[^'\"]+)(?P=kquote)\s*\]\s*=\s*)"
    + _LITERAL_START
)

# Match: {"client_secret": "value"}
_QUOTED_KEY_ASSIGNMENT_PATTERN = re.compile(
    r"(?P<prefix>(?P<kquote>['\"])(?P<key>[^'\"]+)(?P=kquote)\s*:\s*)" + _LITERAL_START
)

_LITERAL_OPENING_PATTERN = re.compile(
    r"(?i:br|rb|fr|rf|tr|rt|[rubft])?(?P<quote>\"\"\"|'''|['\"`])"
)
# Backticks are supported by the existing source redactor for non-Python code.
_BACKTICK_LITERAL_PATTERN = re.compile(r"`(?:\\[\s\S]|[^`\\])*`")

# Match PostgreSQL URLs that include inline credentials:
# - postgresql://user:password@host:5432/dbname
# - postgres://user:password@host/dbname?sslmode=require
_POSTGRES_CREDENTIALS_URL_PATTERN = re.compile(
    r"^postgres(?:ql)?://[^:@/\s]+:[^@/\s]+@[^/\s]+(?:/[^\s]*)?$",
    re.IGNORECASE,
)


def normalize_sensitive_key_name(key: str) -> str:
    normalized_key = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", key)
    normalized_key = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", normalized_key).lower()
    return normalized_key


def is_sensitive_param_key(key: str) -> bool:
    return bool(SENSITIVE_PARAM_KEY_PATTERN.search(normalize_sensitive_key_name(key)))


def is_sensitive_param_value(value: object) -> bool:
    if not isinstance(value, str):
        return False
    return bool(_POSTGRES_CREDENTIALS_URL_PATTERN.match(value.strip()))


def _source_literals(source: str):
    """Yield complete literal/comment spans without evaluating Python source."""
    offsets = [0]
    for line in source.split("\n"):
        offsets.append(offsets[-1] + len(line) + 1)
    stack = []
    previous = None
    try:
        for item in tokenize.generate_tokens(io.StringIO(source).readline):
            prefix, previous = previous, item
            kind = tokenize.tok_name[item.type]
            if kind in ("FSTRING_START", "TSTRING_START"):
                stack.append(item)
                continue
            first = item
            if kind in ("FSTRING_END", "TSTRING_END"):
                first = stack.pop()
            elif item.type not in (tokenize.STRING, tokenize.COMMENT):
                continue
            if (
                item.type == tokenize.STRING
                and prefix is not None
                and prefix.type == tokenize.NAME
                and prefix.string.lower() in ("t", "tr", "rt")
                and prefix.end == item.start
            ):
                # Older tokenizers split template prefixes from their strings.
                # Validate as an f-string so nested quotes cannot cause a partial match.
                try:
                    ast.parse(
                        prefix.string.lower().replace("t", "f") + item.string,
                        mode="eval",
                    )
                except SyntaxError:
                    continue
                first = item._replace(
                    start=prefix.start, string=prefix.string + item.string
                )
            opening = _LITERAL_OPENING_PATTERN.match(first.string)
            if kind == "STRING" and opening and "f" in opening.group(0).lower():
                try:
                    ast.parse(first.string, mode="eval")
                except SyntaxError:
                    continue
            yield (
                offsets[first.start[0] - 1] + first.start[1],
                offsets[item.end[0] - 1] + item.end[1],
                opening.group(0) if opening else "",
                opening.group("quote") if opening else "",
            )
    except (tokenize.TokenError, SyntaxError):
        # Incomplete examples, notebook magics, and non-Python source are allowed.
        pass


def redact_secrets_in_source_code(source_code: str) -> tuple[str, int]:
    patterns = (
        _ENV_ASSIGNMENT_PATTERN,
        _KEY_VALUE_ASSIGNMENT_PATTERN,
        _QUOTED_KEY_ASSIGNMENT_PATTERN,
    )
    matches = sorted(
        (match for pattern in patterns for match in pattern.finditer(source_code)),
        key=lambda match: match.start(),
    )
    literals = list(_source_literals(source_code))
    by_start = {start: (end, opening, quote) for start, end, opening, quote in literals}
    redacted, cursor, count = [], 0, 0
    for match in matches:
        if match.start() < cursor:
            continue
        start = match.end()
        literal = by_start.get(start)
        if literal is None:
            # Preserve redaction inside examples, without consuming their closing quote.
            limit = min(
                (
                    end - len(quote)
                    for begin, end, _, quote in literals
                    if begin < start < end
                ),
                default=len(source_code),
            )
            fragment = source_code[start:limit]
            if fragment.startswith("`"):
                backtick = _BACKTICK_LITERAL_PATTERN.match(fragment)
                literal = (start + backtick.end(), "`", "`") if backtick else None
            else:
                parsed = next(
                    (item for item in _source_literals(fragment) if item[0] == 0), None
                )
                if parsed is not None:
                    literal = (start + parsed[1], parsed[2], parsed[3])
        if literal is None:
            continue
        end, opening, quote = literal
        value = source_code[start + len(opening) : end - len(quote)]
        if not is_sensitive_param_key(
            match.group("key")
        ) and not is_sensitive_param_value(value):
            continue
        replacement = f"{opening}{REDACTED_SECRET_VALUE}{quote}"
        redacted.extend((source_code[cursor:start], replacement))
        count += replacement != source_code[start:end]
        cursor = end
    redacted.append(source_code[cursor:])
    return "".join(redacted), count
