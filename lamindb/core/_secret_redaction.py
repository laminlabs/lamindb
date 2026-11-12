from __future__ import annotations

import ast
import io
import re
import tokenize
from contextlib import contextmanager

__all__ = ["redact"]

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
    r"(?P<prefix>(?P<key>[A-Za-z_][A-Za-z0-9_.\-]*)\s*(?P<op><-|:=|=|:)\s*)"
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


def _normalize_sensitive_key_name(key: str) -> str:
    normalized_key = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", key)
    normalized_key = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", normalized_key).lower()
    return normalized_key


def _is_sensitive_param_key(key: str) -> bool:
    return bool(SENSITIVE_PARAM_KEY_PATTERN.search(_normalize_sensitive_key_name(key)))


def _is_sensitive_param_value(value: object) -> bool:
    if not isinstance(value, str):
        return False
    return bool(_POSTGRES_CREDENTIALS_URL_PATTERN.match(value.strip()))


def _source_literals(source: str, *, join: bool = False):
    """Yield complete literal/comment spans without evaluating Python source."""
    offsets = [0]
    for line in source.split("\n"):
        offsets.append(offsets[-1] + len(line) + 1)
    stack: list[tokenize.TokenInfo] = []
    previous: tokenize.TokenInfo | None = None
    pending: tuple[int, int, str, str] | None = None
    try:
        for item in tokenize.generate_tokens(io.StringIO(source).readline):
            prefix, previous = previous, item
            kind = tokenize.tok_name[item.type]
            if (
                pending is not None
                and not stack
                and item.type not in (tokenize.STRING, tokenize.COMMENT, tokenize.NL)
                and kind not in ("FSTRING_START", "TSTRING_START")
                and not (
                    item.type == tokenize.NAME
                    and item.string.lower() in ("t", "tr", "rt")
                )
            ):
                yield pending
                pending = None
            if kind in ("FSTRING_START", "TSTRING_START"):
                stack.append(item)
                continue
            first = item
            if kind in ("FSTRING_END", "TSTRING_END"):
                first = stack.pop()
                if join and stack:
                    continue
            elif join and stack:
                continue
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
            elif (
                pending is not None
                and prefix is not None
                and prefix.type == tokenize.NAME
            ):
                yield pending
                pending = None
            opening = _LITERAL_OPENING_PATTERN.match(first.string)
            if kind == "STRING" and opening and "f" in opening.group(0).lower():
                try:
                    ast.parse(first.string, mode="eval")
                except SyntaxError:
                    continue
            span = (
                offsets[first.start[0] - 1] + first.start[1],
                offsets[item.end[0] - 1] + item.end[1],
                opening.group(0) if opening else "",
                opening.group("quote") if opening else "",
            )
            if join and opening:
                # Only tokens on the same logical line can concatenate. NL and
                # comments inside brackets are allowed; NEWLINE ends a statement.
                pending = (pending[0], span[1], *pending[2:]) if pending else span
            else:
                yield span
    except (tokenize.TokenError, SyntaxError):
        # Incomplete examples, notebook magics, and non-Python source are allowed.
        pass
    if pending is not None:
        yield pending


def _annotated_assignments(source: str):
    offsets = [0]
    for line in source.split("\n"):
        offsets.append(offsets[-1] + len(line) + 1)
    tokens: list[tokenize.TokenInfo] = []
    try:
        tokens.extend(tokenize.generate_tokens(io.StringIO(source).readline))
    except (tokenize.TokenError, SyntaxError):
        pass
    for index, item in enumerate(tokens[:-1]):
        if item.type != tokenize.NAME or tokens[index + 1].string != ":":
            continue
        depth = 0
        for following in tokens[index + 2 :]:
            if following.type == tokenize.NEWLINE or (
                depth == 0 and following.string in {",", ";"}
            ):
                break
            if following.type != tokenize.OP:
                continue
            if following.string in {"(", "[", "{"}:
                depth += 1
            elif following.string in {")", "]", "}"}:
                depth -= 1
                if depth < 0:
                    break
            elif following.string == "=" and depth == 0:
                start = offsets[following.end[0] - 1] + following.end[1]
                while start < len(source) and source[start].isspace():
                    start += 1
                yield offsets[item.start[0] - 1] + item.start[1], start, item.string
                break


def _redact_secret_assignments(source_code: str, _secrets=None) -> str:
    patterns = (
        _ARGV_LITERAL_PATTERN,
        _ENV_ASSIGNMENT_PATTERN,
        _KEY_VALUE_ASSIGNMENT_PATTERN,
        _QUOTED_KEY_ASSIGNMENT_PATTERN,
    )
    matches = sorted(
        [
            (match.start(), match.end(), match.group("key"))
            for pattern in patterns
            for match in pattern.finditer(source_code)
        ]
        + list(_annotated_assignments(source_code)),
        key=lambda match: match[0],
    )
    literals = list(_source_literals(source_code, join=True))
    by_start = {start: (end, opening, quote) for start, end, opening, quote in literals}
    redacted: list[str] = []
    cursor = 0
    for match_start, start, key in matches:
        if match_start < cursor:
            continue
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
                    (
                        item
                        for item in _source_literals(fragment, join=True)
                        if item[0] == 0
                    ),
                    None,
                )
                if parsed is not None:
                    literal = (start + parsed[1], parsed[2], parsed[3])
        if literal is None:
            continue
        end, opening, quote = literal
        value = source_code[start + len(opening) : end - len(quote)]
        if not _is_sensitive_param_key(key) and not _is_sensitive_param_value(value):
            continue
        if _secrets is not None:
            try:
                expression = ast.parse(
                    "(" + source_code[start:end] + ")", mode="eval"
                ).body
                if isinstance(expression, ast.JoinedStr) and all(
                    isinstance(part, ast.Constant) and isinstance(part.value, str)
                    for part in expression.values
                ):
                    # Static f-strings include decoded escapes and doubled braces.
                    # Never evaluate interpolation expressions to discover a secret.
                    secret = "".join(
                        ast.literal_eval(part) for part in expression.values
                    )
                else:
                    secret = ast.literal_eval(expression)
            except (ValueError, SyntaxError):
                secret = None
            if (
                isinstance(secret, bytes)
                and secret
                and secret != REDACTED_SECRET_VALUE.encode()
            ):
                _secrets.add(repr(secret))
                try:
                    secret = secret.decode("utf-8")
                except UnicodeDecodeError:
                    secret = None
            if isinstance(secret, str) and secret and secret != REDACTED_SECRET_VALUE:
                _secrets.add(secret)
        replacement = f"{opening}{REDACTED_SECRET_VALUE}{quote}"
        redacted.extend((source_code[cursor:start], replacement))
        cursor = end
    redacted.append(source_code[cursor:])
    return "".join(redacted)


_ARGV_LITERAL_PATTERN = re.compile(
    r"[\"']--(?P<key>[\w-]+)[\"']\s*,\s*" + _LITERAL_START
)
_FORMAT_PATTERN = re.compile(
    r"(?<![\w])(?:"
    r"(?:AKIA|ASIA)[A-Z0-9]{16}"
    r"|gh[pousr]_[A-Za-z0-9]{36}"
    r"|github_pat_[A-Za-z0-9_]{22,255}"
    r"|sk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{20,}"
    r"|xox[baprs]-[A-Za-z0-9-]{20,}"
    r"|eyJ[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+"
    r")(?![\w])"
    r"|-----BEGIN (?:[A-Z]+ )?PRIVATE KEY-----[\s\S]*?"
    r"-----END (?:[A-Z]+ )?PRIVATE KEY-----"
)
_PASSWORD_URL_PATTERN = re.compile(
    r"(?P<prefix>[a-zA-Z][a-zA-Z0-9+.-]*://[^\s/:@\"'<>]+:)"
    r"(?P<password>[^\s/@\"'<>]+)(?=@)"
)
_ENV_MAGIC_PATTERN = re.compile(
    r"(?m)^(?P<prefix>[ \t]*%env[ \t]+(?P<key>\w+)[ \t]*(?:=[ \t]*|[ \t]+))(?P<value>[^\n]+)"
)
_SHELL_TOKEN = re.compile(r"(?:[^\s\"']|\"(?:\\.|[^\"\\])*\"|'[^']*')+")


def _redact_argv(text: str) -> str:
    import shlex

    replacements = []
    sensitive = False
    for match in _SHELL_TOKEN.finditer(text):
        token = match.group()
        try:
            words = shlex.split(token)
            word = words[0] if len(words) == 1 else token
        except ValueError:
            word = token
        if sensitive:
            replacements.append((match.start(), match.end(), REDACTED_SECRET_VALUE))
            sensitive = False
        elif word.startswith("--"):
            key, sep, _ = word[2:].partition("=")
            if _is_sensitive_param_key(key):
                if sep:
                    replacements.append(
                        (match.start(), match.end(), f"--{key}={REDACTED_SECRET_VALUE}")
                    )
                else:
                    sensitive = True
    for start, end, replacement in reversed(replacements):
        text = text[:start] + replacement + text[end:]
    return text


def redact(value, kind: str = "text", *, _secrets=None):
    """Sanitize stored copies, never execution inputs. Does not evaluate code.

    Recognizes sensitive keys and known credential formats, not arbitrary secrets
    in images, encoded payloads, or unnamed random strings.
    """
    if isinstance(value, dict):
        return {
            key: REDACTED_SECRET_VALUE
            if _is_sensitive_param_key(str(key)) or _is_sensitive_param_value(item)
            else redact(item, kind)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        result = []
        sensitive = False
        for item in value:
            if sensitive:
                cleaned = REDACTED_SECRET_VALUE
            elif (
                isinstance(item, str)
                and item.startswith("--")
                and "=" in item
                and _is_sensitive_param_key(item[2:].partition("=")[0])
            ):
                cleaned = item.partition("=")[0] + "=" + REDACTED_SECRET_VALUE
            else:
                cleaned = redact(item, kind)
            result.append(cleaned)
            sensitive = (
                isinstance(item, str)
                and item.startswith("--")
                and "=" not in item
                and _is_sensitive_param_key(item[2:])
            )
        return tuple(result) if isinstance(value, tuple) else result
    if not isinstance(value, str):
        return value

    def env_magic(match):
        if not _is_sensitive_param_key(match["key"]):
            return match[0]
        # IPython strips surrounding whitespace but keeps quotes and escapes.
        secret = match["value"].strip()
        if _secrets is not None and secret and secret != REDACTED_SECRET_VALUE:
            _secrets.add(secret)
        return match["prefix"] + REDACTED_SECRET_VALUE

    # Collect the original value before assignment redaction can replace it.
    value = _ENV_MAGIC_PATTERN.sub(env_magic, value)
    value = _redact_secret_assignments(value, _secrets)

    def formats(text):
        text = _FORMAT_PATTERN.sub(REDACTED_SECRET_VALUE, text)
        return _PASSWORD_URL_PATTERN.sub(
            lambda m: m["prefix"] + REDACTED_SECRET_VALUE, text
        )

    if kind == "source":
        # Restrict format matching to literals/comments so token-like identifiers
        # remain valid Python. Replace whole literals, including f-strings.
        parts: list[str] = []
        cursor = 0
        for start, end, opening, quote in sorted(
            _source_literals(value), key=lambda item: (item[0], -item[1])
        ):
            if start < cursor:
                continue
            literal = value[start:end]
            cleaned = formats(literal)
            body = literal[len(opening) : -len(quote)] if quote else ""
            if body.startswith("--") and "=" in body:
                key = body[2:].partition("=")[0]
                if _is_sensitive_param_key(key):
                    cleaned = f"{opening}--{key}={REDACTED_SECRET_VALUE}{quote}"
            if cleaned != literal and quote and not body.startswith("--"):
                cleaned = f"{opening}{REDACTED_SECRET_VALUE}{quote}"
            parts.extend((value[cursor:start], cleaned))
            cursor = end
        return "".join(parts) + value[cursor:]
    return formats(_redact_argv(value))


def _redact_record(record) -> None:
    # Access only loaded fields: partial saves must not fetch deferred fields.
    fields = {
        "lamindb.run": {"cli_args": "argv", "params": "text"},
        "lamindb.transform": {"source_code": "source"},
    }.get(record._meta.label_lower, {})
    if (
        record._meta.label_lower == "lamindb.artifact"
        and record.__dict__.get("kind") == "__lamindb_run__"
    ):
        path = getattr(record, "_local_filepath", None)
        if path is not None:
            _check_run_artifact(path, record.size, record.hash)
    for field, kind in fields.items():
        if field in record.__dict__:
            setattr(record, field, redact(record.__dict__[field], kind))


class _UnsafeRunArtifact(ValueError):
    pass


def _redact_html(text: str, *, _secrets: set[str] | None = None) -> str:
    # Reports are already rendered. Never patch HTML with source-code regexes:
    # a value may span tags, entities, attributes, or an embedded script.
    from html import unescape
    from html.parser import HTMLParser

    class Text(HTMLParser):
        def __init__(self):
            super().__init__(convert_charrefs=True)
            self.parts = []
            self.attributes = []

        def handle_data(self, data):
            self.parts.append(data)

        def handle_starttag(self, tag, attrs):
            self.attributes.extend(value for _, value in attrs if value is not None)

        def handle_comment(self, data):
            self.parts.append(data)

    parsed = Text()
    parsed.feed(text)
    parsed.close()
    # Check both separated and joined nodes (syntax highlighting splits tokens).
    for candidate in (
        unescape(text),
        "".join(parsed.parts),
        "\n".join(parsed.parts),
        *parsed.attributes,
    ):
        if redact(candidate) != candidate or any(
            secret in candidate or repr(secret)[1:-1] in candidate
            for secret in (_secrets or ())
        ):
            raise _UnsafeRunArtifact(
                "skipping run attachment: rendered HTML contains a secret that cannot be safely removed"
            )
    return text


def _read_run_text(path):
    try:
        return path.read_bytes().decode("utf-8")
    except UnicodeError as exc:
        raise _UnsafeRunArtifact(
            "skipping run attachment: content is not supported UTF-8 text"
        ) from exc


def _stage_run_artifact(path):
    """Return an owned, sanitized snapshot and its temporary-directory lifetime."""
    from pathlib import Path
    from tempfile import TemporaryDirectory

    if not isinstance(path, (str, Path)) or "://" in str(path):
        raise _UnsafeRunArtifact(
            "run attachments must be local files to check before upload"
        )
    path = Path(path)
    stage = TemporaryDirectory(prefix="lamin-redacted-")
    target = Path(stage.name) / path.name
    try:
        if path.is_dir():
            target.mkdir()
        paths = list(path.rglob("*")) if path.is_dir() else [path]
        for source in paths:
            if source.is_symlink():
                raise _UnsafeRunArtifact(
                    "skipping run attachment: symbolic links cannot be checked safely"
                )
            if source.is_dir():
                continue
            destination = target / source.relative_to(path) if path.is_dir() else target
            destination.parent.mkdir(parents=True, exist_ok=True)
            text = _read_run_text(source)
            cleaned = (
                _redact_html(text)
                if source.suffix.lower() in {".html", ".htm"}
                else redact(text)
            )
            destination.write_bytes(cleaned.encode("utf-8"))
    except BaseException:
        stage.cleanup()
        raise
    return target, stage


def _check_run_artifact(path, size, hash) -> None:
    """Fail before upload if a snapshot was changed after hashing."""
    from lamindb_setup.core.hashing import hash_dir, hash_file

    clean, stage = _stage_run_artifact(path)
    try:
        hasher = hash_dir if path.is_dir() else hash_file
        if hasher(path)[:2] != (size, hash) or hasher(clean)[:2] != (size, hash):
            raise _UnsafeRunArtifact(
                "run attachment changed after preparation; recreate it before saving"
            )
    finally:
        stage.cleanup()


@contextmanager
def _skip_unsafe_run_artifact():
    from lamindb_setup import logger

    try:
        yield
    except _UnsafeRunArtifact as exc:
        logger.warning(str(exc))


def _redact_notebook(notebook) -> None:
    # Remember literal secrets within this notebook only, to also remove their
    # bare values from printed outputs. Never retain them globally or in metadata.
    secrets: set[str] = set()
    for cell in notebook.cells:
        cell.source = redact(cell.source, "source", _secrets=secrets)

    def clean(value, *, html=False):
        if isinstance(value, dict):
            return {k: clean(v, html=k == "text/html") for k, v in value.items()}
        if isinstance(value, list):
            if html and all(isinstance(part, str) for part in value):
                value = "".join(value)
            else:
                return [clean(v) for v in value]
        if isinstance(value, str):
            from html import escape

            for secret in sorted(secrets, key=len, reverse=True):
                variants = {
                    variant
                    for text in (secret, repr(secret)[1:-1])
                    for variant in (text, escape(text), escape(text, quote=False))
                }
                for variant in sorted(variants, key=len, reverse=True):
                    value = value.replace(variant, REDACTED_SECRET_VALUE)
            if html:
                try:
                    return _redact_html(value, _secrets=secrets)
                except _UnsafeRunArtifact:
                    # Drop only this HTML representation, preserving safe outputs.
                    return f"<pre>{REDACTED_SECRET_VALUE}</pre>"
            return redact(value)
        return value

    for cell in notebook.cells:
        if cell.cell_type == "code":
            cell.outputs = clean(cell.outputs)
