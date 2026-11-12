"""The <pdb>.toml files of the qualification set: Python values as TOML literals, the --set option, the file layout."""

import argparse
import json
import tomllib
from pathlib import Path


def toml_value(value) -> str:
    """Format a Python value as a TOML literal."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, list | tuple):
        return "[" + ",".join(toml_value(v) for v in value) + "]"
    # JSON string escaping is valid for TOML basic strings.
    return json.dumps(str(value))


def parse_set_option(text: str) -> tuple[str, object]:
    """Parse KEY=VALUE; VALUE is read as a TOML literal, else kept as a string."""
    key, sep, raw = text.partition("=")
    if not sep or not key.strip():
        raise argparse.ArgumentTypeError(f"expected KEY=VALUE, got '{text}'")
    try:
        value = tomllib.loads(f"v = {raw}")["v"]
    except tomllib.TOMLDecodeError:
        value = raw
    return key.strip().lstrip("-"), value


def toml_text(pdb: str, config: dict, comments: dict[str, str] | None = None) -> str:
    """The <pdb>.toml of the set: a title line, a blank line, one key=value per line in the order of config.

    comments: {key: text}, written as a `# text` line just above that key."""
    lines = [f"## DockThor qualification set: {pdb}", ""]
    for key, value in config.items():
        if comments and key in comments:
            lines.append(f"# {comments[key]}")
        lines.append(f"{key}={toml_value(value)}")
    return "\n".join(lines) + "\n"


def write_toml(path: Path, pdb: str, config: dict, comments: dict[str, str] | None = None) -> None:
    path.write_text(toml_text(pdb, config, comments))
