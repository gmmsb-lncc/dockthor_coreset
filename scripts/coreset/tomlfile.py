"""The <pdb>.toml files of the qualification set: Python values as TOML literals, the file layout."""

import json
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
