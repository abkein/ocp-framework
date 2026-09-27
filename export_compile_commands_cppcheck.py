import json
import os
import shlex
from pathlib import Path

source = Path("build/compile_commands.json")
destination = Path("build/compile_commands.cppcheck.json")
database = json.loads(source.read_text())
nix_flags = shlex.split(os.environ["NIX_CFLAGS_COMPILE"])

for entry in database:
    arguments = (entry["arguments"] if "arguments" in entry
                else shlex.split(entry["command"]))
    entry["arguments"] = [
        "-I" if arg == "-isystem" else arg
        for arg in arguments + nix_flags
    ]
    entry.pop("command", None)

destination.write_text(json.dumps(database, indent=2) + "\n")
