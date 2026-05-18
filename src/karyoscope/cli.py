"""Top-level command-line interface for KaryoScope.

This module wires the eight v0.1 subcommands into a single ``karyoscope``
entry point. Each subcommand lives in its own module under
``karyoscope.commands`` so that they can grow independently.
"""

from __future__ import annotations

import click

from karyoscope._version import __version__
from karyoscope.commands import (
    annotate,
    bin_cmd,
    centromeres,
    download,
    info,
    karyotype,
    scaffold,
    version,
)

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
    "max_content_width": 100,
}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(
    __version__,
    "-V",
    "--version",
    message="karyoscope %(version)s",
)
def main() -> None:
    """KaryoScope: rapid, alignment-free sequence annotation for the pangenome era.

    Run a subcommand with ``--help`` to see its options, e.g.:

    \b
        karyoscope annotate --help

    For more information, see https://github.com/barthel-lab/KaryoScope.
    """


# Register subcommands. The order here determines the order in `--help`.
main.add_command(download.cmd, name="download")
main.add_command(annotate.cmd, name="annotate")
main.add_command(scaffold.cmd, name="scaffold")
main.add_command(bin_cmd.cmd, name="bin")
main.add_command(centromeres.cmd, name="centromeres")
main.add_command(karyotype.cmd, name="karyotype")
main.add_command(info.cmd, name="info")
main.add_command(version.cmd, name="version")


if __name__ == "__main__":
    main()
