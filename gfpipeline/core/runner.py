"""External tool invocation wrapper.

Provides ToolRunner which encapsulates subprocess calls with support for
dry_run mode, verbose logging, and structured error reporting.
"""

import subprocess
from typing import Optional

from gfpipeline.core.exceptions import ExternalToolError
from gfpipeline.core.logger import get_logger

log = get_logger(__name__)


class ToolRunner:
    """Encapsulates subprocess calls for external bioinformatics tools.

    Args:
        dry_run: If True, print commands without executing them.
        verbose: If True, log full command lines at DEBUG level.
    """

    def __init__(self, dry_run: bool = False, verbose: bool = False) -> None:
        self.dry_run = dry_run
        self.verbose = verbose

    def run(
        self,
        cmd: list[str],
        cwd: Optional[str] = None,
    ) -> subprocess.CompletedProcess:
        """Execute an external command.

        Args:
            cmd: Command as a list of strings.
            cwd: Working directory for the subprocess (optional).

        Returns:
            CompletedProcess instance on success.

        Raises:
            ExternalToolError: If the command exits with a non-zero return code.
        """
        cmd_str = " ".join(str(c) for c in cmd)

        if self.verbose or self.dry_run:
            log.debug("CMD: %s", cmd_str)

        if self.dry_run:
            log.info("[dry-run] %s", cmd_str)
            return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

        log.debug("Running: %s", cmd_str)
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            tool = cmd[0] if cmd else "unknown"
            raise ExternalToolError(
                tool=tool,
                cmd=cmd,
                returncode=result.returncode,
                stderr=result.stderr,
            )

        if self.verbose and result.stdout:
            log.debug("STDOUT: %s", result.stdout.strip())

        return result

    def run_shell(
        self,
        cmd: str,
        cwd: Optional[str] = None,
    ) -> subprocess.CompletedProcess:
        """Execute a shell string command (for pipelines and complex scenarios).

        Args:
            cmd: Shell command string.
            cwd: Working directory for the subprocess (optional).

        Returns:
            CompletedProcess instance on success.

        Raises:
            ExternalToolError: If the command exits with a non-zero return code.
        """
        if self.verbose or self.dry_run:
            log.debug("SHELL CMD: %s", cmd)

        if self.dry_run:
            log.info("[dry-run] %s", cmd)
            return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

        log.debug("Running shell: %s", cmd)
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise ExternalToolError(
                tool="shell",
                cmd=[cmd],
                returncode=result.returncode,
                stderr=result.stderr,
            )

        if self.verbose and result.stdout:
            log.debug("STDOUT: %s", result.stdout.strip())

        return result
