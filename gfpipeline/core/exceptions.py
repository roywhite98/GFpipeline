"""Exception hierarchy for gfpipeline."""


class PipelineError(Exception):
    """Base class for all pipeline errors. Carries stage name and reason."""


class ConfigError(PipelineError):
    """Configuration file is missing required fields or has invalid format."""


class ToolNotFoundError(PipelineError):
    """A required external tool is not executable."""


class StageInputError(PipelineError):
    """A required input file for a stage does not exist."""


class ExternalToolError(PipelineError):
    """An external tool exited with a non-zero return code."""

    def __init__(self, tool: str, cmd: list, returncode: int, stderr: str = ""):
        self.tool = tool
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr
        msg = (
            f"Tool '{tool}' failed (exit {returncode}).\n"
            f"Command: {' '.join(str(c) for c in cmd)}\n"
            f"Stderr: {stderr}"
        )
        super().__init__(msg)


class ApiError(PipelineError):
    """NCBI API request failed."""

    def __init__(self, status_code: int | None = None, reason: str = ""):
        self.status_code = status_code
        self.reason = reason
        if status_code is not None:
            msg = f"API error (HTTP {status_code}): {reason}"
        else:
            msg = f"API error: {reason}"
        super().__init__(msg)
