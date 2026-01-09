import subprocess
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def run_command(cmd, description: str, stdout=None, stderr=subprocess.PIPE):
    """
    Run a shell command with error handling and logging.
    
    Args:
        cmd: List of command arguments (e.g., ["samtools", "view", "file.bam"])
        description: Human-readable description of the command
        stdout: Where to redirect stdout (default: None = inherit)
        stderr: Where to redirect stderr (default: PIPE)
    
    Raises:
        subprocess.CalledProcessError: If command fails
    """
    cmd_str = " ".join(str(c) for c in cmd)
    logger.info(f"Running: {description} (command: {cmd_str})")
    
    try:
        result = subprocess.run(
            cmd,
            stdout=stdout,
            stderr=stderr,
            check=True,
            text=True
        )
        if stderr and result.stderr:
            logger.debug(f"Command stderr: {result.stderr}")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {description}")
        logger.error(f"Error output: {e.stderr}")
        raise

def validate_inputs(paths: list, types: list):
    """
    Validate input paths are of the correct type (file/dir).
    
    Args:
        paths: List of paths to validate
        types: List of types ("file" or "dir") corresponding to paths
    
    Raises:
        FileNotFoundError: If path does not exist
        ValueError: If path is not the correct type
    """
    if len(paths) != len(types):
        raise ValueError("Paths and types lists must have the same length")
    
    for path, typ in zip(paths, types):
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Path not found: {path}")
        
        if typ == "file" and not path.is_file():
            raise ValueError(f"Expected file, got directory: {path}")
        elif typ == "dir" and not path.is_dir():
            raise ValueError(f"Expected directory, got file: {path}")
