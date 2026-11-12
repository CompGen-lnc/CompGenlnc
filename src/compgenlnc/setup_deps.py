import importlib.util
import subprocess
import sys

import torch


def get_pyg_wheel_url() -> str:
    """Get the URL for pyg_lib wheel for the current PyTorch version."""
    torch_version = torch.__version__.split("+")[0]
    patches = torch_version.split(".")
    torch_main_version = f"{patches[0]}.{patches[1]}.0"

    if torch.cuda.is_available() and torch.version.cuda:
        cuda_tag = f"cu{torch.version.cuda.replace('.', '')}"
    else:
        cuda_tag = "cpu"

    base_url = "https://data.pyg.org/whl/"
    return f"{base_url}torch-{torch_main_version}+{cuda_tag}.html"


def install_pyg_lib() -> None:
    """Install pyg_lib for the current PyTorch version."""
    try:
        importlib.util.find_spec("pyg_lib")
    except ImportError:
        url = get_pyg_wheel_url()
        cmd = [sys.executable, "-m", "pip", "install", "pyg_lib", "-f", url]
        print(f"Ejecutando: {' '.join(cmd)}")
        subprocess.check_call(cmd)


def install_deps() -> None:
    """Install the rest of the dependencies for the package."""
    install_pyg_lib()
