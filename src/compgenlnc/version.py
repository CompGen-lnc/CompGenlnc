from importlib.metadata import version

try:
    __version__ = version("compgenlnc")
except Exception:
    __version__ = "Unknown"

__all__ = ["__version__"]