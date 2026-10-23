from importlib.metadata import version

try:
    __version__ = version("compgenlnc")
except:
    __version__ = "Unknown"

__all__ = ["__version__"]