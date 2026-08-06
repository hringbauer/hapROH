import logging
import psutil

def print_memory_usage(logger:None|logging.Logger=None):
    process = psutil.Process()
    if logger is not None:
        logger.debug(f"Memory usage: {process.memory_info().rss:,}")
    else:
        print(f"Memory usage: {process.memory_info().rss:,}")