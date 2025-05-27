import psutil
import os
import logging

def memory_usage_mb():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024**2

def log_memory(tag=""):
    mem = memory_usage_mb()
    logging.info(f"{tag} | Memory usage: {mem:.2f} MB")