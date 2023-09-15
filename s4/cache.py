"""Defining cache folders for storing entries from the Materials Project."""
import os

CACHE_DIR = os.path.join(os.path.expanduser("~"), ".s4_cache")

if not os.path.exists(CACHE_DIR):
    os.makedirs(CACHE_DIR, exist_ok=True)
