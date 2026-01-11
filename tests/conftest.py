import os

# Avoid joblib multiprocessing warnings in test environment.
os.environ.setdefault("JOBLIB_MULTIPROCESSING", "0")
