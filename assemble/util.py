import time
import math

def current_millis():
    """returns the current time in milliseconds"""
    return int(math.floor(time.time() * 1000))
