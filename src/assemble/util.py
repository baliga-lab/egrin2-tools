import time
import math

def current_millis():
    """returns the current time in milliseconds"""
    return int(math.floor(time.time() * 1000))


def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [alist[i * length // wanted_parts: (i + 1) * length // wanted_parts]
            for i in range(wanted_parts)]
