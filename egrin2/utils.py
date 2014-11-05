import numpy as np
import os
import pandas as pd
import weblogolib as wl
import cStringIO
import sqlite3

from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def system( cmd ):
    print cmd
    tmp = os.popen( cmd ).read()
    return tmp

def stop():
    raise 'STOPPING ON PURPOSE!'

