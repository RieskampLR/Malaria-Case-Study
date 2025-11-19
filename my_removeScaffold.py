# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 13:32:33 2025

@author: lea
"""

import sys
from pathlib import Path
import re

genome_file = Path(sys.argv[1])
ouput_file = Path(sys.argv[2])

seq_dictionary = {}

with open (genome_file, "r"):
    for line in genome_file:
        if "length=" in line:
            words = re.split(r'[\r=]', line)
            if words [3] > 3000:
                seq_dictionary[line] = line + 1
    

with open (output_file, "w"):
    for i in seq_dictionary:
        output_file.write(i)
