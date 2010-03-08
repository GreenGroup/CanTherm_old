"""
Copyright (c) 2002-2009 William H. Green and the CanTherm Team

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


#!/usr/bin/env python

import os
import sys
from Harmonics import *

file = sys.argv[1]


Kcos=[]
Ksin =[]
harmonic = Harmonics(5,Kcos,Ksin)
harmonic.fitPotential(file)

for i in range(5):
    print harmonic.Kcos[i], harmonic.Ksin[i]

for i in range(5):
    print harmonic.Kcos[i], " ", "%3.1f"%(float(i+1))
for i in range(5):
    print harmonic.Ksin[i], " ", "%3.1f"%(float(i+1))




