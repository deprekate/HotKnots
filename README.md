# HotKnots
Fast, and best, DNA/RNA folding algorithm

This was my attempt to get an RNA/DNA secondary structure program that included pseudoknot
structures working as a Python C Extenstion.  

I tried several different algoritms/packages, and Hotknots was not only the fastest, but
also the only one that actually did pseudoknots.

I got it mostly working, but have not gottten around to making the parameter files as 
internal variables.  Having external data files within python PIP packages is a hassle
due to operating environments being virtual and not having actual file system paths.

To install:
```
pip install hotknots
```
or
```
git clone https://github.com/deprekate/HotKnots.git
pip install HotKnots/
```



To use on the command line:
```
echo AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC | hotknots.py
```
and the output should be the sequence, followed by the structure and mfe of the best folding
```
AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC
...((((((([[[[[.)))))))......]]]]]... -9.883
```


To import and use in other python code, you need to import the package, and then find out where
it is installed, so that it can find the various parameter files.  This is also when you can specify
which model and parameters to use:
```
import os
from hotknots import hotknots as hk
# initialize everything first
path = os.path.dirname(hk.__file__)
hk.initialize("DP", os.path.join(path, "parameters_DP09.txt" ) , os.path.join(path,"multirnafold.conf"), os.path.join(path,"pkenergy.conf") )

print(hk.fold("AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC", "DP"))
```
