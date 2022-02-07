# HotKnots
Fast, and best, DNA/RNA folding algorithm


This was my attempt to get an RNA/DNA secondary structure program that included pseudoknot
structures working as a Python C Extenstion.  

I tried several different algoritms/packages, and Hotknots not only was the fastest, but
also the only that actually did pseudoknots.

I got it mostly working, but have not gottten around to making the parameter files as 
internal variables.  Having extenal data files with python PIP pacakges is a hassle
due to operating environments being virtual and not having actual file system paths.

To install:
```
https://github.com/deprekate/HotKnots.git
pip install HotKnots --user
```



To use on the command line:
```
echo AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC | python3 hotknots.py -p params/parameters_DP09.txt
```
and the output should be the sequence, followed by the structure and mfe of the best folding
```
AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC
...((((((([[[[[.)))))))......]]]]]... -9.883
```


To import and use in other python code:
```
import HotKnots as hk

model = "CC"
params = "params/parameters_CC09.txt"
hk.initialize(model, params)

print(hk.fold("AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC", model))
```
