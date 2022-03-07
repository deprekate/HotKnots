# HotKnots
Fast, and best, DNA/RNA folding algorithm


This was my attempt to get an RNA/DNA secondary structure program that included pseudoknot
structures working as a Python C Extenstion.  

I tried several different algoritms/packages, and Hotknots was not only the fastest, but
also the only that actually did pseudoknots.

I got it mostly working, but have not gottten around to making the parameter files as 
internal variables.  Having extenal data files with python PIP pacakges is a hassle
due to operating environments being virtual and not having actual file system paths.
Currently the params folder has to be in the same folder that you run HotKnots

To install:
```
git clone https://github.com/deprekate/HotKnots.git
pip install HotKnots/ --user
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
it is installed, so it can find the various parameter files.  This is also when you can specify
which model and paramters to use:
```
import HotKnots as hk
# initialize everything first
params = os.path.dirname(hk.__file__)
hk.initialize( "DP", os.path.join(params,"parameters_DP09.txt") , os.path.join(params,"multirnafold.conf"), os.path.join(params,"pkenergy.conf") )

print(hk.fold("AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC", "DP"))
```
