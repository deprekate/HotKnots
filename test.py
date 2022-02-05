import faulthandler; faulthandler.enable()
import sys

import HotKnots as hk

model = "CC"
params = "params/parameters_CC09.txt"
hk.initialize(model, params)


print(hk.fold("UUUGCCCUGAAACUGGCGCGUGAGAUGGGGCGACCCGACUGGCGUGCCAU", model))

print(hk.fold("AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC", model))

print(hk.fold("GGCGCGTGAGATGGGGCGACCCGACTGGCGTGCCA", model))

