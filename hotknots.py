#!/usr/bin/env python3
import os
import sys
import fileinput
import argparse
from argparse import RawTextHelpFormatter
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import fileinput
from termcolor import colored
from importlib import reload

sys.path.pop(0)
from hotknots import hotknots as hk


def test():
	print(hk.fold("UUUGCCCUGAAACUGGCGCGUGAGAUGGGGCGACCCGACUGGCGUGCCAU", model))
	print(hk.fold("AACCCCUGCUGAAUAAAGCGGGGAAUAACUAUUCUAC", model))
	print(hk.fold("GGCGCGTGAGATGGGGCGACCCGACTGGCGTGCCA", model))

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	#parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-m', '--model', type=str, default='DP', choices=['DP', 'CC', 'RE'], help='The model to use [DP]')
	#parser.add_argument('-p', '--params', required=True, type=is_valid_file, help="this is the path to the parameter file")
	#parser.add_argument('-1', '--config_file', required=True, type=is_valid_file, help="this is the path to the config file")
	#parser.add_argument('-2', '--config_filepk', required=True, type=is_valid_file, help="this is the path to the configPK file")
	args = parser.parse_args()


	# initialize everything first
	params = os.path.dirname(hk.__file__)
	#for line in fileinput.input():
	for line in args.infile:
		# then run each sequence through
		print(line.rstrip())
		hk.initialize( args.model, os.path.join(params,"parameters_DP09.txt") , os.path.join(params,"multirnafold.conf"), os.path.join(params,"pkenergy.conf") )
		seq,mfe = hk.fold( line.rstrip().upper() , args.model )
		print(seq, colored(mfe, 'green') )

