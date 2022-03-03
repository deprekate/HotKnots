#!/usr/bin/env python3
import io
import os
import sys
import gzip
import faulthandler; faulthandler.enable()
import sys
import fileinput
import argparse
from argparse import RawTextHelpFormatter
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 


def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def read_fasta(filepath, base_trans=str.maketrans('','')):
    contigs_dict = dict()
    name = seq = ''

    lib = gzip if filepath.endswith(".gz") else io
    with lib.open(filepath, mode="rb") as f:
        for line in f:
            if line.startswith(b'>'):
                contigs_dict[name] = seq
                name = line[1:].decode("utf-8").split()[0]
                seq = ''
            else:
                seq += line.decode("utf-8").rstrip().lower()
        contigs_dict[name] = seq.translate(base_trans)
    if '' in contigs_dict: del contigs_dict['']

    assert contigs_dict, "No DNA sequence found in the infile"

    return contigs_dict

def rev_comp(seq):
    seq_dict = {'a':'t','t':'a','g':'c','c':'g',
				'A':'T','T':'A','G':'C','C':'G',
                'n':'n',
                'r':'y','y':'r','s':'s','w':'w','k':'m','m':'k',
                'b':'v','v':'b','d':'h','h':'d'}
    return "".join([seq_dict[base] for base in reversed(seq)])

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	#parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-w', '--window_size', type=int, default=30, help='The size of the window [30]')
	args = parser.parse_args()

	contigs = read_fasta(args.infile)

	for header,sequence in contigs.items():
		sequence = sequence.upper().replace('M','A').replace('S','A').replace('R','A')
		for i in range(0, len(sequence)-args.window_size-2, 3):
			#print(sequence[i:i+args.window_size])
			print(rev_comp(sequence[i:i+args.window_size]))





