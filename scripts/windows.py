import io
import sys
import gzip

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
                'n':'n',
                'r':'y','y':'r','s':'s','w':'w','k':'m','m':'k',
                'b':'v','v':'b','d':'h','h':'d'}
    return "".join([seq_dict[base] for base in reversed(seq)])

contigs = read_fasta(sys.argv[1])

for header,sequence in contigs.items():
	for i in range(0, len(sequence), 3):
		print(sequence[i:i+50])
		#print(rev_comp(sequence[i:i+50]))
