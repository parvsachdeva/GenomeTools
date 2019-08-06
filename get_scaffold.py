import argparse
import re

# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
	help="name of the fasta file")
ap.add_argument("-s", "--scaffold", required=True,
	help="name of the scaffold to fetch")
ap.add_argument("-o", "--output", required=True,
	help="name of output fasta file")
args = vars(ap.parse_args())

def get_scaffold(file_read, file_write, scaf):
	fh1 = open(file_read, 'r')
	i=0
	fh2=open(file_write, 'w')
	tmp=''; chk=0
	for line in fh1:
		if re.search('>', line):
			chk=0
		if re.search(scaf, line):
			chk=1
		if chk==1:
			tmp=tmp+line
	fh2.write(tmp)
#######


if args['scaffold']:
	print('Outputting the scaffold.')
	get_scaffold(args['input'], args['output'], args['scaffold'])
else :
	print('ERROR : Check options')
#seq_reverse(fh, 'out.fa')
#seq_complement(fh, 'out2.fa')

#seq_complement(args['input'], args['output'])