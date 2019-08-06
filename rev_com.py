import argparse
import re

# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
	help="name of the fasta file")
ap.add_argument("-c", "--choice", required=True,
	help="0 for reverse; 1 for complement; 2 for reverse complement")
ap.add_argument("-o", "--output", required=True,
	help="name of output fasta file")
args = vars(ap.parse_args())

def seq_reverse(file_read, file_write):
	fh = open(file_read, 'r')
	i=0
	fh2=open(file_write, 'w')
	rev_tmp=''
	head=''
	for line in fh:
		#print(line)
		if re.search('>', line):
			#print(head+rev_tmp[-2::-1])
			if i>1:
				fh2.write(head+rev_tmp[-2::-1]+'\n')
			rev_tmp=''
			head=line
			start=1
		else:
			rev_tmp=rev_tmp+line
		i+=1
	else:
		#print(head+rev_tmp[-2::-1])
		fh2.write(head+rev_tmp[-2::-1]+'\n')

def seq_complement(file_read, file_write):
	fh = open(file_read, 'r')
	i=0
	fh2=open(file_write, 'w')
	com_tmp=''
	head=''
	for line in fh:
		#print(line)
		if re.search('>', line):
			if i>1:
				#print(head+com_tmp.lstrip())
				fh2.write(head+com_tmp)				
			com_tmp=''
			head=line
			start=1
		else:
			for ch in line:
				if ch=='A' or ch=='a':
					com_tmp=com_tmp+'T'
				elif ch=='T' or ch=='t':
					com_tmp=com_tmp+'A'
				elif ch=='G' or ch=='g':
					com_tmp=com_tmp+'C'
				elif ch=='C' or ch=='c':
					com_tmp=com_tmp+'G'
				else:
					com_tmp=com_tmp+ch
			#rev_tmp=rev_tmp+line
		i+=1
	else:
		#print(head+com_tmp.lstrip())
		fh2.write(head+com_tmp)

#print(args['option'])
#print(args['choice'])
if args['choice']=='0':
	print('Outputting the reverse of sequences.')
	seq_reverse(args['input'], args['output'])
elif args['choice']=='1':
	print('Outputting the complement of sequences.')
	seq_complement(args['input'], args['output'])
elif args['choice']=='2':
	print('Outputting the reverse complement of sequences.')
	seq_reverse(args['input'], 'tmp')
	seq_complement('tmp', args['output'])
	#system('rm tmp')
else :
	print('ERROR : Check options')
#seq_reverse(fh, 'out.fa')
#seq_complement(fh, 'out2.fa')

#seq_complement(args['input'], args['output'])