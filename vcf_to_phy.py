import re
import gzip
import os
import argparse

parser = argparse.ArgumentParser(
		description="Turn VCF into phy.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

parser.add_argument('--infile', type=str, default=None,
			help='Gzipped infile to convert')

parser.add_argument('--thin', type=int, default = 1,
			help='Further thinning')

args = parser.parse_args()

# set up the input and output files
infile = args.infile
f = gzip.open(infile, 'r')
out_fa = re.sub('.vcf.gz', '_thin%s.phy' % args.thin, infile)
o = open(out_fa, 'w')

seq = {}

for ix, l in enumerate(f):
	l = l.decode('ascii')
	if re.search('CHROM', l):
		d = re.split('\t', l.rstrip())
		inds = d[9:]
		for ind in inds:
			seq[ind] = ''
	elif not re.search('^#', l):
		# thin it further!
		if ix % args.thin == 0:
			d = re.split('\t', l.rstrip())
	
			snps = {}
			alleles = [d[3]] + re.split(',', d[4])
			for ix, a in enumerate(alleles):
				snps[str(ix)] = a
			snps['.'] = 'N'

			genos = [re.search('^(\S)', x).group(1) for x in d[9:]]
		
			for ind, geno in zip(inds, genos):
				seq[ind] += snps[geno]
f.close()

o.write('%s %s\n' % (len(seq), len(seq[list(seq.keys())[0]])))
for ind, s in seq.items():
	o.write('%s %s\n' % (ind, s))
o.close()
