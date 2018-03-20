import re
import gzip
import os

dir = '/home/ubuntu/'

in_vcf = os.path.join(dir, 'Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.vcf.gz')
f = gzip.open(in_vcf, 'r')

out_g = os.path.join(dir, 'Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.geno')
out_s = os.path.join(dir, 'Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.snp')
out_i = os.path.join(dir, 'Massoko_Dryad_VCF_final_subset_noIndels_maf05_thinned1K.ind')
o1 = open(out_g, 'w')
# total hack the snp data because human data
o2 = open(out_s, 'w')
o3 = open(out_i, 'w')

alleles = {	'0|0': '0', '0|1': '1', '1|1': '2', 
		'0/0': '0', '0/1': '1', '1/1': '2', '1|0': '1'}

gen_length = 850000000
max_len = 0
cur_len = 0
cur_chr = 'NA'

snp = 1
snps = {}
for l in f:
	l = l.decode('ascii')
	if re.search('CHROM', l):
		d = re.split('\t', l.rstrip())
		inds = d[9:]
		for ind in inds:
			group = re.sub('_[^_]+$', '', ind)
			o3.write('%s\tU\t%s\n' % (ind, group))
		o3.close()
	elif not re.search('^#', l):
		d = re.split('\t', l.rstrip())
	
		if d[0] != cur_chr:
			cur_chr = d[0]	
			cur_len += max_len
	
		max_len = int(d[1])
		cur_chr = d[0]

		genos = [re.search('^(\S\S\S)', x).group(1) for x in d[9:]]
		if len(set(genos)) > 1:
			genos = [alleles[x] if x in alleles else '9' for x in genos] 		
			o1.write('%s\n' % ''.join(genos))
			o2.write('rs%s\t1\t%.3f\t%s\n' % 
				(snp, (int(d[1]) + cur_len) / float(gen_length), int(d[1]) + cur_len))
			snp += 1
f.close()
o1.close()
