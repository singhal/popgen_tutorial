import re
import gzip
import os
import argparse

parser = argparse.ArgumentParser(
                description="Summarize LD file.",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
                )

parser.add_argument('--infile', type=str, default=None,
                        help='VCFtools LD file to convert')

parser.add_argument('--win', type=int, default = 1,
                        help='The window across which to summarize LD patterns')

args = parser.parse_args()

f = open(args.infile, 'r')
out = args.infile + '_summary_window%s' % args.win
o = open(out, 'w')

r2 = {}
header = f.next()
for l in f:
	d = re.split('\s+', l.rstrip())
	dist = int(d[2]) - int(d[1])
	dist = int(dist / float(args.win)) * args.win

	if dist not in r2:
		r2[dist] = {'val': 0, 'num': 0}
	r2[dist]['val'] += float(d[4])
	r2[dist]['num'] += 1
f.close()

o.write('distance\tnum_comparisons\tr2\n')
for dist in sorted(r2.keys()):
	avg = r2[dist]['val'] / float(r2[dist]['num'])
	o.write('%s\t%s\t%.4f\n' % (dist, r2[dist]['num'], avg))
o.close()
	
