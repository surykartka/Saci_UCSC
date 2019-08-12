import os

gff = 'combined_features.gff'
out = 'combined_features.bb'

id2name = {}
for line in open(gff):
	tab = line.strip().split('\t')
	if len(tab) < 8:
		continue
	desc = tab[8]
	desc_spl = desc.split(';')
	if desc.startswith('ID') and 'Name=' in desc:
		id2name[desc.split('ID=')[1].split(';')[0]] = desc.split('Name=')[1].split(';')[0]

with open(out, 'w') as f:
	for line in open(gff):
		tab = line.strip().split('\t')
		if len(tab) < 8:
			continue
		desc = tab[8]
		if tab[2] == 'CDS':
			col = '31,119,180'
			xid = desc.split('Parent=')[1].split(';')[0]
			print(tab[0], tab[3], tab[4], id2name[xid], 1000, tab[6], tab[3], tab[4], col, sep='\t', file=f)
		elif tab[2] in ('transcript', 'tRNA', 'rRNA'):
			col = '44,160,44'
			xid = desc.split('Parent=')[1].split(';')[0]
			print(tab[0], tab[3], tab[4], id2name[xid], 1000, tab[6], tab[3], tab[4], col, sep='\t', file=f)
		elif 'TSS' in tab[2]:
			col = '255,127,14'
			gene = desc.split('Gene=')[1].split(';')[0]
			print(tab[0], tab[3], tab[4], 'TSS', 1000, tab[6], tab[3], tab[4], col, sep='\t', file=f)
		elif 'terminator' in tab[2]:
			col = '214,39,40'
			gene = desc.split('Gene=')[1].split(';')[0]
			print(tab[0], tab[3], tab[4], 'term', 1000, tab[6], tab[3], tab[4], col, sep='\t', file=f)

os.system('sort -k1,1 -k2,2n %s > tmp2' % out)
os.system('bedToBigBed -type=bed9 -tab tmp2 chrom.sizes %s' % out)
os.system('rm tmp2')