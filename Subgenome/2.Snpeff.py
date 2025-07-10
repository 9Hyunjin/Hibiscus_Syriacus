import pandas as pd
import glob
import os
import numpy as np
from scipy.stats import chi2_contingency
from collections import defaultdict
from statsmodels.stats.multitest import multipletests

impactEffects = {'missense_variant', 'stop_gained', 'stop_lost', 'start_lost', 'start_gained',
	'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'exon_loss_variant'}

for file in glob.glob('Snpeff_DAG*.vcf'):
	if '_filtered' not in file:
		outFile = file.replace('.vcf', '_filtered.vcf')
		with open(file, 'r') as f:
			header = []
			lines = []
			for line in f:
				if line.startswith('#'):
					header.append(line)
				else:
					if 'ANN=' in line:
						infoField = [x for x in line.strip().split('\t') if 'ANN=' in x]
						if not infoField:
							continue
						annField = infoField[0].split('ANN=')[-1].split(';')[0]
						annotations = annField.split(',')
						for ann in annotations:
							parts = ann.split('|')
							if len(parts) > 1 and parts[1] in impactEffects:
								lines.append(line)
								break
		with open(outFile, 'w') as out:
			for h in header:
				out.write(h)
			for l in lines:
				out.write(l)

sampleMeta = pd.read_csv('Sample_Info.csv')
sampleToGroup = dict(zip(sampleMeta['SampleName'], sampleMeta['Region']))

effectSet = {'missense_variant', 'stop_gained', 'stop_lost', 'start_lost',
	'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant'}

for file in glob.glob('Snpeff_DAG*_filtered.vcf'):
	records = []
	with open(file) as f:
		for line in f:
			if line.startswith('#'):
				if line.startswith('#CHROM'):
					header = line.strip().split('\t')
					samples = header[9:]
				continue
			cols = line.strip().split('\t')
			chrom, pos, ref, alt, info = cols[0], int(cols[1]), cols[3], cols[4], cols[7]
			genotypes = cols[9:]
			if 'ANN=' not in info:
				continue
			annPart = info.split('ANN=')[1].split(';')[0]
			annotations = annPart.split(',')
			genes = set()
			types = set()
			for ann in annotations:
				fields = ann.split('|')
				if len(fields) > 4:
					if fields[1] in effectSet:
						types.add(fields[1])
						genes.add(fields[3])
			if not genes:
				continue
			for gene in genes:
				groupGenotypes = defaultdict(lambda: {'hom_ref': 0, 'het': 0, 'hom_alt': 0})
				for gt, sample in zip(genotypes, samples):
					group = sampleToGroup.get(sample)
					if not group or gt == './.':
						continue
					gt = gt.replace('|', '/')
					alleles = gt.split('/')
					if alleles == ['0', '0']:
						groupGenotypes[group]['hom_ref'] += 1
					elif set(alleles) == {'0', '1'}:
						groupGenotypes[group]['het'] += 1
					elif alleles == ['1', '1']:
						groupGenotypes[group]['hom_alt'] += 1
				if len(groupGenotypes) >= 2:
					table = []
					for g in groupGenotypes:
						counts = groupGenotypes[g]
						table.append([counts['hom_ref'], counts['het'], counts['hom_alt']])
					try:
						chi2, pval, _, _ = chi2_contingency(table)
						records.append({'Gene': gene, 'CHROM': chrom, 'POS': pos,
							'Effect': ','.join(types), 'PValue': pval})
					except:
						continue
	df = pd.DataFrame(records)
	if not df.empty:
		df['-log10(P)'] = df['PValue'].apply(lambda x: -np.log10(x) if x > 0 else 0)
		outName = file.replace('.vcf', '_genotypelevel_stats_chisq.csv')
		df.to_csv(outName, index=False)

for file in glob.glob('*_genotypelevel_stats_chisq.csv'):
	df = pd.read_csv(file)
	if 'PValue' not in df.columns:
		continue
	fdr_results = multipletests(df['PValue'], method='fdr_bh')
	df['FDR'] = fdr_results[1]
	significant = df[df['FDR'] < 0.01]
	outName = file.replace('.csv', '_fdr_sig.csv')
	significant.to_csv(outName, index=False)

