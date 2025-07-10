import pandas as pd
import numpy as np

subgenomeMap = {'Triplet set1_A': {1, 3}, 'Triplet set1_B': {6}, 'Triplet set2_A': {9, 14},
	'Triplet set2_B': {4}, 'Triplet set3_B': {8, 23, 33},
	'Triplet set4_A': {10, 13, 17}, 'Triplet set5_B': {18, 21, 26}, 'Triplet set6_A': {2, 30, 31},
	'Quadruplet set1_A': {27, 34}, 'Quadruplet set1_B': {36, 40}, 'Quadruplet set2_A': {15, 25},
	'Quadruplet set2_B': {32, 39}, 'Quintuplet set1_A': {38, 41, 42},
	'Quintuplet set1_B': {5, 20}, 'Quintuplet set2_A': {22, 28, 29},
	'Quintuplet set2_B': {35, 37}, 'Sextuplet set1_A': {16, 19, 24},
	'Sextuplet set1_B': {7, 11, 12}}

crossPairMap = {'Triplet set1 A vs B': ('Triplet set1_A', 'Triplet set1_B'),
	'Triplet set 2 A vs B': ('Triplet set2_A', 'Triplet set2_B'),
	'Quadruplet set 1 A vs B': ('Quadruplet set1_A', 'Quadruplet set1_B'),
	'Quadruplet set 2 A vs B': ('Quadruplet set2_A', 'Quadruplet set2_B'),
	'Quintuplet set 1 A vs B': ('Quintuplet set1_A', 'Quintuplet set1_B'),
	'Quintuplet set 2 A vs B': ('Quintuplet set2_A', 'Quintuplet set2_B'),
	'Sextuplet 1 A vs B': ('Sextuplet set1_A', 'Sextuplet set1_B')}

colDf = []
with open('Baek_Baek.collinearity', 'r') as file:
	currentBlock = None
	for line in file:
		line = line.strip()
		if line.startswith("# Alignment"):
			currentBlock = line.split()[2].strip(":")
		elif line and not line.startswith('#'):
			parts = line.split()
			id1, id2 = parts[0], parts[2]
			colDf.append((currentBlock, id1, id2))
colDf = pd.DataFrame(colDf, columns=["blockId", "id1", "id2"])

ksDf = pd.read_csv('Baek_Baek.collinearity.ks', sep='\t')
mergedDf = pd.merge(ksDf, colDf, on=["id1", "id2"], how='inner')

def extractChr(geneId):
	return int(geneId.split('_')[1].split('g')[0])
mergedDf['Chr1'] = mergedDf['id1'].apply(extractChr)
mergedDf['Chr2'] = mergedDf['id2'].apply(extractChr)

blockKsIntra = []
for name, chrs in subgenomeMap.items():
	df = mergedDf[(mergedDf['Chr1'].isin(chrs)) & (mergedDf['Chr2'].isin(chrs))]
	avgKs = df.groupby('blockId')['ks_YN00'].mean().reset_index()
	avgKs['Group'] = name
	blockKsIntra.append(avgKs)
blockKsIntra = pd.concat(blockKsIntra, ignore_index=True)
blockKsIntra = blockKsIntra[blockKsIntra['ks_YN00'] > 0]

blockKsInter = []
for label, (grp1, grp2) in crossPairMap.items():
	chrs1 = subgenomeMap[grp1]
	chrs2 = subgenomeMap[grp2]
	df = mergedDf[((mergedDf['Chr1'].isin(chrs1)) & (mergedDf['Chr2'].isin(chrs2))) |
		((mergedDf['Chr1'].isin(chrs2)) & (mergedDf['Chr2'].isin(chrs1)))]
	avgKs = df.groupby('blockId')['ks_YN00'].mean().reset_index()
	avgKs['Group'] = label
	blockKsInter.append(avgKs)
blockKsInter = pd.concat(blockKsInter, ignore_index=True)
blockKsInter = blockKsInter[blockKsInter['ks_YN00'] > 0]

