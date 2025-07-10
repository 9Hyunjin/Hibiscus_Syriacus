import os
import glob
import pandas as pd

def loadGeneSymbols(symbolFile):
	with open(symbolFile, 'r') as inf:
		inf.readline()
		return {line.strip().split('\t')[0]: line.strip().split('\t')[2] for line in inf}

def loadRepresentativeArch(repFile):
	representative = {}
	with open(repFile, 'r') as rep:
		for line in rep:
			sLine = line.strip().split('\t')
			if sLine[1] != 'N/A' or sLine[2] != 'N/A':
				key = '\t'.join(sLine[1:])
				representative.setdefault(key, set()).add(sLine[0].split('.')[0])
	return {k: list(v) for k, v in representative.items()}

def processSpeciesFiles(representative, geneSymbols):
	for file in glob.glob('*Filled*.txt'):
		species = file.split('Filled_')[1].split('.')[0]
		if species == 'Arabidopsis_thaliana':
			continue

		with open(file, 'r') as inf, \
			open(f'{species}_process.csv', 'w') as outFull, \
			open(f'{species}_process_noDomain.csv', 'w') as outNoDomain, \
			open(f'{species}_process_notmatched_with_A.thaliana.csv', 'w') as outNoMatch:

			header = 'GeneID,Domain architecture(Domain + Family),Domain architecture(Domain only),Arabidopsis_Gene_Identifier,Gene_symbol\n'
			outFull.write(header)
			outNoDomain.write(header)
			outNoMatch.write('GeneID,Domain architecture(Domain + Family),Domain architecture(Domain only)\n')

			for line in inf:
				sLine = line.strip().split('\t')
				key = '\t'.join(sLine[1:])
				if sLine[1] == 'N/A' and sLine[2] == 'N/A':
					outNoDomain.write(f"{sLine[0]},{','.join(sLine[1:])},{','.join(sLine[1:])}\n")
				elif key in representative:
					genes = representative[key]
					symbols = [geneSymbols[g] for g in genes]
					if len(genes) > 1:
						outFull.write(f"{sLine[0]},{','.join(sLine[1:])},{'|'.join(genes)},{'|'.join(symbols)}\n")
					else:
						outFull.write(f"{sLine[0]},{','.join(sLine[1:])},{genes[0]},{symbols[0]}\n")
				else:
					outNoMatch.write(f"{sLine[0]},{','.join(sLine[1:])}\n")

def assignDagIdentifiers(fileList):
	dagDict = {}
	def assignDag(row):
		key = (row['Domain architecture(Domain + Family)'], row['Domain architecture(Domain only)'])
		if key not in dagDict:
			dagDict[key] = f"DAG{len(dagDict)+1}"
		return dagDict[key]

	for file in fileList:
		data = pd.read_csv(file)
		data['DAG'] = data.apply(assignDag, axis=1)
		outFile = file.replace('_process.csv', '_DAG.csv')
		data[['GeneID', 'Domain architecture(Domain + Family)', 'Domain architecture(Domain only)', 'DAG']].to_csv(outFile, index=False)

def combineDagFiles():
	inputFiles = glob.glob('*_DAG.csv')
	allData = []
	for file in inputFiles:
		species = os.path.splitext(file)[0].split("_DAG")[0]
		speciesName = species.replace('Hibiscus_syriacus_', '') if 'Hibiscus_syriacus' in species else species
		df = pd.read_csv(file, names=["GeneID", "Domain_architecture_Family", "Domain_architecture_Domain", "DAG"], header=0)
		df["Species"] = speciesName
		allData.append(df)

	combined = pd.concat(allData)
	groupGene = combined.groupby(["DAG", "Species"])["GeneID"].apply(lambda x: ','.join(x)).reset_index()
	groupGeneCount = combined.groupby(["DAG", "Species"])["GeneID"].count().reset_index(name="count")

	groupGene.pivot(index="DAG", columns="Species", values="GeneID").fillna("").reset_index().to_csv("DAG_ortho.tsv", sep="\t", index=False)
	groupGeneCount.pivot(index="DAG", columns="Species", values="count").fillna(0).astype(int).reset_index().to_csv("DAG_ortho_geneCount.tsv", sep="\t", index=False)

def main():
	geneSymbols = loadGeneSymbols('Arabidopsis_thaliana.TAIR10.pep.representative.fa.gene.symbol.txt')
	representative = loadRepresentativeArch('Filled_Arabidopsis_thaliana.all.domain_architecture.txt')
	processSpeciesFiles(representative, geneSymbols)

	processFiles = glob.glob('*_process.csv')
	assignDagIdentifiers(processFiles)
	combineDagFiles()

if __name__ == "__main__":
	main()

