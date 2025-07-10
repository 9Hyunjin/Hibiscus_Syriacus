import pandas as pd

def loadDagDict(dagFile, speciesList):
	dagDf = pd.read_csv(dagFile, sep='\t')
	dagDf = dagDf.set_index('DAG')
	dagDict = {}
	for specie in speciesList:
		for dag, genes in zip(dagDf.index, dagDf[specie]):
			for gene in str(genes).split(','):
				gene = gene.strip()
				if gene and gene != 'nan':
					dagDict[gene] = dag
	return dagDict

def splitOrthogroupsByDag(orthogroupFile, dagDict, speciesList, outputFile):
	orthogroupsDf = pd.read_csv(orthogroupFile, sep='\t')
	resultRows = []

	for _, row in orthogroupsDf.iterrows():
		newOrthogroups = {}
		for specie in speciesList:
			genes = str(row[specie]).split(',')
			for gene in genes:
				gene = gene.strip()
				if gene:
					dag = dagDict.get(gene, 'DAGNO')
					newOrthogroup = row['Orthogroup'] + dag
					if newOrthogroup not in newOrthogroups:
						newOrthogroups[newOrthogroup] = {s: [] for s in speciesList}
					newOrthogroups[newOrthogroup][specie].append(gene)

		for newOrthogroup, data in newOrthogroups.items():
			newRow = {'Orthogroup': newOrthogroup}
			for specie, geneList in data.items():
				if geneList:
					newRow[specie] = ', '.join(geneList)
				else:
					newRow[specie] = 'NA'
			resultRows.append(newRow)

	resultDf = pd.DataFrame(resultRows, columns=['Orthogroup'] + speciesList)
	resultDf.to_csv(outputFile, index=False, sep='\t')

def main():
	speciesList = ['Baekdansim', 'Gangneung', 'Gossypium_herbaceum', 'Gossypium_hirsutum', 'Gossypium_raimondii', 'Theobroma_cacao']
	dagDict = loadDagDict('DAG_ortho.tsv', speciesList)
	splitOrthogroupsByDag('Orthogroups.tsv', dagDict, speciesList, 'Orthogroups_withDAG.csv')

if __name__ == "__main__":
	main()

