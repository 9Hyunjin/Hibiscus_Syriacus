import pandas as pd

df = pd.read_csv("Orthogroups_withDAG.csv", sep="\t")

subgenomeA = {1, 3, 9, 14, 10, 13, 17, 2, 30, 31, 27, 34, 15, 25, 38, 41, 42, 22, 28, 29, 16, 19, 24}
subgenomeB = {6, 4, 8, 23, 33, 18, 21, 26, 36, 40, 32, 39, 5, 20, 35, 37, 7, 11, 12}

def countSubgenome(genes):
	if pd.isna(genes):
		return pd.Series([0, 0, 0])
	geneList = genes.split(", ")
	total = len(geneList)
	aCount = 0
	bCount = 0
	for gene in geneList:
		try:
			chromNum = int(gene.split("_")[1].split("g")[0])
			if chromNum in subgenomeA:
				aCount += 1
			elif chromNum in subgenomeB:
				bCount += 1
		except:
			continue
	return pd.Series([total, aCount, bCount])

dfCounts = df["Baekdansim"].apply(countSubgenome)
dfCounts.columns = ["totalGeneCount", "subgenomeAGeneCount", "subgenomeBGeneCount"]

dfResult = pd.concat([df["Orthogroup"], dfCounts], axis=1)
dfResult = dfResult[dfResult["totalGeneCount"] > 0]
dfResult.to_csv("BaekdansimCountsFiltered.csv", index=False)

print(dfResult.head())

mainDf = pd.read_csv("BaekdansimCountsFiltered.csv")
flDagDf = pd.read_csv("DAG_TotalFL_Domain.csv", header=None)
crgDagDf = pd.read_csv("DAG_CRG_Domain.csv")
drDagDf = pd.read_csv("DAG_Drought_Domain.csv")
etDagDf = pd.read_csv("DAG_Ethylene_Domain.csv")

flDags = set(flDagDf.iloc[:, -1].astype(str).str.strip().unique())
crgDags = set(crgDagDf["DAG"].astype(str).str.strip().unique())
drDags = set(drDagDf["DAG"].astype(str).str.strip().unique())
etDags = set(etDagDf["DAG"].astype(str).str.strip().unique())

def filterByDags(df, dagSet):
	def extractDag(x):
		parts = x.split("DAG")
		if len(parts) == 2 and parts[1].isdigit():
			return "DAG" + parts[1]
		return None

	return df[df["Orthogroup"].apply(lambda x: extractDag(x) in dagSet)]

flResult = filterByDags(mainDf, flDags)
crgResult = filterByDags(mainDf, crgDags)
drResult = filterByDags(mainDf, drDags)
etResult = filterByDags(mainDf, etDags)
rabResult = mainDf[mainDf["Orthogroup"].str.contains("DAG4252")]

rabResult.to_csv("Baekdansim_RabGTPase_SubgenomeCount.csv", index=False)
flResult.to_csv("Baekdansim_FL_SubgenomeCount.csv", index=False)
crgResult.to_csv("Baekdansim_CRG_SubgenomeCount.csv", index=False)
drResult.to_csv("Baekdansim_Drought_SubgenomeCount.csv", index=False)
etResult.to_csv("Baekdansim_Ethylene_SubgenomeCount.csv", index=False)

