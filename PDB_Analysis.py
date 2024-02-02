from Bio.PDB import Selection
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import parse_pdb_header
from Bio import Align
from Bio.PDB import NeighborSearch
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import random
import shutil


# this function gives out the valid chain and the corresponding amino acid sequence 
#input: .ent pdb file
#output: a dictionary{validChainId: aa sequence} and a dictionary for original chain info
def chain_basic_info(pdbfile):
	print('1 works')
	parser = PDBParser(PERMISSIVE=1)
	print('2works')
	structure_id = str(pdbfile)
	print('3works')
	filename = str(pdbfile)
	print('4works')
	structure = parser.get_structure(structure_id, filename)
	print('5works')
	print(structure)

	#get the amino acids sequence for each chain
	model = structure[0]
	residueTotal = [] #each element in this list is a long string, which stands for all the residues in one chain
	chainList = []#a list containing the identities of all chains 
	originalChainDic = {} #a dictionary, the key is the chain, the corresponding value is its amino acid sequence.
	aaDic = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','LYS':'K','ARG':'R',
						'HIS':'H','ASP':'D',"ASN":'N','GLU':'E','GLN':'Q','PRO':'P','TRP':'W','PHE':'F',
						'TYR':'Y','MET':'M','CYS':'C'}
	for model in structure:
		for chain in model:
			chainList.append(chain) 
			resChain = str() #store the residues from the same chains, e,g. all residues from chain A
			for residue in chain:
				res = str(residue.get_resname()) #get the name for each chain
				if res in aaDic: # we only consider residues that are amino acids
					res = aaDic[res] #change the amino acid names into their abbreviation
					resChain += res
			residueTotal.append(resChain)

	for i, chain in enumerate(chainList):
		originalChainDic[chainList[i]] = residueTotal[i]

	validChainDic = {} #store the id of valid chains as the key, and the corresponding value is its amino acid sequence
	for i in range(len(chainList)):
		if originalChainDic[chainList[i]] == '' or len(originalChainDic[chainList[i]]) < 20:
			continue
		else:
			validChainDic[chainList[i]] = originalChainDic[chainList[i]]
	return validChainDic, originalChainDic

#this function give out the similar chains list in a list:[[similar chains]]
#input: a dictionary contains chain id and its corresponding sequence information
#conditions: if two chains share more that 95% of their amino acids, these two chains are regarded as similar
def similar_chain(ChainDictionary):
	chainList = list(ChainDictionary.keys()) #get all the valid chains 
	aligner = Align.PairwiseAligner()
	aligner.match_score = 1.0
	similarChainList = []
	alignedChain = []#store the chains that have already been aligned
	for i in range(len(ChainDictionary)):
		if chainList[i] in alignedChain:
			continue
		else:
			similarChainForChaini = [] #all the chains that are similar to chain i
			similarChainForChaini.append(chainList[i])
			alignedChain.append(chainList[i])
			for j in range(i+1,len(ChainDictionary)):
				aligner.mode = "local"
				score = aligner.score(ChainDictionary[chainList[i]],ChainDictionary[chainList[j]])
				lengthCom = min(len(ChainDictionary[chainList[i]]),len(ChainDictionary[chainList[j]]))/max(len(ChainDictionary[chainList[i]]),len(ChainDictionary[chainList[j]]))
				similarity = score/min(len(ChainDictionary[chainList[i]]),len(ChainDictionary[chainList[j]]))
				if similarity >= 0.95 and lengthCom >= 0.95:
					similarChainForChaini.append(chainList[j])
					alignedChain.append(chainList[j])
		if len(similarChainForChaini) > 1:
			similarChainList.append(similarChainForChaini)
	print(similarChainList)
	return similarChainList


#this function is to figure out the interaction inside a given protein
#input: a dictionary that contains the chain information of a pdbfile, a distance cutoff which decides whether the two atoms are connected, 
#       and the number of connected residues that decides whether the two chains are connected
#output: a list that contains the interacting chains[[chain i, chain j], [chain k, chain m],...]

def connect_chain(chains,cutOff,numberOfConnection):
    connectChain = []
    connectStrength = {} #keys: connect chain; values: calculated energy
    chainList = list(chains.keys())
    aaDic = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','LYS':'K','ARG':'R',
                        'HIS':'H','ASP':'D',"ASN":'N','GLU':'E','GLN':'Q','PRO':'P','TRP':'W','PHE':'F',
                        'TYR':'Y','MET':'M','CYS':'C'}
    energyDic = {'CYS-CYS':5.44, 'CYS-MET':4.99, 'CYS-PHE':5.80, 'CYS-ILE':5.50, 'CYS-LEU':5.83,'CYS-VAL':4.96,'CYS-TRP':4.95,'CYS-TYR':4.16,
                    'CYS-ALA':3.57,'CYS-GLY':3.16,'CYS-THR':3.11,'CYS-SER':2.86,'CYS-ASN':2.59,'CYS-GLN':2.85,'CYS-ASP':2.41,'CYS-GLU':2.27,'CYS-HIS':3.60,
                    'CYS-ARG':2.57,'CYS-LYS':1.95,'CYS-PRO':3.07, 'MET-MET':5.46,'MET-PHE':6.56,'MET-ILE':6.02,'MET-LEU':6.41,'MET-VAL':5.32,'MET-TRP':5.55,
                    'MET-TYR':4.91,'MET-ALA':3.94,'MET-GLY':3.39,'MET-THR':3.51,'MET-SER':3.03,'MET-ASN':2.95,'MET-GLN':3.30,'MET-ASP':2.57,'MET-GLU':2.89,
                    'MET-HIS':3.98,'MET-ARG':3.12,'MET-LYS':2.48,'MET-PRO':3.56,'PHE-PHE':7.26,'PHE-ILE':6.84,'PHE-LEU':7.28,'PHE-VAL':6.29,'PHE-TRP':6.16,
                    'PHE-TYR':5.66,'PHE-ALA':4.81,'PHE-GLY':4.13,'PHE-THR':4.28,'PHE-SER':4.02,'PHE-ASN':3.75,'PHE-GLN':4.10,'PHE-ASP':3.48,'PHE-GLU':3.56,
                    'PHE-HIS':4.77,'PHE-ARG':3.98,'PHE-LYS':3.36,'PHE-PRO':4.25,'ILE-ILE':6.54,'ILE-LEU':7.04,'ILE-VAL':6.05,'ILE-TRP':5.78,'ILE-TYR':5.25,
                    'ILE-ALA':4.58,'ILE-GLY':3.78,'ILE-THR':4.03,'ILE-SER':3.52,'ILE-ASN':3.24,'ILE-GLN':3.67,'ILE-ASP':3.17,'ILE-GLU':3.27,'ILE-HIS':4.14,
                    'ILE-ARG':3.63,'ILE-LYS':3.01,'ILE-PRO':3.76,'LEU-LEU':7.37,'LEU-VAL':6.48,'LEU-TRP':6.14,'LEU-TYR':5.67,'LEU-ALA':4.91,'LEU-GLY':4.16,
                    'LEU-THR':4.34,'LEU-SER':3.92,'LEU-ASN':3.74,'LEU-GLN':4.04,'LEU-ASP':3.40,'LEU-GLU':3.59,'LEU-HIS':4.54,'LEU-ARG':4.03,'LEU-LYS':3.37,
                    'LEU-PRO':4.20,'VAL-VAL':5.52,'VAL-TRP':5.18,'VAL-TYR':4.62,'VAL-ALA':4.04,'VAL-GLY':3.38,'VAL-THR':3.46,'VAL-SER':3.05,'VAL-ASN':2.83,
                    'VAL-GLN':3.07,'VAL-ASP':2.48,'VAL-GLU':2.67,'VAL-HIS':3.58,'VAL-ARG':3.07,'VAL-LYS':2.49,'VAL-PRO':3.32,'TRP-TRP':5.06,'TRP-TYR':4.66,
                    'TRP-ALA':3.82,'TRP-GLY':3.42,'TRP-THR':3.22,'TRP-SER':2.99,'TRP-ASN':3.07,'TRP-GLN':3.11,'TRP-ASP':2.84,'TRP-GLU':2.99,'TRP-HIS':3.98,
                    'TRP-ARG':3.41,'TRP-LYS':2.69,'TRP-PRO':3.73,'TYR-TYR':4.17,'TYR-ALA':3.36,'TYR-GLY':3.01,'TYR-THR':3.01,'TYR-SER':2.78,'TYR-ASN':2.76,
                    'TYR-GLN':2.97,'TYR-ASP':2.76,'TYR-GLU':2.79,'TYR-HIS':3.52,'TYR-ARG':3.16,'TYR-LYS':2.60,'TYR-PRO':3.19,'ALA-ALA':2.72,'ALA-GLY':2.31,
                    'ALA-THR':2.32,'ALA-SER':2.01,'ALA-ASN':1.84,'ALA-GLN':1.89,'ALA-ASP':1.70,'ALA-GLU':1.51,'ALA-HIS':2.41,'ALA-ARG':1.83,'ALA-LYS':1.31,
                    'ALA-PRO':2.03,'GLY-GLY':2.24,
                    'GLY-THR':2.08,'GLY-SER':1.82,'GLY-ASN':1.74,'GLY-GLN':1.66,'GLY-ASP':1.59,'GLY-GLU':1.22,'GLY-HIS':2.15,'GLY-ARG':1.72,'GLY-LYS':1.15,
                    'GLY-PRO':1.87,'THR-THR':2.12,'THR-SER':1.96,'THR-ASN':1.88,'THR-GLN':1.90,'THR-ASP':1.80,'THR-GLU':1.74,'THR-HIS':2.42,'THR-ARG':1.90,'THR-LYS':1.31,
                    'THR-PRO':1.90,'SER-SER':1.67,'SER-ASN':1.58,'SER-GLN':1.49,'SER-ASP':1.63,'SER-GLU':1.48,'SER-HIS':2.11,'SER-ARG':1.62,'SER-LYS':1.05,
                    'SER-PRO':1.57,'ASN-ASN':1.68,'ASN-GLN':1.71,'ASN-ASP':1.68,'ASN-GLU':1.51,'ASN-HIS':2.08,'ASN-ARG':1.64,'ASN-LYS':1.21,
                    'ASN-PRO':1.53,'GLN-GLN':1.54,'GLN-ASP':1.46,'GLN-GLU':1.42,'GLN-HIS':1.98,'GLN-ARG':1.80,'GLN-LYS':1.29,
                    'GLN-PRO':1.73,'ASP-ASP':1.21,'ASP-GLU':1.02,'ASP-HIS':2.32,'ASP-ARG':2.29,'ASP-LYS':1.68,
                    'ASP-PRO':1.33,'GLU-GLU':0.91,'GLU-HIS':2.15,'GLU-ARG':2.27,'GLU-LYS':1.80,
                    'GLU-PRO':1.26,'HIS-HIS':3.05,'HIS-ARG':2.16,'HIS-LYS':1.35,
                    'HIS-PRO':2.25,'ARG-ARG':1.55,'ARG-LYS':0.59,
                    'ARG-PRO':1.70,'LYS-LYS':0.12,
                    'LYS-PRO':0.97,'PRO-PRO':1.75, 'MET-CYS':4.99, 'PHE-CYS':5.80, 'ILE-CYS':5.50, 'LEU-CYS':5.83,'VAL-CYS':4.96,'TRP-CYS':4.95,'TYR-CYS':4.16,
                    'ALA-CYS':3.57,'GLY-CYS':3.16,'THR-CYS':3.11,'SER-CYS':2.86,'ASN-CYS':2.59,'GLN-CYS':2.85,'ASP-CYS':2.41,'GLU-CYS':2.27,'HIS-CYS':3.60,
                    'ARG-CYS':2.57,'LYS-CYS':1.95,'PRO-CYS':3.07,'PHE-MET':6.56,'ILE-MET':6.02,'LEU-MET':6.41,'VAL-MET':5.32,'TRP-MET':5.55,
                    'TYR-MET':4.91,'ALA-MET':3.94,'GLY-MET':3.39,'THR-MET':3.51,'SER-MET':3.03,'ASN-MET':2.95,'GLN-MET':3.30,'ASP-MET':2.57,'GLU-MET':2.89,
                    'HIS-MET':3.98,'ARG-MET':3.12,'LYS-MET':2.48,'PRO-MET':3.56,'ILE-PHE':6.84,'LEU-PHE':7.28,'VAL-PHE':6.29,'TRP-PHE':6.16,
                    'TYR-PHE':5.66,'ALA-PHE':4.81,'GLY-PHE':4.13,'THR-PHE':4.28,'SER-PHE':4.02,'ASN-PHE':3.75,'GLN-PHE':4.10,'ASP-PHE':3.48,'GLU-PHE':3.56,
                    'HIS-PHE':4.77,'ARG-PHE':3.98,'LYS-PHE':3.36,'PRO-PHE':4.25,'LEU-ILE':7.04,'VAL-ILE':6.05,'TRP-ILE':5.78,'TYR-ILE':5.25,
                    'ALA-ILE':4.58,'GLY-ILE':3.78,'THR-ILE':4.03,'SER-ILE':3.52,'ASN-ILE':3.24,'GLN-ILE':3.67,'ASP-ILE':3.17,'GLU-ILE':3.27,'HIS-ILE':4.14,
                    'ARG-ILE':3.63,'LYS-ILE':3.01,'PRO-ILE':3.76,'VAL-LEU':6.48,'TRP-LEU':6.14,'TYR-LEU':5.67,'ALA-LEU':4.91,'GLY-LEU':4.16,
                    'THR-LEU':4.34,'SER-LEU':3.92,'ASN-LEU':3.74,'GLN-LEU':4.04,'ASP-LEU':3.40,'GLU-LEU':3.59,'HIS-LEU':4.54,'ARG-LEU':4.03,'LYS-LEU':3.37,
                    'PRO-LEU':4.20,'TRP-VAL':5.18,'TYR-VAL':4.62,'ALA-VAL':4.04,'GLY-VAL':3.38,'THR-VAL':3.46,'SER-VAL':3.05,'ASN-VAL':2.83,
                    'GLN-VAL':3.07,'ASP-VAL':2.48,'GLU-VAL':2.67,'HIS-VAL':3.58,'ARG-VAL':3.07,'LYS-VAL':2.49,'PRO-VAL':3.32,'TYR-TRP':4.66,
                    'ALA-TRP':3.82,'GLY-TRP':3.42,'THR-TRP':3.22,'SER-TRP':2.99,'ASN-TRP':3.07,'GLN-TRP':3.11,'ASP-TRP':2.84,'GLU-TRP':2.99,'HIS-TRP':3.98,
                    'ARG-TRP':3.41,'LYS-TRP':2.69,'PRO-TRP':3.73,'ALA-TYR':3.36,'GLY-TYR':3.01,'THR-TYR':3.01,'SER-TYR':2.78,'ASN-TYR':2.76,
                    'GLN-TYR':2.97,'ASP-TYR':2.76,'GLU-TYR':2.79,'HIS-TYR':3.52,'ARG-TYR':3.16,'LYS-TYR':2.60,'PRO-TYR':3.19,'GLY-ALA':2.31,
                    'THR-ALA':2.32,'SER-ALA':2.01,'ASN-ALA':1.84,'GLN-ALA':1.89,'ASP-ALA':1.70,'GLU-ALA':1.51,'HIS-ALA':2.41,'ARG-ALA':1.83,'LYS-ALA':1.31,
                    'PRO-ALA':2.03,
                    'THR-GLY':2.08,'SER-GLY':1.82,'ASN-GLY':1.74,'GLN-GLY':1.66,'ASP-GLY':1.59,'GLU-GLY':1.22,'HIS-GLY':2.15,'ARG-GLY':1.72,'LYS-GLY':1.15,
                    'PRO-GLY':1.87,'SER-THR':1.96,'ASN-THR':1.88,'GLN-THR':1.90,'ASP-THR':1.80,'GLU-THR':1.74,'HIS-THR':2.42,'ARG-THR':1.90,'LYS-THR':1.31,
                    'PRO-THR':1.90,'ASN-SER':1.58,'GLN-SER':1.49,'ASP-SER':1.63,'GLU-SER':1.48,'HIS-SER':2.11,'ARG-SER':1.62,
                    'LYS-SER':1.05,
                    'PRO-SER':1.57,'GLN-ASN':1.71,'ASP-ASN':1.68,'GLU-ASN':1.51,'HIS-ASN':2.08,'ARG-ASN':1.64,'LYS-ASN':1.21,
                    'PRO-ASN':1.53,'ASP-GLN':1.46,'GLU-GLN':1.42,'HIS-GLN':1.98,'ARG-GLN':1.80,'LYS-GLN':1.29,
                    'PRO-GLN':1.73,'GLU-ASP':1.02,'HIS-ASP':2.32,'ARG-ASP':2.29,'LYS-ASP':1.68,
                    'PRO-ASP':1.33,'HIS-GLU':2.15,'ARG-GLU':2.27,'LYS-GLU':1.80,
                    'PRO-GLU':1.26,'ARG-HIS':2.16,'LYS-HIS':1.35,
                    'PRO-HIS':2.25,'LYS-ARG':0.59,
                    'PRO-ARG':1.70,
                    'PRO-LYS':0.97}

    ns = []
    for i in range(len(chains)):
        ns.append(NeighborSearch(list(chainList[i].get_atoms())))

    int_chain = []

    for i in range(len(chainList)):
        for j in range(i+1, len(chainList)):
            ns2 = ns[j]
            pairRes=[]
            ekey_ref=[]
            interacting_chains = []
            for residue1 in chainList[i]:
                for atom1 in residue1:
                    closest_atoms2 = ns2.search(atom1.get_coord(), 3.5)
                    for atom2 in closest_atoms2:
                        residue_name1 = residue1.get_resname()
                        residue_name2 = atom2.get_parent().get_resname()
                        unique_int = str(residue1.get_id()[1])+'-'+str(atom2.get_parent().get_id()[1])
                        if residue_name1 in aaDic and residue_name2 in aaDic and unique_int not in ekey_ref:
                            ekey = residue_name1+'-'+residue_name2
                            ekey_ref.append(unique_int)
                            pairRes.append(ekey)
                            interacting_chains.append((chainList[i], chainList[j]))


            if len(interacting_chains) >= cutOff:
                interacting_chains = interacting_chains[0]
                int_chain.append(interacting_chains)
                strength = 0
                for r in pairRes:
                    strength = strength + energyDic[r]
                connectStrength[tuple(interacting_chains)]=strength
    print(int_chain)
    print(connectStrength)
    return int_chain, connectStrength

#this function tells the kinds of connection in the "connectChain". A:B or A:A/B:B
#input: similarChainList, e.g [[chain A, chain B, chian C],[Chain D, Chain E]];connectChain[[chainA, chainB],[chainB, chain C]...]
#output: the number of total connection; the number of hetero A:B bonds; the number of homo A:A/B:B/C:C bonds
def bonds_genre(similarChainList, connectChain):
	bondsTypes = []
	totalConnect = len(connectChain)
	numberofHeteroBonds = 0
	numberofHomoBonds = 0
	homoChainPair = []
	heteChainPair = []
	for i, connection in enumerate(connectChain):
		judge = 0
		for j in range(len(similarChainList)):
			if connection[0] in similarChainList[j] and connection[1] in similarChainList[j]:
				judge += 1
			else:
				continue
		if judge == 0:
			numberofHeteroBonds += 1
			heteChainPair.append(connection)
		else:
			numberofHomoBonds += 1
			homoChainPair.append(connection)
	bondsTypes.append(totalConnect)
	bondsTypes.append(numberofHomoBonds)
	bondsTypes.append(numberofHeteroBonds)
	
	return bondsTypes,homoChainPair,heteChainPair

allHetero = 0
twoHomoTwoHetero = 0
twoHomoaTwohomob = 0
threeHomoOneHetero = 0
allHomo = 0
totalProtein = 0
Struc = []
AAAA_total = 0
AAAA_homo = 0
AAAA_hete = 0
AAAB_total = 0
AAAB_homo = 0
AAAB_hete = 0
AABB_total = 0
AABB_homo = 0
AABB_hete = 0
AABC_total = 0
AABC_homo = 0
AABC_hete = 0
ABCD_total = 0
ABCD_homo = 0
ABCD_hete = 0
circleStruc_ABCD = 0
chainStruc_ABCD =0
circleStruc_AABC = 0
chainStruc_AABC =0
circleStruc_AAAB = 0
chainStruc_AAAB =0
circleStruc_AABB = 0
chainStruc_AABB =0
circleStruc_AAAA = 0
chainStruc_AAAA =0
AAAA_homoStrengthDistri = []
AAAA_heteStrengthDistri = []
AAAB_homoStrengthDistri = []
AAAB_heteStrengthDistri = []
AABB_homoStrengthDistri = []
AABB_heteStrengthDistri = []
AABC_homoStrengthDistri = []
AABC_heteStrengthDistri = []
ABCD_homoStrengthDistri = []
ABCD_heteStrengthDistri = []

connectLessThanThreeFile = []
invalidFile = []
ConnectSummary = []
ConnectStrengthSummary = []
cutOff = float(sys.argv[1]) 
numberOfConnection = int(sys.argv[2])
rate1 = float(sys.argv[3])
# path = '/Users/cmdb/Dropbox/Mac/Desktop/pdb'
path = '/home/hsohail1/JohnsonLab/4chains'
files = os.listdir(path)
#print(files)
filenumber = len(files)
picknumber1 = int(filenumber*rate1)
sample1 = random.sample(files,picknumber1)
print (sample1)

for i, file in enumerate(sample1):
	print(i)
	try:
		print(file)
		validChainDic = chain_basic_info(file)[0]
		similarChainList = similar_chain(validChainDic)
		chainNumber = len(validChainDic)
		print(chainNumber)
		if chainNumber == 4:
			conn = connect_chain(validChainDic,cutOff,numberOfConnection)
			print(conn)
			connectChain = conn[0]
			if len(connectChain) > 2:
				StrengthChain = conn[1]			
				bondsTypes = bonds_genre(similarChainList,connectChain)
				homoChainPair = bondsTypes[1]
				heteChainPair = bondsTypes[2]
				print(bondsTypes)
				totalProtein += 1
				ConnectSummary.append(connectChain)
				ConnectStrengthSummary.append(StrengthChain)
				print(totalProtein)
				if len(similarChainList) == 0:
					allHetero += 1
					if len(connectChain) >=4:
						circleStruc_ABCD += 1
					else:
						chainStruc_ABCD += 1
					ABCD_total = ABCD_total + bondsTypes[0][0]
					ABCD_homo = ABCD_homo + bondsTypes[0][1]
					ABCD_hete = ABCD_hete + bondsTypes[0][2]
					for chainp in connectChain:
						if chainp in homoChainPair:
							ABCD_homoStrengthDistri.append(StrengthChain[tuple(chainp)])
						elif chainp in heteChainPair:
							ABCD_heteStrengthDistri.append(StrengthChain[tuple(chainp)])
						else:
							print("uncatolized pair")
							print(chainp)
				if len(similarChainList) == 1 and len(similarChainList[0]) == 2:
					twoHomoTwoHetero += 1
					if len(connectChain) >=4:
						circleStruc_AABC += 1
					else:
						chainStruc_AABC += 1
					AABC_total = AABC_total + bondsTypes[0][0]
					AABC_homo = AABC_homo + bondsTypes[0][1]
					AABC_hete = AABC_hete + bondsTypes[0][2]
					for chainp in connectChain:
						if chainp in homoChainPair:
							AABC_homoStrengthDistri.append(StrengthChain[tuple(chainp)])
						elif chainp in heteChainPair:
							AABC_heteStrengthDistri.append(StrengthChain[tuple(chainp)])
						else:
							print("uncatolized pair")
							print(chainp)
				if len(similarChainList) == 1 and len(similarChainList[0]) == 3:
					threeHomoOneHetero += 1
					if len(connectChain) >=4:
						circleStruc_AAAB += 1
					else:
						chainStruc_AAAB += 1
					AAAB_total = AAAB_total + bondsTypes[0][0]
					AAAB_homo = AAAB_homo + bondsTypes[0][1]
					AAAB_hete = AAAB_hete + bondsTypes[0][2]
				if len(similarChainList) == 1 and len(similarChainList[0]) == 4:
					allHomo += 1
					if len(connectChain) >=4:
						circleStruc_AAAA += 1
					else:
						chainStruc_AAAA += 1
					AAAA_total = AAAA_total + bondsTypes[0][0]
					AAAA_homo = AAAA_homo + bondsTypes[0][1]
					AAAA_hete = AAAA_hete + bondsTypes[0][2]
					for chainp in connectChain:
						if chainp in homoChainPair:
							AAAA_homoStrengthDistri.append(StrengthChain[tuple(chainp)])
						elif chainp in heteChainPair:
							AAAA_heteStrengthDistri.append(StrengthChain[tuple(chainp)])
						else:
							print("uncatolized pair")
							print(chainp)
				if len(similarChainList) == 2:
					twoHomoaTwohomob += 1
					if len(connectChain) >=4:
						circleStruc_AABB += 1
					else:
						chainStruc_AABB += 1
					AABB_total = AABB_total + bondsTypes[0][0]
					AABB_homo = AABB_homo + bondsTypes[0][1]
					AABB_hete = AABB_hete + bondsTypes[0][2]
					for chainp in connectChain:
						if chainp in homoChainPair:
							AABB_homoStrengthDistri.append(StrengthChain[tuple(chainp)])
						elif chainp in heteChainPair:
							AABB_heteStrengthDistri.append(StrengthChain[tuple(chainp)])
						else:
							print("uncatolized pair")
							print(chainp)
			else:
				connectLessThanThreeFile.append(file)
				continue

	except Exception as ex:
		print(ex)
		invalidFile.append(file)
	

print("The invalid files that have less than three bonds or less than four chains are: ")
print(connectLessThanThreeFile)
print("The unexplained files are:")
print(invalidFile)


AAAABondTypes = [AAAA_total,AAAA_homo,AAAA_hete]
AAABBondTypes = [AAAB_total,AAAB_homo,AAAB_hete]
AABBBondTypes = [AABB_total,AABB_homo,AABB_hete]
AABCBondTypes = [AABC_total,AABC_homo,AABC_hete]
ABCDBondTypes = [ABCD_total,ABCD_homo,ABCD_hete]
distributionList = []
distributionList.append(allHomo/totalProtein)
distributionList.append(twoHomoaTwohomob/totalProtein)
distributionList.append(twoHomoTwoHetero/totalProtein)
distributionList.append(threeHomoOneHetero/totalProtein)
distributionList.append(allHetero/totalProtein)
print("components-AAAA-AABB-AABC-AAAB-ABCD-total")
print(allHomo)
print(twoHomoaTwohomob)
print(twoHomoTwoHetero)
print(threeHomoOneHetero)
print(allHetero)
print(totalProtein)
print("Bond types: AAAA-AAAB-AABB-AABC-ABCD, total/homo/hetero")
print(AAAABondTypes)
print(AAABBondTypes)
print(AABBBondTypes)
print(AABCBondTypes)
print(ABCDBondTypes)
print(AAAA_homoStrengthDistri)
print(AAAA_heteStrengthDistri)
print(AAAB_homoStrengthDistri)
print(AAAB_heteStrengthDistri)
print(AABB_homoStrengthDistri)
print(AABB_heteStrengthDistri)
print(AABC_homoStrengthDistri)
print(AABC_heteStrengthDistri)
print(ABCD_homoStrengthDistri)
print(ABCD_heteStrengthDistri)

#plot the distrobution of repeated or unique subunits
x = [0, 1, 2, 3, 4]
fig, ax = plt.subplots()
ax.bar(x,distributionList)
ax.set_xticks(x,['AAAA','AABB','AABC','AAAB','ABCD'])
ax.set_ylabel("Frequency")
ax.set_xlabel("Subunits identity")
plt.show()

#plot the overall structure distribution of chain structure or circular structure
x_s = [0,1,2,3,4,5,6,7,8,9]
Struc.append(circleStruc_AAAA)
Struc.append(chainStruc_AAAA)
Struc.append(circleStruc_AABB)
Struc.append(chainStruc_AABB)
Struc.append(circleStruc_AABC)
Struc.append(chainStruc_AABC)
Struc.append(circleStruc_AAAB)
Struc.append(chainStruc_AAAB)
Struc.append(circleStruc_ABCD)
Struc.append(chainStruc_ABCD)
print("structure_AAAA-AABB-AABC-AAAB-ABCD/circle-chain")
print(Struc)
fig, ax = plt.subplots()
ax.bar(x_s,Struc)
ax.set_xticks(x_s,['Circle_AAAA','Chain_AAAA','Circle_AABB','Chain_AABB','Circle_AABC','Chain_AABC','Circle_AAAB','Chain_AAAB','Circle_ABCD',
	'Chain_ABCD'])
plt.show()


fig, (ax1,ax2) = plt.subplots(1,2) 
fig.suptitle('bonds strength for AAAA')
ax1.hist(AAAA_homoStrengthDistri, density = True)
ax1.set_title('homobonds')
ax2.hist(AAAA_heteStrengthDistri, density = True)
ax2.set_title('heterbonds')
plt.show()

fig, (ax1,ax2) = plt.subplots(1,2) 
fig.suptitle('bonds strength for AAAB')
ax1.hist(AAAB_homoStrengthDistri, density = True)
ax1.set_title('homobonds')
ax2.hist(AAAB_heteStrengthDistri, density = True)
ax2.set_title('heterbonds')
plt.show()

fig, (ax1,ax2) = plt.subplots(1,2) 
fig.suptitle('bonds strength for AABB')
ax1.hist(AABB_homoStrengthDistri, density = True)
ax1.set_title('homobonds')
ax2.hist(AABB_heteStrengthDistri, density = True)
ax2.set_title('heterbonds')
plt.show()

fig, (ax1,ax2) = plt.subplots(1,2) 
fig.suptitle('bonds strength for AABC')
ax1.hist(AABC_homoStrengthDistri, density = True)
ax1.set_title('homobonds')
ax2.hist(AABC_heteStrengthDistri, density = True)
ax2.set_title('heterbonds')
plt.show()

fig, (ax1,ax2) = plt.subplots(1,2) 
fig.suptitle('bonds strength for ABCD')
ax1.hist(ABCD_homoStrengthDistri, density = True)
ax1.set_title('homobonds')
ax2.hist(ABCD_heteStrengthDistri, density = True)
ax2.set_title('heterbonds')
plt.show()

# with  open("connection_summary.txt","w") as f:
# 	f.writelines(ConnectSummary)
# 	f.writelines('\n')
# 	f.writelines(ConnectStrengthSummary)

print("the connect summary is:")
print(ConnectSummary)
print("the strength summary is:")
print(ConnectStrengthSummary)