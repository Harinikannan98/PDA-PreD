#!/opt/websites/anaconda/envs/tf37/bin/python
'''#!/usr/bin/env python3.7'''
import cgi, os
from functools import total_ordering
import sys
import time
import shutil
import subprocess
import math
import numpy as np
from os import path
#import ssl
#ssl._create_default_https_context = ssl._create_unverified_context
import uuid
import cgitb;
cgitb.enable()
import timeit
start = timeit.default_timer()
#print ('Content-Type: text/html\r\n')
#print('/r/n')
import sys
sys.stdout.flush()

#os.chdir("vossvolvox-master/")
#print(os.listdir())


#p = subprocess.Popen("./bin -i 1aay.pdb >ooo",shell=True,stdout=subprocess.PIPE)
#p = subprocess.Popen("./bin.exe ", "-i "," apo_1ais.pdb","> ooo",shell=True,stdout=subprocess.PIPE)
#time.sleep(5)
#h=open("ooo").readlines()
#print(h)
#print(i)
#form = cgi.FieldStorage()
do=input("Enter the type of input: 'pdb-id' or 'pdb-file'")
if do=="pdb-id":
	inn=input("Enter the PDB-ID and protein/DNA/chains (optional) as 1AAY-A-B or 1AAY--B or 1AAY-A")
	inn=inn.split("-")
	if len(inn)>1:
		pdb_id=inn[0]
		chain=inn[1]
	if (len(inn)>2):
		dchain=inn[2]
	else:
		pdb_id=inn[0]
		chain=""
		dchain=""
	pdb_id=pdb_id.lower()
if do=='pdb-file':
	inn=input("Enter the PDB file name in current folder")
	pdb_file=inn
	pdb_id=""
	chain=""
	dchain=""
dna_strand=input("Please select if the DNA is single stranded or double stranded (type ds or ss)")
#dna_strand=form.getvalue('dna_strand')


#os.system("rm -r pd_res_*")
if __name__ == '__main__':
	
	#print ('Content-Type: text/html\r\n')
	#print('/r/n')
	#os.system("chmod +x pdpredict.py")
	'''
	redirectURL = "pdpredict.py"
	print ('Content-Type: text/html\r\n',flush=True)
	print ('<html>',flush=True)
	print ('  <head>',flush=True)
	print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL,flush=True)
	'''
	#print ('	<title>You are going to be redirected</title>')
	#print ("Use this url to obtain data: <a href='{}'>Result_page</a>".format(redirectURL)
	#print ('Content-Type: text/html\r\n')
	
	#print(pdb_id+'/r/n')
	#print ('<html><head><title>Input error</title><body>PDB ID is missing</body></html>')
	#exit()
	'''g=open("index.txt").readlines()
	for gg in g:
		print(gg,flush=True)
		#gg=gg.rstrip()
		#print ("""{}""").format(gg)
		#print("""print (\"\"\"{}\"\"\")""".format(gg))
	print("calculating features....",flush=True)
	#print(fileitem,flush=True)
	#print(pdb_id,flush=True)
	h=open("footer.txt").readlines()
	for hh in h:
		#hh=hh.rstrip()
		print (hh,flush=True)
		#print("""print (\"\"\"{}\"\"\")""".format(hh))
		#print("\n")
	'''
	#os.system("find -mmin +20 -type d -exec rmdir !('tmp') > files_older_last_deleted 2>&1")
	#find dir_* -mmin +34 -type d -exec rm -r {} \;
	#if chain=="":
	#	print("chain,dchain")
	#pdb_id='1aay'
	#chain='A'
	#dchain='B'
	#fileitem=''
	#dna_strand='ss'
	#print(chain,dchain,pdb_id,fileitem,dna_strand)
	#print ('Content-Type: text/html\r\n')
	#print('/r/n')
	if pdb_id!="" and chain!="" and dchain!="" and do=='pdb-id' and dna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=chain
		dchain=dchain
		method=1
	elif pdb_id!="" and chain!="" and dchain=="" and do=='pdb-id' and dna_strand!="": 
		pdb_id_up=pdb_id.upper()
		chain=chain
		dchain=""
		method=1
	elif pdb_id!="" and chain=="" and dchain!="" and do=='pdb-id' and dna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=""
		dchain=dchain
		method=1
	elif pdb_id!="" and chain=="" and dchain=="" and do=='pdb-id' and dna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=""
		dchain=""
		print("in")
		method=1
	elif pdb_id==""and do=='pdb-id':
		print ('Content-Type: text/html\r\n')
		#print(fileitem,flush=True)
		print('/r/n')
		print ('<html><head><title>Input error</title><body>PDB ID is missing</body></html>')
		exit()
	elif dna_strand=="":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Input error</title><body>Select DNA strand</body></html>')
		exit()
	elif pdb_id=="" and do=='pdb-file':
		method=2
		#print(fileitem,flush=True)
		#print (fileitem.filename)
		#if form['pdbf'].filename:
		#	fn = os.path.basename(form['pdbf'].filename)
		#	open('tmp/' + fn, 'wb').write(form['pdbf'].file.read())
	elif pdb_id=="" and chain!="" and dchain=="" and do=='pdb-file':
		#fileitem = form['pdbf']
		chain=chain
		#lchain=form.getvalue('chain')
		method=2
		#print (fileitem.filename)
		#if form['pdbf'].filename:
		#	fn = os.path.basename(form['pdbf'].filename)
		#	open('tmp/' + fn, 'wb').write(form['pdbf'].file.read())
	else:
		#print ('Content-Type: text/html\r\n')
		#print('/r/n')
		print ('<html><head><title>Input error</title><body>No input found</body></html>')
		exit()
	if dna_strand=='ds':
		model=input("Please enter the structural classification of protein (all-alpha/all-beta/alpha-beta/other)")
		if model=="all-alpha":
			model="alpha"
		if model=="beta":
			model="beta"
		if model=="alpha-beta":
			model="alphabeta"
		if model=="other" or model=="others":
			model="other"
		model1=input("Please enter the Functional classification of protein (Regulatory/other")
		if model1=="Regulatory" or model1=="regulatory":
			model1="reg"
		if model1=="other" or model1=="others":
			model1="nreg"
		#print(model) #model1='reg'
		#print(model1)
	if dna_strand=='ss':
		model1=""
		model=""
	#model=""
	#model1=""
	model2=dna_strand
	#model=""	#model2='ss' 
	#model1=""
	if model=="None" and model2!="ss":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Please select the structural classfication </title><body>Please select the structural classfication of the protein</body></html>')
		exit()
	if model1=="None" and model2!="ss":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Please select the Functional classfication </title><body>Please select the Functional classfication of the protein </body></html>')
		exit()
	if model2=="None":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Please select the  dna strand</title><body>please select type of DNA</body></html>')
		exit()
	
	if model=='alpha':
		modeltype="All Alpha protein"
	elif model=='beta':
		modeltype="All Beta protein"
	elif model=='alphabeta':
		modeltype="Alpha Beta protein"
		#print("i")
	elif model=='other':
		modeltype="other strucutral protein"
		modeltype1=""
	if model1=='reg':
		modeltype1="Regulatory"
	if model1=='nreg':
		modeltype1="Not Regulatory"
	elif model=='nreg':
		modeltype1="Not regulatory"
	if model2=='ds':
		modeltype2="Double strand"
	elif model2=='ss':
		modeltype2="Single strand"
		modeltype1=''
		modeltype=''
	#print(model+model2+model1)
#########################################################################################################	Single strand	 ###########################################################################################################
#########################################################################################################	Single strand	  ###########################################################################################################
#########################################################################################################	Single strand	  ###########################################################################################################
#########################################################################################################	Single strand	  ###########################################################################################################
	import shutil
	import wget
	import glob
	import Bio.PDB as bpdb
	from Bio.PDB import is_aa
	from Bio.PDB import PDBParser, PDBIO, Select
	import urllib
	import os
	import numpy as np
	import re
	import pandas as pd
	import math
	from Bio import PDB
	import warnings
	from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
	pdb_id_up=pdb_id
	print(pdb_id_up)
	print(method)
	if method==1:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		randname="pd_res_"+uuid.uuid4().hex
		path = os.path.join(dir_path, randname)
		os.mkdir(path,0o777)
		#os.system("chmod -R 777 {}".format(path))
		#os.system("chmod 777 bin")
		#os.system("chmod +x bin")
		os.chdir(path)
		#print(os.listdir())
		#print(path)
		#shutil.copyfile("../1aay.pdb", "1aay.pdb")
		#os.system("wget 'https://files.rcsb.org/download/{}.pdb'".format(pdb_id_up))
		#os.system("wget 'https://files.rcsb.org/download/{}.pdb1'".format(pdb_id_up))
		#print(pdb_id_up)
		filename=wget.download('https://files.rcsb.org/download/{}.pdb1'.format(pdb_id_up))
		filename1=wget.download('https://files.rcsb.org/download/{}.pdb'.format(pdb_id_up))
	elif method==2:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		randname=uuid.uuid4().hex
		path = os.path.join(dir_path, randname)
		os.mkdir(path,0o777)
		os.chdir(path)
		os.system("scp ../"+pdb_file+" input.pdb")
		pdb_id_up='input'
	#shutil.copyfile("../clean_pdb.py", "clean_pdb.py")
	#os.system(r"python3 clean_pdb.py {}.pdb".format(pdb_id_up))
	'''Process PDB'''
	#redirectURL = "%s/result.py" % randname
	#print ('Content-Type: text/html\r\n')
	#print ('<html>')
	#print ('  <head>')
	#print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
	#print('	<title>You are going to be redirected</title>')
	class ProtSelect(Select):
		warnings.simplefilter('ignore', PDBConstructionWarning)
		warnings.simplefilter('ignore', FutureWarning)
		def accept_residue(self, residue):
			if not is_aa(residue, standard=True):
				res = residue.id[0]
				if not res == "W":
					return True 
			else:
				return False
	class ProtSelect1(Select):
		def accept_residue(self, residue):
			warnings.simplefilter('ignore', PDBConstructionWarning)
			warnings.simplefilter('ignore', FutureWarning)
			if is_aa(residue, standard=True):
				return True 
			else:
				return False
	
	parser = PDBParser()
	structure = parser.get_structure(pdb_id_up, pdb_id_up+".pdb")
	modelll = structure[0]
	io = bpdb.PDBIO()
	io.set_structure(modelll)
	io.save('dna_'+pdb_id_up+'.pdb', ProtSelect())
	io.save('apo_'+pdb_id_up+'.pdb', ProtSelect1())
	if dchain!="":
		dchains=dchain.split(",")
		class ChainSelect(Select):
			def __init__(self, schain):
				self.schain = schain
			def accept_schain(self, schain):
				if schain.get_id() in self.schain:
					return 1
				else:
					return 0
		p = PDBParser(PERMISSIVE=1)	   
		structure = p.get_structure('dna'+pdb_id_up+'.pdb', 'dna_'+pdb_id_up+'.pdb')
		io_w_no_h = PDBIO()
		io_w_no_h.set_structure(structure)
		io_w_no_h.save('dna_'+pdb_id_up+'.pdb', ChainSelect(dchains))	
	if chain!="":
		chains=chain.split(",")
		class ChainSelect(Select):
			def __init__(self, schain):
				self.schain = schain
			def accept_schain(self, schain):
				if schain.get_id() in self.schain:
					return 1
				else:
					return 0
		p = PDBParser(PERMISSIVE=1)	   
		structure = p.get_structure('apo_'+pdb_id_up+'.pdb', 'apo_'+pdb_id_up+'.pdb')
		io_w_no_h = PDBIO()
		io_w_no_h.set_structure(structure)
		io_w_no_h.save('apo_'+pdb_id_up+'.pdb', ChainSelect(chains))		
	if dchain!="" or chain!="":	
		chains=chain.split(",")
		dchains=dchain.split(",")
		allchain=chains+dchains
		print(allchain)
		class ChainSelect(Select):
			def __init__(self, schain):
				self.schain = schain
			def accept_schain(self, schain):
				if schain.get_id() in self.schain:
					return 1
				else:
					return 0
		p = PDBParser(PERMISSIVE=1)
		structure = p.get_structure(pdb_id_up+'.pdb', pdb_id_up+'.pdb')
		io_w_no_h = PDBIO()
		io_w_no_h.set_structure(structure)
		io_w_no_h.save(pdb_id_up+'.pdb', ChainSelect(allchain))
	
	#sys.stdout.write("Content-type: text/plain\r\n\r\n")
	
	#shutil.copyfile("apo_"+pdb_id_up+".pdb", "../vossvolvox-master/apo_"+pdb_id_up+".pdb")
	#os.chdir("../vossvolvox-master/")
	#os.system("chmod 777 bin")
	#os.system("chmod +x bin")
	#os.system("chmod 777 apo_"+pdb_id_up+".pdb")
	#p = subprocess.Popen(["./bin ", "-i "," apo_",pdb_id_up,".pdb","> ooo"],shell=True,stdout=subprocess.PIPE)
	#print (p.stdout.read(), flush=True)
	#p.wait()
	#os.system("scp "+path+"/apo_"+pdb_id_up+".pdb .")
	#p4=subprocess.Popen("./vossvolvox-master/bin -i apo_"+pdb_id_up+".pdb> ooo", shell=True)
	#p4.wait()
	
	#os.chdir("../")
	#os.system("rm ooo.py")
	#p=subprocess.Popen("./vol_ex.sh 1hvn",shell=True)
		
	
	print("redirected")
	shutil.copyfile("../foldx", "foldx")
	shutil.copyfile("../rotabase.txt", "rotabase.txt")
	shutil.copyfile("../naccess", "naccess")
	shutil.copyfile("../bin","bin")
	shutil.copyfile("../style4.css", "style4.css")
	shutil.copyfile("../clean_pdb.py", "clean_pdb.py")
	#os.system(r"python3 clean_pdb.py {}.pdb".format(pdb_id_up))
	shutil.copyfile("../index.txt", "index.txt")
	shutil.copyfile("../footer.txt", "footer.txt")
	shutil.copyfile("../dssp","dssp")
	shutil.copyfile("../dna_4vdrch.csv", "dna_4vdrch.csv")
	shutil.copyfile("../aa_20vdrch.csv", "aa_20vdrch.csv")
	shutil.copyfile("../potential_res.csv","potential_res.csv")
	os.system("chmod -R +x {}".format("foldx"))
	os.system("chmod -R +x {}".format("naccess"))
	os.system("chmod -777 {}".format("naccess"))
	os.system("chmod -R +x {}".format("bin"))
	os.system("chmod -R +x {}".format("dssp"))
	os.system("chmod -R 777 {}".format(path))
	print ('Content-Type: text/html\r\n')
		#print(fileitem,flush=True)
	print('/r/n')
	'''
	shutil.copyfile("apo_"+pdb_id_up+".pdb","../vossvolvox-master/apo_"+pdb_id_up+".pdb")
	os.system("chmod -R 777 ../vossvolvox-master/apo_"+pdb_id_up+".pdb")
	os.chdir("../vossvolvox-master/src/")
	p3=subprocess.Popen("make vol>simpl",shell=True)
	p3.wait()
	os.chdir("../")
	os.system("chmod -R +x bin")
	os.system("chmod -R 777 bin")
	p4=subprocess.Popen("./bin -i apo_1noy.pdb > volume2",shell=True)
	p4.wait()
	os.system("chmod -R 777 volumeplease")
	shutil.copyfile("volumeplease","../volumeplease")
	os.chdir(path)
	shutil.copyfile("../volumeplease","volumeplease")
	#shutil.copyfile("../bin","bin.exe")
	#os.system("chmod -R +x {}".format("bin.exe"))
	
	os.system("./bin -i apo_1gwv.pdb ")
	p4.wait()
	with open("volppp") as fil:
		for f in fil.readlines():
			print(f)
			har=0
	print(har)
	'''
	'''
	redirectURL = "%s/result.py" % randname
	print ('Content-Type: text/html\r\n')
	print ('<html>')
	print ('  <head>')
	print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
	print('	<title>You are going to be redirected</title>')
	'''
	if model2=='ss':	
		with open("result.txt","w") as resultout:
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[0]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['DA','DG','DC','DT']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									if (distance<= 10):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('dna_4vdrch.csv')
			elec=[]
			atom_count=[]	
			
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					time.sleep(5)
					print(" ")
					sys.stdout.flush()
					print(" ")
					df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
					if len(df_inter.columns)>3:
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'interaction_energy.csv')
						print(" ")
						atoms_involve = inst2.f8_interaction_type(df_inter)
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						time.sleep(5)
						print(" ")
						sys.stdout.flush()
						print(" ")
						elec.append(energy_df['ele_energy'].sum())
						atom_count.append(df_at['count'].sum())
						print(" ")
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						#energy_dict = inst2.f9_energy_div(energy_df)
						#df_en=pd.Series(energy_dict).to_frame()
						#df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
			
			'''
			for f in glob.glob(pdb_id_up+'*interaction_energy.csv'):
				df=pd.read_csv(f)
				print(df.columns)
				
				print(elec)
			for f in glob.glob(pdb_id_up+'*_atom_count.csv'):
				df=pd.read_csv(f)
				print(df.columns)
				
				print(atom_count)
			'''

			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb="+pdb_id_up+".pdb --complexWithDNA=true", stdout=subprocess.PIPE, shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			fold_dict={}
			fold_dict={}
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				print(f)
				with open(f) as file:
					lis=file.readlines()
					print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				print(fold_dict)
			fold_energy_ion=fold_dict['energy Ionisation']
			fold_vdw=fold_dict['Van der Waals']
			print(fold_energy_ion)
			print(fold_vdw)
			print(elec)
			predval=-5.86507871e+00*float(fold_energy_ion)+3.11970945e+00*float(sum(elec))+8.84807513e-04*float(sum(atom_count))+4.96708922e-01*float(fold_vdw)-8.738514288740806
			predval="%.2f" % predval
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval)+" ± 0.24")
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
##########################################################################################					Double-all-alpha-reg		  ##########################################################################################
##########################################################################################					Double-all-alpha-reg			   ####alpha######################################################################################
##########################################################################################					Double-all-alpha-reg			   ##########################################################################################
##########################################################################################					Double-all-alpha-reg			   ##########################################################################################
	#print(model2,model,model1)
	if model2=='ds' and model=='alpha' and model1 =='reg':
		with open("result.txt","w") as resultout:
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[0]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['DA','DG','DC','DT']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									if (distance<= 10):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('dna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			cp_count=[]
			tot=[]
			mcmc=[]
			bind=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					if len(df_inter.columns)>3:
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'interaction_energy.csv')
						atoms_involve = inst2.f8_interaction_type(df_inter)
						time.sleep(5)
						print(" ")
						sys.stdout.flush()
						print(" ")
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						energy_dict = inst2.f9_energy_div(energy_df)
						df_en=pd.Series(energy_dict).to_frame()
						#df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
						#print(df_en)
						for p0 in df_inter['prot_resno'].tolist():
							res_bind_count.append(p0)
							bind.append(i+"_"+str(int(p0)))
						tot.append(df_at['count'].sum())
						if 'CP' in df_at.index:
							cp_count.append(df_at.loc[['CP']].values)
						mcmc.append(df_en.loc[['mc_mc']].values)
					
			print(res_bind_count)
			bind=list(set(bind))
			res_bind_uni=list(set(res_bind_count))
			time.sleep(5)
			print(" ")
			sys.stdout.flush()
			print(" ")
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb="+pdb_id_up+".pdb --complexWithDNA=true", shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				print(f)
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				print(fold_dict)
			total_residue=fold_dict['Number of Residues']
			bind_per=len(res_bind_uni)/int(total_residue)
			print(bind_per)
#******************************************************************************bind less than 0.29*****************************************************************	
			if bind_per>0.29:
				shift,tilt=[],[]
				method_bind=1
				helix_dipole=fold_dict['helix dipole']
				k=0
				c=0
				with open(r"apo_bind_"+pdb_id_up+".pdb","w") as file1:
					with open(r"apo_"+pdb_id_up+".pdb") as file:
						for i in file.readlines():
							#print(i[0:26])
							if i[0:4]=="ATOM":
								chnu=i[21:26].split()
								cha=chnu[0]
								number=chnu[1]
								#print(chain)
								#print(number)
								print(cha+"_"+number)
								print(bind[0:3])
								if cha+"_"+number in bind:
									file1.write(i)
				'''
				shutil.copyfile("apo_"+pdb_id_up+".pdb", "../apo_"+pdb_id_up+".pdb")
				os.chdir("../")
				os.system("chmod 777 bin")
				os.system("chmod +x bin")
				p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".pdb> ooo",shell=True)
				p4.wait()
				#p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".pdb> ooo",shell=True)
				#p4.wait()
				os.chdir(path)
				shutil.copyfile("../ooo", "ooo")
				'''
				shutil.copyfile("../3vvv.py","3vvv.py")
				shutil.copyfile("../vol_out.txt","vol_out.txt")
				p=subprocess.Popen("/opt/websites/anaconda/bin/python3.7 3vvv.py apo_bind"+pdb_id_up,shell=True)
				p.wait()
				with open("vol_out.txt") as file:
					for vol in file.readlines():
						vol1=vol[0]					
						surface1=vol[1]
				#shutil.copyfile("apo_"+pdb_id_up+".pdb", "../apo_"+pdb_id_up+".pdb")
				from Bio.PDB import PDBParser
				from Bio.PDB.DSSP import DSSP
				p = PDBParser()
				structure = p.get_structure("apo_"+pdb_id_up+".pdb", "apo_"+pdb_id_up+".pdb")
				modelp = structure[0]
				dssp = DSSP(modelp, "apo_"+pdb_id_up+".pdb")
				#print(dssp)
				from Bio.PDB.DSSP import dssp_dict_from_pdb_file
				dssp_tuple = dssp_dict_from_pdb_file("apo_"+pdb_id_up+".pdb")
				dssp_dict = dssp_tuple[0]
				print(bind)
				dssp_sec=[]
				data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
				for ds in bind:
					dssp_sec.append(data[dssp_dict[(ds.split("_")[0], (' ',int(ds.split("_")[1]), ' '))][1]])
				print(dssp_sec)					
				#dssp=dssp_sec
				#i[11:13].strip()+"_"+kk[1]
				#a_key = list(dssp.keys())
				#print(dssp[a_key])
				
				'''				
				os.chdir("../")
				p3=subprocess.Popen("./dssp -i apo_"+pdb_id_up+".pdb -o "+pdb_id_up+".dssp",shell=True)
				p3.wait()
				os.chdir(path)		
				shutil.copyfile("../"+pdb_id_up+".dssp", pdb_id_up+".dssp")		
				k=1
				
				data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
				dssp=[]
				with open(pdb_id_up+".dssp") as file:
					for i in file.readlines():
							if "#" in i:
									k=0
							if "!" in i:
									chain.append("-")
									ds.append("-")
							if k==0:
									if not "#" in i and not "!" in i:
											j=i.split(" ")
											#print( i[13:15].strip()+":"+i[0:11].strip())
											kk=i[0:11].strip()
											kk=[str for str in kk.split() if str.strip()]
											#print(kk)
											#print(f)
											#print(i[11:13].strip())
											if i[11:13].strip()+"_"+kk[1] in bind:
													dssp.append(data[i[16:18].strip()])
													chain=i[11:13].strip()+":"+kk[1]
													ds=data[i[16:18].strip()]
				'''
				dssp=dssp_sec
				tot=len(dssp)
				he=dssp.count("helix")
				be=dssp.count("Beta")
				co=dssp.count("coil")
				percent_coil=100*(co/tot)
				import urllib.request
				urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.pars",pdb_id_up+"_bp_step.pars")
				l=1
				k=0
				with open(pdb_id_up+"_bp_step.pars") as file:
					
					for j in file.readlines():
						if "***local step parameters***" in j:
							k=1
						if k==1:
							if not "#" in j:
								j=j.strip()
								#print(j)
								s=[str for str in j.split(" ") if str.strip()]
								shift.append(float(s[1]))
								tilt.append(float(s[4]))
								l=0
					if l==1:
						print(i+"error:single_base_step")
					shift1=sum(shift)/len(shift)
					tilt1=sum(tilt)/len(tilt)
				predval=-3.17737626e+00*float(helix_dipole)+-6.75059750e-04*float(vol1)+1.47194228e-01*float(tilt1)+1.43387174e-03*float(surface1)+-9.81150741e-02*float(percent_coil)+4.09567809e+00*float(shift1)-8.063297881608307
				predval="%.2f" % predval
				if pdb_id_up=='input':
					resultout.write("User input")
				else:
					resultout.write(pdb_id_up)
				resultout.write("\n")
				#print (chain)
				prot_chain=",".join(prot_chain)
				resultout.write(str(prot_chain))
				resultout.write("\n")
				#print (binlig5)
				binlig5=",".join(dna_chain)
				resultout.write(str(binlig5))
				resultout.write("\n")
				resultout.write(str(predval)+" ± 1.22")
				disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
				resultout.write("\n")
				resultout.write(str(disass))
			#****volume_bs and surface_bs*********
#******************************************************************************bind greater than 0.29*****************************************************************	
			if bind_per<0.29:
				method_bind=2
				fold_tor=fold_dict['torsional clash']
				cp1=sum(cp_count)/sum(tot)
				mcmc_energy=sum(mcmc)
				shutil.copyfile("../3vvv.py","3vvv.py")
				shutil.copyfile("../vol_out.txt","vol_out.txt")
				p=subprocess.Popen("/opt/websites/anaconda/bin/python3.7 3vvv.py apo_"+pdb_id_up,shell=True)
				p.wait()
				with open("vol_out.txt") as file:
					for vol in file.readlines():
						vol1=vol[0]					
						surface1=vol[1]	
				print(vol1)				
				double_shift,double_tilt=[],[]
				urllib.request.urlretrieve(r"http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".outs",pdb_id_up+"_summary.txt")
				k=0
				l=1
				with open(pdb_id_up+"_summary.txt") as file:
					for j in file.readlines():
						if "Local base-pair step parameters" in j or "Local base step parameters" in j:
							k=1
						if k==1:
							if "ave" in j or "a/a" in j:
								j=j.strip()
								s=[str for str in j.split(" ") if str.strip()]
								double_shift=s[1]
								double_tilt=s[4]
								k=0
								l=0
					if l==1:
						print("error-3DNA-double")	
				predval=2.26314218e+02*float(cp1)+2.31840985e-05*float(vol1)+-2.73724094e-01*float(double_tilt)+6.76226587e-02*float(mcmc_energy)+2.93204980e-01*float(fold_tor)+2.17691420e+00*float(double_shift)-8.063297881608307
				predval="%.2f" % predval
				if pdb_id_up=='input':
					resultout.write("User input")
				else:
					resultout.write(pdb_id_up)
				resultout.write("\n")
				#print (chain)
				prot_chain=",".join(prot_chain)
				resultout.write(str(prot_chain))
				resultout.write("\n")
				#print (binlig5)
				binlig5=",".join(dna_chain)
				resultout.write(str(binlig5))
				resultout.write("\n")
				resultout.write(str(predval)+" ± 0.94")
				disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
				resultout.write("\n")
				resultout.write(str(disass))
##########################################################################################					Double-all-alpha-nreg		  ##########################################################################################
##########################################################################################					Double-all-alpha-nreg		  ##########################################################################################
##########################################################################################					Double-all-alpha-nreg		  ##########################################################################################
##########################################################################################					Double-all-alpha-nreg		  ##########################################################################################
	if model2=='ds' and model=='alpha' and model1 =='nreg':
		with open("result.txt","w") as resultout:
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[0]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['DA','DG','DC','DT']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									if (distance<= 10):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('dna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			cp_count=[]
			mcmc=[]
			charged,nonpolar,polar=[],[],[]
			bind=[]
			for i in prot_chain:
				for j in dna_chain:
					####p5=subprocess.Popen(df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j))
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					time.sleep(5)
					print(" ")
					sys.stdout.flush()
					print(" ")
					df_inter = pd.DataFrame(df_inter)
					df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
					aa={'non-polar':['GLY','ALA','VAL','LEU','MET','ILE','PHE','TYR','TRP'],'polar':['SER','THR','CYS','PRO','ASN','GLN'],'charged':['LYS','ARG','HIS','ASP','GLU']}
					df_c=df_inter[['prot_res','prot_resno']].drop_duplicates(ignore_index=True)
					df_c.drop_duplicates(ignore_index=True)
					df_c=df_c.replace(aa['non-polar'],'non-polar')
					df_c=df_c.replace(aa['polar'],'polar')
					df_c=df_c.replace(aa['charged'],'charged')
					print(df_c)
					df_new=df_c.groupby(df_c.iloc[:,0])["prot_res"].count()
					print(df_new)
					x=df_new.iloc[0:3].to_list()
					charged.append(x[0])
					nonpolar.append(x[1])
					polar.append(x[2])
					
					for oo in list(set(df_inter['prot_resno'])):
						bind.append(i+"_"+str(int(oo)))
			from Bio.PDB import PDBParser
			from Bio.PDB.DSSP import DSSP
			p = PDBParser()
			structure = p.get_structure("apo_"+pdb_id_up+".pdb", "apo_"+pdb_id_up+".pdb")
			modelp = structure[0]
			dssp = DSSP(modelp, "apo_"+pdb_id_up+".pdb")
			#print(dssp)
			from Bio.PDB.DSSP import dssp_dict_from_pdb_file
			dssp_tuple = dssp_dict_from_pdb_file("apo_"+pdb_id_up+".pdb")
			dssp_dict = dssp_tuple[0]
			print(bind)
			dssp_sec=[]
			data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
			for ds in bind:
				dssp_sec.append(data[dssp_dict[(ds.split("_")[0], (' ',int(ds.split("_")[1]), ' '))][1]])
			print(dssp_sec)					
			#dssp=dssp_sec
			#i[11:13].strip()+"_"+kk[1]
			#a_key = list(dssp.keys())
			#print(dssp[a_key])
			
			'''				
			os.chdir("../")
			p3=subprocess.Popen("./dssp -i apo_"+pdb_id_up+".pdb -o "+pdb_id_up+".dssp",shell=True)
			p3.wait()
			os.chdir(path)		
			shutil.copyfile("../"+pdb_id_up+".dssp", pdb_id_up+".dssp")		
			k=1
			
			data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
			dssp=[]
			with open(pdb_id_up+".dssp") as file:
				for i in file.readlines():
						if "#" in i:
								k=0
						if "!" in i:
								chain.append("-")
								ds.append("-")
						if k==0:
								if not "#" in i and not "!" in i:
										j=i.split(" ")
										#print( i[13:15].strip()+":"+i[0:11].strip())
										kk=i[0:11].strip()
										kk=[str for str in kk.split() if str.strip()]
										#print(kk)
										#print(f)
										#print(i[11:13].strip())
										if i[11:13].strip()+"_"+kk[1] in bind:
												dssp.append(data[i[16:18].strip()])
												chain=i[11:13].strip()+":"+kk[1]
												ds=data[i[16:18].strip()]
			'''
			dssp=dssp_sec
			tot=len(dssp)
			he=dssp.count("helix")
			be=dssp.count("Beta")
			co=dssp.count("coil")
			percent_beta=100*(be/tot)
			print(percent_beta)
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb="+pdb_id_up+".pdb --complexWithDNA=true", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				print(f)
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				print(fold_dict)
			tilt=[]
			k=0
			c=0
			import urllib.request
			urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.pars",pdb_id_up+"_bp_step.pars")
			l=1
			k=0
			with open(pdb_id_up+"_bp_step.pars") as file:
				for j in file.readlines():
					print(j)
					if "***local step parameters***" in j:
						k=1
					if k==1:
						if not "#" in j:
							j=j.strip()
							#print(j)
							s=[str for str in j.split(" ") if str.strip()]
							tilt.append(float(s[4]))
							l=0
				if l==1:
					print(i+"error:single_base_step")
			tilt1=sum(tilt)/len(tilt)
			torsional_clash=fold_dict['torsional clash']
			vander_clash=fold_dict['Van der Waals clashes']
			#***charged***
			charged=sum(charged)

			#***percentage beta***
			print(torsional_clash, percent_beta,tilt1, vander_clash, charged)
			predval=1.52177601*float(torsional_clash)+0.11620005*float(percent_beta)+-0.11936982*float(tilt1)+-0.1088115*float(vander_clash)+0.00958888*float(charged)-11.232446301368846
			predval="%.2f" % predval
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval)+" ± 1.11")
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
			print(disass)
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
	if model2=='ds' and model=='beta' and model1 =='reg':
		with open("result.txt","w") as resultout:
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[0]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['DA','DG','DC','DT']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									if (distance<= 10):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			class ProtSelect(Select):
				warnings.simplefilter('ignore', PDBConstructionWarning)
				warnings.simplefilter('ignore', FutureWarning)
				def accept_residue(self, residue):
					if not is_aa(residue, standard=True):
						res = residue.id[0]
						if not res == "W":
							return True 
					else:
						return False
			class ProtSelect1(Select):
				def accept_residue(self, residue):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					if is_aa(residue, standard=True):
						return True 
					else:
						return False
			
			parser = PDBParser()
			structure = parser.get_structure(pdb_id_up, pdb_id_up+".pdb")
			modell = structure[0]
			io = bpdb.PDBIO()
			io.set_structure(modell)
			io.save('dna_'+pdb_id_up+'.pdb', ProtSelect())
			io.save('apo_'+pdb_id_up+'.pdb', ProtSelect1())
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('dna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			scmc=[]
			tot=[]
			bind=[]
			cn_count=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					time.sleep(5)
					print(" ")
					sys.stdout.flush()
					print(" ")
					df_inter = pd.DataFrame(df_inter)
					if len(df_inter.columns)>3:
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'interaction_energy.csv')
						for oo in list(set(energy_df['prot_resno'])):
							bind.append(i+"_"+str(int(oo)))				
						atoms_involve = inst2.f8_interaction_type(df_inter)
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						tot.append(df_at['count'].sum())
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						if 'CN' in df_at.index:
							cn_count.append(df_at.loc[['CN']].values)
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			#os.system("chmod -R +x {}".format(foldx))
			bind=list(set(bind))
			with open(r"apo_bind"+pdb_id_up+".pdb","w") as file1:
				with open(r"apo_"+pdb_id_up+".pdb") as file:
					for i in file.readlines():
						if i[0:4]=="ATOM":
							print(i)
							chain=i[21:26].split(" ")[0]
							print(chain)
							number=i[21:26].split(" ")[1]
							if chain+"_"+number in bind:
								file1.write(i)
			
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb="+pdb_id_up+".pdb --complexWithDNA=true", shell=True)
			p2.wait()
			
			fold_par=[]
			fold_val=[]
			cn_count=[]
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					#print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
			#os.system("chmod -R +x {}".format(naccess))
			
			p3=subprocess.Popen("./naccess "+pdb_id_up+".pdb -h", shell=True)
			p3.wait()
			rsa={}
			with open(pdb_id_up+".rsa") as file:
				for rows in file.readlines():
					if rows[0:3]=='RES':
						if not rows.split(" ")[2] in rsa:
							rsa[rows.split(" ")[2]]={}
						#print(rows.split(" "))
						rows=[str for str in rows.split(" ") if str.strip()]
						if not rows[2].isalpha():
								new_rows=rows[0:2]
								new_rows.append(rows[2][0])
								new_rows.append(rows[2][1:])
								for kkkk in rows[3:]:
									new_rows.append(kkkk)
								rows=new_rows						
						if not rows[2] in rsa:
							rsa[rows[2]]={}
						rsa[rows[2]][rows[3]]=float(rows[4])
			print(rsa)
			pdb_new=[]
			rsa_final=[]
			res={}
			rsa_dna=[]
			pdb_1=[]
			new_rsa=[]
			new_rsa_1={}
			new_rsa_1=[]
			for i in glob.glob("*interaction_energy.csv"):
					if len(rsa)>0:
						if not i.split("_")[1] in res:
							res[i.split("_")[1]]=[]
						if not i.split("_")[2] in res:
							res[i.split("_")[2]]=[]
						df=pd.read_csv(i)
						if 'prot_resno' in df.columns:
						 	for kk in range(len(df['prot_resno'])):
						 		if not str(int(df['prot_resno'][kk])) in res[i.split("_")[1]]:
						 			res[i.split("_")[1]].append(str(int(df['prot_resno'][kk])))
						 	for kk in range(len(df['na_resno'])):
						 		if not str(int(df['prot_resno'][kk])) in res[i.split("_")[2]]:
						 			res[i.split("_")[2]].append(str(int(df['prot_resno'][kk])))
						for j in res[i.split("_")[1]]:
							if str(j) in rsa[i.split("_")[1][0]]:
								new_rsa.append(rsa[i.split("_")[1][0]][str(j)])
							else:
								new_rsa.append(0)
						
						for j in res[i.split("_")[2]]:
							print(i)
							#print(j)
							print(rsa[x])
							if str(j) in rsa[i.split("_")[2][0]]:
								new_rsa_1.append(rsa[i.split("_")[2][0]][str(j)])
							else:
								new_rsa_1.append(0)
			rsa_protein=sum(new_rsa)
			rsa_dna=sum(new_rsa_1)
			aa = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M','PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
			df=pd.read_csv(r"potential_res.csv")
			c=df['pd'].tolist()
			d=df['potential_3.5'].tolist()
			pote={}
			pot=0
			append1=[]
			int_dict={}
			for i in glob.glob("*interaction_energy.csv"):
				print(i)
				#int_dict = {}
				df_inter=pd.read_csv(i)
				if len(df_inter.columns)>3:
					df_inter=df_inter.loc[df_inter['distance']<=3.5]
					df_inter=df_inter[['prot_res','prot_resno','na_res','na_resno']]
					#df_inter=df_inter.drop_duplicates()
					append1.append(df_inter)
				df_inter = pd.concat(append1)
				
			print(df_inter)
			df_inter=df_inter.drop_duplicates()
			print(df_inter)
			for index, dfrow in df_inter.iterrows():
				in_t = dfrow["prot_res"]+""+dfrow["na_res"]
				#print(in_t)
				if in_t in int_dict:
					int_dict[in_t] += 1
				else:
					int_dict.update({in_t:1})
			print(int_dict)
			for k in int_dict:
				print(pot)
				print(c)
				pot=pot+(int_dict[k]*d[c.index(k)])
			print(pot)
			
			cis_bond=fold_dict['cis_bond']
			
			cn1=sum(cn_count)/sum(tot)
			print(cn1, rsa_dna,rsa_protein,pot,cis_bond)
			predval=-3.57057422e+01*float(cn1)+-4.63008263e-04 *float(rsa_dna)+1.88277137e-04*float(rsa_protein)+1.99053887e-01*float(pot)+3.96877459e+02*float(cis_bond)-4.072444679861142
			predval="%.2f" % predval
			#if pdb_id_up=='input':
			#	resultout.write("User input")
			#else:
			#	resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval)+" ± 1.04")
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
			
##########################################################################################					Double-all-beta-nreg		  ##########################################################################################
##########################################################################################					Double-all-beta-nreg		  ##########################################################################################
##########################################################################################					Double-all-beta-nreg		  ##########################################################################################
##########################################################################################					Double-all-beta-nreg		  ##########################################################################################
	if model2=='ds' and model=='beta' and model1 =='nreg':
		with open("result.txt","w") as resultout:
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[0]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['DA','DG','DC','DT']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									if (distance<= 10):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('dna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			scmc=[]
			bind=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					time.sleep(5)
					print(" ")
					sys.stdout.flush()
					print(" ")
					if len(df_inter.columns)>3:	
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'_interaction_energy.csv')
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						atoms_involve = inst2.f8_interaction_type(df_inter)
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						for oo in list(set(df_inter['prot_resno'])):
							bind.append(i+"_"+str(int(oo)))
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
			print(res_bind_count)
			bind=list(set(bind))
			res_bind_uni=list(set(res_bind_count))
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb="+pdb_id_up+".pdb --complexWithDNA=true", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			cn_count=[]
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				print(f)
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
			Intra_clash_1=fold_dict['IntraclashesGroup1']
			#os.chdir("../")
			time.sleep(5)
			print(" ")
			sys.stdout.flush()
			print(" ")

			from Bio.PDB import PDBParser
			from Bio.PDB.DSSP import DSSP
			p = PDBParser()
			structure = p.get_structure("apo_"+pdb_id_up+".pdb", "apo_"+pdb_id_up+".pdb")
			modelp = structure[0]
			dssp = DSSP(modelp, "apo_"+pdb_id_up+".pdb")
			#print(dssp)
			from Bio.PDB.DSSP import dssp_dict_from_pdb_file
			dssp_tuple = dssp_dict_from_pdb_file("apo_"+pdb_id_up+".pdb")
			dssp_dict = dssp_tuple[0]
			print(bind)
			dssp_sec=[]
			data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
			for ds in bind:
				dssp_sec.append(data[dssp_dict[(ds.split("_")[0], (' ',int(ds.split("_")[1]), ' '))][1]])
			print(dssp_sec)					
			#dssp=dssp_sec
			#i[11:13].strip()+"_"+kk[1]
			#a_key = list(dssp.keys())
			#print(dssp[a_key])
			
			'''				
			os.chdir("../")
			p3=subprocess.Popen("./dssp -i apo_"+pdb_id_up+".pdb -o "+pdb_id_up+".dssp",shell=True)
			p3.wait()
			os.chdir(path)		
			shutil.copyfile("../"+pdb_id_up+".dssp", pdb_id_up+".dssp")		
			k=1
			
			data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
			dssp=[]
			with open(pdb_id_up+".dssp") as file:
				for i in file.readlines():
						if "#" in i:
								k=0
						if "!" in i:
								chain.append("-")
								ds.append("-")
						if k==0:
								if not "#" in i and not "!" in i:
										j=i.split(" ")
										#print( i[13:15].strip()+":"+i[0:11].strip())
										kk=i[0:11].strip()
										kk=[str for str in kk.split() if str.strip()]
										#print(kk)
										#print(f)
										#print(i[11:13].strip())
										if i[11:13].strip()+"_"+kk[1] in bind:
												dssp.append(data[i[16:18].strip()])
												chain=i[11:13].strip()+":"+kk[1]
												ds=data[i[16:18].strip()]
			'''
			dssp=dssp_sec
			tot=len(dssp_sec)
			he=dssp.count("helix")
			be=dssp.count("Beta")
			co=dssp.count("coil")
			percent_beta=100*(be/tot)
			print(percent_beta)	
			# polar asa
			
			p3=subprocess.Popen("./naccess "+pdb_id_up+".pdb", shell=True)
			p3.wait()
			rsa={}
			pol={'N':'polar','O':'polar','C':'nonpolar','S':'nonpolar'}
			i=pdb_id_up+".asa"
			if ".asa" in i and not "dna" in i and not "apo" in i:
				rsa[i.split(".")[0]]={}
				rsa[i.split(".")[0]]['polar']=0
				rsa[i.split(".")[0]]['nonpolar']=0
				with open(i) as file:
					for rows in file.readlines():
						if rows[0:4]=='ATOM':
							rows=[str for str in rows.split(" ") if str.strip()]
							if rows[4]+"_"+rows[5] in bind:			
								if rows[2][0] in pol:	
									rsa[i.split(".")[0]][pol[rows[2][0]]]=rsa[i.split(".")[0]][pol[rows[2][0]]]+float(rows[-2])
			polar_asa=[]					
			pdbapo=[]
			for i,j in rsa.items():
				polar_asa=j['polar']
				nonpolar_asa=j['nonpolar']
			
			
			time.sleep(5)
			print(" ")
			sys.stdout.flush()
			print(" ")
			os.system(r"rm *.asa")
			os.system(r"rm *.rsa")
			os.system(r"rm *.log")
			p3=subprocess.Popen("./naccess "+pdb_id_up+".pdb", shell=True)
			p3.wait()
			rsa={}
			with open(pdb_id_up+".rsa") as file:
				for rows in file.readlines():
					if rows[0:3]=='RES':
						if not rows.split(" ")[2] in rsa:
							rsa[rows.split(" ")[2]]={}
						#print(rows.split(" "))
						rows=[str for str in rows.split(" ") if str.strip()]
						if not rows[2].isalpha():
								new_rows=rows[0:2]
								new_rows.append(rows[2][0])
								new_rows.append(rows[2][1:])
								for kkkk in rows[3:]:
									new_rows.append(kkkk)
								rows=new_rows						
						if not rows[2] in rsa:
							rsa[rows[2]]={}
						rsa[rows[2]][rows[3]]=float(rows[4])
			res={}
			rsa_dna=[]
			new_rsa_1={}
			new_rsa_1=[]
			print(rsa)
			for i in glob.glob("*_atom_interaction.csv"):
					if len(rsa)>0:
						if not i.split("_")[1] in res:
							res[i.split("_")[1]]=[]
						if not i.split("_")[2] in res:
							res[i.split("_")[2]]=[]
						df=pd.read_csv(i)
						if 'prot_resno' in df.columns:
						 	for kk in range(len(df['na_resno'])):
						 		if not str(int(df['na_resno'][kk])) in res[i.split("_")[2]]:
						 			res[i.split("_")[2]].append(str(int(df['na_resno'][kk])))
						
						print(res)
						for j in res[i.split("_")[2]]:
							print(i)
							print(j)
							if str(j) in rsa[i.split("_")[2]]:
								new_rsa_1.append(rsa[i.split("_")[2]][str(j)])
							else:
								new_rsa_1.append(0)
			rsa_dna=sum(new_rsa_1)#complex rsa dna
			print(rsa_dna)
			time.sleep(5)
			print(" ")
			sys.stdout.flush()
			print(" ")
			p3=subprocess.Popen("./naccess dna_"+pdb_id_up+".pdb ", shell=True)
			p3.wait()
			rsa={}
			i="dna_{}.pdb".format(pdb_id_up)
			rsa[i.split(".")[0].split("_")[1]]={}
			
			with open("dna_"+pdb_id_up+".rsa") as file:
				for rows in file.readlines():
					#print(rows)
					if rows[0:3]=='RES':
						if not rows.split(" ")[2] in rsa[i.split(".")[0].split("_")[1]]:
							rsa[i.split(".")[0].split("_")[1]][rows.split(" ")[2]]={}
						rows=[str for str in rows.split(" ") if str.strip()]
						#print(rows)
						if not rows[2].isalpha():
							new_rows=rows[0:2]
							new_rows.append(rows[2][0])
							new_rows.append(rows[2][1:])
							for kkkk in rows[3:]:
								new_rows.append(kkkk)
							rows=new_rows
						#print(rows)						
						if not rows[2] in rsa[i.split(".")[0].split("_")[1]]:
							rsa[i.split(".")[0].split("_")[1]][rows[2]]={}
						rsa[i.split(".")[0].split("_")[1]][rows[2]][rows[3]]=float(rows[4])
						print(rsa)		
			res={}
			print(rsa)
			new_rsa_1={}
			for i in glob.glob("*_atom_interaction.csv"):
				x=i.split(".")[0]
				print(x)
				if x.split("_")[0] in rsa:
					print(x)
					if len(rsa[x.split("_")[0]])>0:
						print("i")
						if not i.split(".")[0].split("_")[1] in res:
							res[i.split(".")[0].split("_")[1]]=[]
						if not i.split(".")[0].split("_")[2] in res:
							res[i.split(".")[0].split("_")[2]]=[]
						df=pd.read_csv(i)
						print(res)
						if 'na_resno' in df.columns:
						 	for kk in range(len(df['na_resno'])):
						 		if not str(int(df['na_resno'][kk])) in res[i.split(".")[0].split("_")[2]]:
						 			res[i.split(".")[0].split("_")[2]].append(str(int(df['na_resno'][kk])))
						print(res)
						if not x.split("_")[0] in new_rsa_1:
						 	new_rsa_1[x.split("_")[0]]=[]
						print(new_rsa_1)
						for j in res[i.split(".")[0].split("_")[2]]:
							print(j)
							print(new_rsa_1)
							if str(j) in rsa[i.split(".")[0].split("_")[0]][i.split(".")[0].split("_")[2]]:
								new_rsa_1[x.split("_")[0]].append(rsa[i.split(".")[0].split("_")[0]][i.split(".")[0].split("_")[2]][str(j)])
							else:
								new_rsa_1[x.split("_")[0]].append(0)
			apo_rsa_dna=sum(new_rsa_1[x.split("_")[0]])
			print(apo_rsa_dna)
			del_dna=apo_rsa_dna-rsa_dna
			print(Intra_clash_1, polar_asa,nonpolar_asa,percent_beta,del_dna)
			
			predval=-0.17314104722099818*float(Intra_clash_1)+-0.015946034158376036 *float(polar_asa)+0.01208447830360034*float(nonpolar_asa)+0.07759484513559246*float(percent_beta)+-0.00010177351498883673*float(del_dna)-9.949351054556551
			predval="%.2f" % predval
			print(predval)
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval)+" ± 0.61")
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
##########################################################################################					Double-alpha-beta-reg		  ##########################################################################################
##########################################################################################					Double-alpha-beta-reg		 ##########################################################################################
##########################################################################################					Double-alpha-beta-reg		  ##########################################################################################
##########################################################################################					Double-alpha-beta-reg		  ##########################################################################################
	
	if model2=='ds' and model=='alphabeta' and model1 =='reg':	
		
		with open("result.txt","w") as resultout:
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[0]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['DA','DG','DC','DT']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						return int_df1, flag
					v=1
					for prot_res in prot_chain:
						print(prot_res)
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									if (distance<= 10):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('dna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			tot=[]
			mcmc=[]
			catom=0
			print(prot_chain)
			print(dna_chain)
			bind=[]
			scmc=[]
			sn_count=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					time.sleep(5)
					print(" ")
					sys.stdout.flush()
					print(" ")
					df_inter = pd.DataFrame(df_inter)
					if len(df_inter.columns)>3:
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'interaction_energy.csv')
						atoms_involve = inst2.f8_interaction_type(df_inter)
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						time.sleep(5)
						print("calculating features", flush="True")
						sys.stdout.flush()
						print(" ")
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						energy_dict = inst2.f9_energy_div(energy_df)
						df_en=pd.Series(energy_dict).to_frame()
						df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
						time.sleep(5)
						print(" ")
						sys.stdout.flush()
						print(" ")
						for p0 in df_inter['prot_resno'].tolist():
							res_bind_count.append(p0)
							bind.append(i+"_"+str(int(p0)))
						mcmc.append(df_en.loc[['mc_mc']].values)
						scmc.append(df_en.loc[['sc_mc']].values)
						for ppp in list(df_at.index.values):
							if ppp[0]=='C':
								catom=catom+df_at.loc[[ppp]].values
						#print(catom)
						tot.append(df_at['count'].sum())
						if 'SN' in df_at.index:
							sn_count.append(df_at.loc[['SN']].values)
	
			#print(res_bind_count)
			print(sn_count)
			print(scmc)
			#print(mcmc)
			time.sleep(5)
			print(" ")
			sys.stdout.flush()
			print(" ")
			res_bind_uni=list(set(res_bind_count))
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb="+pdb_id_up+".pdb --complexWithDNA=true", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				print(f)
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				print(fold_dict)
			total_residue=fold_dict['Number of Residues']
			bind_per=len(res_bind_uni)/int(total_residue)
			print(bind_per)
			#******************************************************************************bind greater than 0.25*****************************************************************	
			if bind_per>0.25:
				elec=fold_dict['Electrostatics']
				print(mcmc)
				mcmc1=sum(mcmc)[0]
				print(mcmc1)
				sn1=sum(sn_count)/sum(tot)
				p3=subprocess.Popen("./naccess "+pdb_id_up+".pdb -h", shell=True)
				p3.wait()
				rsa={}
				with open(pdb_id_up+".rsa") as file:
					for rows in file.readlines():
						if rows[0:3]=='RES':
							if not rows.split(" ")[2] in rsa:
								rsa[rows.split(" ")[2]]={}
							#print(rows.split(" "))
							rows=[str for str in rows.split(" ") if str.strip()]
							if not rows[2].isalpha():
									new_rows=rows[0:2]
									new_rows.append(rows[2][0])
									new_rows.append(rows[2][1:])
									for kkkk in rows[3:]:
										new_rows.append(kkkk)
									rows=new_rows						
							if not rows[2] in rsa:
								rsa[rows[2]]={}
							rsa[rows[2]][rows[3]]=float(rows[4])
				print(rsa)
				res={}
				rsa_dna=[]
				new_rsa=[]
				for i in glob.glob("*interaction_energy.csv"):
						if len(rsa)>0:
							if not i.split("_")[1] in res:
								res[i.split("_")[1]]=[]
							df=pd.read_csv(i)
							if 'prot_resno' in df.columns:
								for kk in range(len(df['prot_resno'])):
									if not df['prot_resno'][kk] in res[i.split("_")[1]]:
										res[i.split("_")[1]].append(df['prot_resno'][kk])
								if str(j) in rsa[i.split("_")[1][0]]:
									new_rsa.append(rsa[i.split("_")[1][0]][str(j)])
								else:
									new_rsa.append(0)
				rsa_protein=sum(new_rsa)
				print(rsa_protein)
				print("////////")
				p3=subprocess.Popen("./naccess "+pdb_id_up+".pdb", shell=True)
				p3.wait()
				x=pdb_id_up+".asa"
				rsa={}
				rsa_apo={}
				pol={'N':'polar','O':'polar','C':'nonpolar','S':'nonpolar'}
				i=pdb_id_up+".asa"
				if ".asa" in i and not "dna" in i and not "apo" in i:
					rsa[i.split(".")[0]]={}
					rsa[i.split(".")[0]]['polar']=0
					rsa[i.split(".")[0]]['nonpolar']=0
					with open(i) as file:
						for rows in file.readlines():
							if rows[0:4]=='ATOM':
								rows=[str for str in rows.split(" ") if str.strip()]
								if rows[4]+"_"+rows[5] in bind:			
									if rows[2][0] in pol:	
										rsa[i.split(".")[0]][pol[rows[2][0]]]=rsa[i.split(".")[0]][pol[rows[2][0]]]+float(rows[-2])
				print("done")
				time.sleep(5)
				print(" ")
				sys.stdout.flush()
				print(" ")
				p3=subprocess.Popen("./naccess apo_"+pdb_id_up+".pdb", shell=True)
				p3.wait()	
				time.sleep(5)
				print(" ")
				sys.stdout.flush()
				print(" ")
				i="apo_"+pdb_id_up+".asa"
				if ".asa" in i and not "dna" in i and "apo" in i:	
					with open(i) as file:
						#i=i.split("_")[1]
						rsa_apo[i.split(".")[0]]={}
						rsa_apo[i.split(".")[0]]['polar']=0
						rsa_apo[i.split(".")[0]]['nonpolar']=0
						for rows in file.readlines():
							if rows[0:4]=='ATOM':
								rows=[str for str in rows.split(" ") if str.strip()]
								if rows[4]+"_"+rows[5] in bind:			
									if rows[2][0] in pol:
										rsa_apo[i.split(".")[0]][pol[rows[2][0]]]=rsa_apo[i.split(".")[0]][pol[rows[2][0]]]+float(rows[-2])	
				print("done")
				pdbid=[]
				polar_asa=[]					
				nonpolar_asa=[]
				polar_apo_asa=[]
				nonpolar_apo_asa=[]
				pdbapo=[]
				time.sleep(5)
				print(" ")
				sys.stdout.flush()
				print(" ")
				for i,j in rsa.items():
					nonpolar_asa=j['nonpolar']
				for i,j in rsa_apo.items():
					nonpolar_apo_asa=j['nonpolar']
				print(nonpolar_apo_asa)
				print(nonpolar_asa)
				non_polar_asa_diff=nonpolar_apo_asa-nonpolar_asa
				#print(res_depth)
				print(catom)	
				#print(rsa_protein)
				print(non_polar_asa_diff)
				print(elec)
				print(mcmc1)
				print(sn1)
				'''
				def Convertst(string): 
					list1=[] 
					list1[:0]=string 
					return list1
				'''  
				from Bio.PDB.ResidueDepth import ResidueDepth
				from Bio.PDB.PDBParser import PDBParser
				from Bio.PDB.ResidueDepth import get_surface
				from Bio.PDB.ResidueDepth import min_dist
				from Bio.PDB.ResidueDepth import residue_depth
				'''				
				binres=binres5
				chains=[i.split("_")[1] for i in binres]
				for rs in chains:
					break
				res_dep=0
				print (rs)
				print (Convertst(rs))
				'''
				for j in prot_chain:
					#resnum=[int(i.split("_")[0][3:]) for i in binres if j==i.split("_")[1]]
					parser = PDBParser()
					structure = parser.get_structure("apo_"+pdb_id_up+".pdb", "apo_"+pdb_id_up+".pdb")
					model = structure[0]
					rd = ResidueDepth(model)
					surface = get_surface(model)
					rschain = model["{}".format(j)]
					resdepth=0
					ax=[]
					for a in bind:
						res = rschain[int(a.split("_")[1])]
						rd = residue_depth(res, surface)
						ax.append(rd)
					res_dep+=sum(ax)
				print(res_dep)	
				time.sleep(5)
				print(" ")
				sys.stdout.flush()
				print(" ")	
				print(catom, rsa_protein,non_polar_asa_diff,elec,mcmc1,sn1,res_dep)
				predval=-4.81026515e-04*float(catom)+2.25202998e-04*float(rsa_protein)+1.14774185e-03*float(non_polar_asa_diff)+6.78619000e-01*float(elec)+-7.45660148e-02*float(mcmc1)+-8.69409943e+02*float(sn1)+-8.561940404506892 +4.04921293e-04*float(res_dep)
				predval="%.2f" % predval
				if pdb_id_up=='input':
					resultout.write("User input")
				else:
					resultout.write(pdb_id_up)
				resultout.write("\n")
				#print (chain)
				prot_chain=",".join(prot_chain)
				resultout.write(str(prot_chain))
				resultout.write("\n")
				#print (binlig5)
				binlig5=",".join(dna_chain)
				resultout.write(str(binlig5))
				resultout.write("\n")
				resultout.write(str(predval)+" ± 1.01")
				disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
				resultout.write("\n")
				resultout.write(str(disass))

			#******************************************************************************bind less than 0.25*****************************************************************	
			if bind_per<=0.25:
					
					inter_clash_2=fold_dict['IntraclashesGroup2']
					scmc_energy=sum(scmc)
					urllib.request.urlretrieve(r"http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".out",pdb_id_up+"_summary.txt")
					k=0
					l=1
					with open(pdb_id_up+"_summary.txt") as file:
						for j in file.readlines():
							if "Local base-pair step parameters" in j or "Local base step parameters" in j:
								k=1
							if k==1:
								if "ave" in j or "a/a" in j:
									j=j.strip()
									s=[str for str in j.split(" ") if str.strip()]
									double_rise=s[3]
									double_twist=s[6]
									k=0
									l=0
						if l==1:
							print("error-3DNA-double")
					single_roll,single_twist=[],[]
					method_bind=1
					helix_dipole=fold_dict['helix dipole']
					k=0
					c=0
					l=1
					urllib.request.urlretrieve(r"http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.pars",pdb_id_up+"_bp_step.pars")
					with open(pdb_id_up+"_bp_step.pars") as file:
						for j in file.readlines():
							if "***local step parameters***" in j:
								k=1
							if k==1:
								if not "#" in j:
									j=j.strip()
									#print(j)
									s=[str for str in j.split(" ") if str.strip()]
									single_roll.append(float(s[5]))
									single_twist.append(float(s[6]))
									l=0
						if l==1:
							print(i+"error:single_base_step")
						roll1=sum(single_roll)/len(single_roll)
						twist1=sum(single_twist)/len(single_twist)
					print(inter_clash_2,double_rise,double_twist,roll1,twist1,scmc_energy,inter_clash_2)
					
					predval=-0.006348557963034027*float(inter_clash_2)+9.491377687782258 *float(double_rise)+0.1352193249119195*float(double_twist)+-0.062267425641198096*float(scmc_energy[0])+0.16710017174843358*float(roll1)+0.13799382696308116*float(twist1)-51.69011957297221
					predval="%.2f" % predval
					if pdb_id_up=='input':
						resultout.write("User input")
					else:
						resultout.write(pdb_id_up)
					resultout.write("\n")
					#print (chain)
					prot_chain=",".join(prot_chain)
					resultout.write(str(prot_chain))
					resultout.write("\n")
					#print (binlig5)
					binlig5=",".join(dna_chain)
					resultout.write(str(binlig5))
					resultout.write("\n")
					resultout.write(str(predval)+" ± 0.89")
					disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
					resultout.write("\n")
					resultout.write(str(disass))

##########################################################################################					Double-alpha-beta-nreg		  ##########################################################################################
##########################################################################################					Double-alpha-beta-nreg		 ##########################################################################################
##########################################################################################					Double-alpha-beta-nreg		  ##########################################################################################
##########################################################################################					Double-alpha-beta-nreg		  ##########################################################################################
	if model2=='ds' and model=='alphabeta' and model1 =='nreg':		
		with open("result.txt","w") as resultout:
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[0]
					#print(prot_chain.get_residues())
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['DA','DG','DC','DT']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						return int_df1, flag
					v=1
					#print(prot_chain.get_residues())
					#print(pp)
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									print(" ")
									distance = prot_atoms-rna_atoms
									if (distance<= 10):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						print(" ")
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							#if not 'H' in in_t:
							int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					print(" ")
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('dna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			bind=[]
			scmc=[]
			co_count=[]
			cn_count=[]
			tot=[]
			inter_energy=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					time.sleep(5)
					print(" ")
					sys.stdout.flush()
					print(" ")
					if len(df_inter.columns)>3:
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						print(" ")
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'interaction_energy.csv')
						if "total_energy" in energy_df:
							inter_energy.append(energy_df['total_energy'].sum())
							atoms_involve = inst2.f8_interaction_type(df_inter)
							df_at=pd.Series(atoms_involve).to_frame()
							df_at.columns=['count']
							time.sleep(5)
							print(" ")
							sys.stdout.flush()
							print(" ")
							df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
							energy_dict = inst2.f9_energy_div(energy_df)
							df_en=pd.Series(energy_dict).to_frame()
							df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
							for p0 in df_inter['prot_resno'].tolist():
								res_bind_count.append(p0)
								bind.append(i+"_"+str(int(p0)))
							scmc.append(df_en.loc[['sc_mc']].values)
							tot.append(df_at['count'].sum())
							if 'CO' in df_at.index:
								co_count.append(df_at.loc[['CO']].values)
							if 'CN' in df_at.index:
								cn_count.append(df_at.loc[['CN']].values)
							print(df_at)
							print(df_at.loc[['CN']])
			print(res_bind_count)
			co1=sum(co_count)/sum(tot)
			cn1=sum(cn_count)/sum(tot)
			bind=list(set(bind))
			res_bind_uni=list(set(res_bind_count))
			print(bind)
			print(co_count,cn_count)
			#print(res_bind_count)
			time.sleep(5)
			print(" ")
			sys.stdout.flush()
			print(" ")
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb="+pdb_id_up+".pdb --complexWithDNA=true", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				print(f)
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				print(fold_dict)
			total_residue=fold_dict['Number of Residues']
			bind_per=len(res_bind_uni)/int(total_residue)
			print(bind_per)		
			if bind_per>0.18:
				aa = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M','PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
				df=pd.read_csv(r"potential_res.csv")
				c=df['pd'].tolist()
				print(c)
				d=df['potential_6'].tolist()
				pote={}
				pot=0
				append1=[]
				int_dict={}
				for i in glob.glob("*interaction_energy.csv"):
					print(i)
					#int_dict = {}
					df_inter=pd.read_csv(i)
					if len(df_inter.columns)>3:
						df_inter=df_inter.loc[df_inter['distance']<=6]
						#df_inter=df_inter[['prot_res','prot_resno','na_res','na_resno']]
						#df_inter=df_inter.drop_duplicates()
						append1.append(df_inter)
					df_inter = pd.concat(append1)
				print(df_inter)
				df_inter=df_inter.drop_duplicates(['prot_res','prot_resno','na_res','na_resno'])#,'na_resno'])
				df_inter=df_inter.drop_duplicates()
				print(df_inter)
				for index, dfrow in df_inter.iterrows():
					in_t = dfrow["prot_res"]+" "+dfrow["na_res"]
					print(in_t)
					if in_t in int_dict:
						int_dict[in_t] += 1
					else:
						int_dict.update({in_t:1})
				print(int_dict)
				for k in int_dict:
					print(pot)
					pot=pot+(int_dict[k]*d[c.index(k)])
				print(pot)
				time.sleep(5)
				print(" ")
				sys.stdout.flush()
				print(" ")
				#shutil.copyfile("apo_"+pdb_id_up+".pdb", "../apo_"+pdb_id_up+".pdb")
				#os.chdir("../")
				from Bio.PDB import PDBParser
				from Bio.PDB.DSSP import DSSP
				p = PDBParser()
				structure = p.get_structure("apo_"+pdb_id_up+".pdb", "apo_"+pdb_id_up+".pdb")
				modelp = structure[0]
				dssp = DSSP(modelp, "apo_"+pdb_id_up+".pdb")
				#print(dssp)
				from Bio.PDB.DSSP import dssp_dict_from_pdb_file
				dssp_tuple = dssp_dict_from_pdb_file("apo_"+pdb_id_up+".pdb")
				dssp_dict = dssp_tuple[0]
				print(bind)
				time.sleep(5)
				print(" ")
				sys.stdout.flush()
				print(" ")
				dssp_sec=[]
				data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
				for ds in bind:
					dssp_sec.append(data[dssp_dict[(ds.split("_")[0], (' ',int(ds.split("_")[1]), ' '))][1]])
				print(dssp_sec)	
				'''
				p3=subprocess.Popen("./dssp -i apo_"+pdb_id_up+".pdb -o "+pdb_id_up+".dssp",shell=True)
				p3.wait()
				#os.system("rm apo_"+pdb_id_up+".pdb")
				#os.chdir(path)		
				#shutil.copyfile("../"+pdb_id_up+".dssp", pdb_id_up+".dssp")	
				os.system("rm ../apo_"+pdb_id_up+".dssp")	
				k=1
				data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','':'coil'}
				dssp=[]
				with open(pdb_id_up+".dssp") as file:
					for i in file.readlines():
						if "#" in i:
							k=0
						if "!" in i:
							chain.append("-")
							ds.append("-")
						if k==0:
							if not "#" in i and not "!" in i:
								j=i.split(" ")
								kk=i[0:11].strip()
								kk=[str for str in kk.split() if str.strip()]
								if i[11:13].strip()+"_"+kk[1] in bind:
										dssp.append(data[i[16:18].strip()])
										chain=i[11:13].strip()+":"+kk[1]
										ds=data[i[16:18].strip()]
				'''
				dssp=dssp_sec
				tot=len(dssp_sec)
				he=dssp.count("helix")
				be=dssp.count("Beta")
				co=dssp.count("coil")
				percent_coil=100*(co/tot)
				print(percent_coil)	
				# polar asa
				time.sleep(5)
				print(" ")
				sys.stdout.flush()
				print(" ")
				p3=subprocess.Popen("./naccess "+pdb_id_up+".pdb -h", shell=True)
				p3.wait()
				rsa={}
				pol={'N':'polar','O':'polar','C':'nonpolar','S':'nonpolar'}
				i=pdb_id_up+".asa"
				if ".asa" in i and not "dna" in i and not "apo" in i:
					rsa[i.split(".")[0]]={}
					rsa[i.split(".")[0]]['polar']=0
					rsa[i.split(".")[0]]['nonpolar']=0
					with open(i) as file:
						for rows in file.readlines():
							if rows[0:4]=='ATOM':
								rows=[str for str in rows.split(" ") if str.strip()]
								if rows[4]+"_"+rows[5] in bind:			
									if rows[2][0] in pol:	
										rsa[i.split(".")[0]][pol[rows[2][0]]]=rsa[i.split(".")[0]][pol[rows[2][0]]]+float(rows[-2])
				polar_asa=[]					
				pdbapo=[]
				for i,j in rsa.items():
					polar_asa=j['polar']
					nonpolar_asa=j['nonpolar']
				if pdb_id_up=="2hhu":
					polar_asa=0
					nonpolar_asa=0
				#complex_rsa_dna
				p3=subprocess.Popen("./naccess "+pdb_id_up+".pdb -h", shell=True)
				p3.wait()
				rsa={}
				with open(pdb_id_up+".rsa") as file:
					for rows in file.readlines():
						if rows[0:3]=='RES':
							if not rows.split(" ")[2] in rsa:
								rsa[rows.split(" ")[2]]={}
							#print(rows.split(" "))
							rows=[str for str in rows.split(" ") if str.strip()]
							if not rows[2].isalpha():
									new_rows=rows[0:2]
									new_rows.append(rows[2][0])
									new_rows.append(rows[2][1:])
									for kkkk in rows[3:]:
										new_rows.append(kkkk)
									rows=new_rows						
							if not rows[2] in rsa:
								rsa[rows[2]]={}
							rsa[rows[2]][rows[3]]=float(rows[4])
				res={}
				rsa_dna=[]
				new_rsa_1={}
				new_rsa_1=[]
				for i in glob.glob("*interaction_energy.csv"):
						if len(rsa)>0:
							if not i.split("_")[1] in res:
								res[i.split("_")[1]]=[]
							if not i.split("_")[2] in res:
								res[i.split("_")[2]]=[]
							df=pd.read_csv(i)
							if 'prot_resno' in df.columns:
							 	for kk in range(len(df['na_resno'])):
							 		if not str(int(df['na_resno'][kk])) in res[i.split("_")[2]]:
							 			res[i.split("_")[2]].append(str(int(df['na_resno'][kk])))
							
							for j in res[i.split("_")[2]]:
								print(i)
								#print(j)
								print(rsa[x])
								if str(j) in rsa[i.split("_")[2][0]]:
									new_rsa_1.append(rsa[i.split("_")[2][0]][str(j)])
								else:
									new_rsa_1.append(0)
				rsa_dna=sum(new_rsa_1)#complex rsa dna
				#polar asa#
				k=0
				l=1
				urllib.request.urlretrieve(r"http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".outs",pdb_id_up+"_summary.txt")
				with open(pdb_id_up+"_summary.txt") as file:
					for j in file.readlines():
						if "Local base-pair step parameters" in j or "Local base step parameters" in j:
							k=1
						if k==1:
							if "ave" in j or "a/a" in j:
								j=j.strip()
								s=[str for str in j.split(" ") if str.strip()]
								double_roll=s[5]
								k=0
								l=0
					if l==1:
						print("error-3DNA-double")
				single_slide,single_roll=[],[]
				k=0
				c=0
				l=1
				urllib.request.urlretrieve(r"http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.pars",pdb_id_up+"_bp_step.pars")
				with open(pdb_id_up+"_bp_step.pars") as file:
					for j in file.readlines():
						if "***local step parameters***" in j:
							k=1
						if k==1:
							if not "#" in j:
								j=j.strip()
								#print(j)
								s=[str for str in j.split(" ") if str.strip()]
								single_slide.append(float(s[2]))
								single_roll.append(float(s[5]))
								l=0
					if l==1:
						print(i+"error:single_base_step")
					single_slide=sum(single_slide)/len(single_slide)
					roll1=sum(single_roll)/len(single_roll)
					print(pot, percent_coil,double_roll, roll1, single_slide, rsa_dna,polar_asa)
					predval=1.64427115e-02*float(pot)+5.22101243e-02*float(percent_coil)+2.97321993e-01*float(double_roll)-1.85630247e-01*float(roll1)+-2.02462842e+00*float(single_slide)+2.98758486e-04*float(rsa_dna)+-1.10121795e-03*float(polar_asa)-11.64571271247264
				predval="%.2f" % predval
				if pdb_id_up=='input':
					resultout.write("User input")
				else:
					resultout.write(pdb_id_up)
				resultout.write("\n")
				#print (chain)
				prot_chain=",".join(prot_chain)
				resultout.write(str(prot_chain))
				resultout.write("\n")
				#print (binlig5)
				binlig5=",".join(dna_chain)
				resultout.write(str(binlig5))
				resultout.write("\n")
				resultout.write(str(predval)+" ± 1.53")
				disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
				resultout.write("\n")
				resultout.write(str(disass))
#******************************************************************************bind less than 0.25*****************************************************************		
			if bind_per<=0.18:
				#co1=sum(co_count)
				#cn1=sum(cn_count)
				print(cn1)
				fold_van=fold_dict['Van der Waals']
				fold_sol=fold_dict['Solvation Polar']
				fold_ene=fold_dict['energy Ionisation']
				fold_back=fold_dict['Backbone Hbond']
				inter1=sum(inter_energy)
				print(inter1)
				print(fold_back)
				print(fold_ene)
				print(cn1)
				print(co1)
				predval=50.4776234*float(co1[0])+0.49956648*float(fold_van)+-0.27422254*float(fold_sol)-0.17749198*float(inter1)+11.72404768*float(fold_ene)+33.32381432*float(cn1[0])-0.15405359*float(fold_back)-24.157688458017326
				predval="%.2f" % predval
				if pdb_id_up=='input':
					resultout.write("User input")
				else:
					resultout.write(pdb_id_up)
				resultout.write("\n")
				#print (chain)
				prot_chain=",".join(prot_chain)
				resultout.write(str(prot_chain))
				resultout.write("\n")
				#print (binlig5)
				binlig5=",".join(dna_chain)
				resultout.write(str(binlig5))
				resultout.write("\n")
				resultout.write(str(predval)+" ± 1.32")
				disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
				resultout.write("\n")
				resultout.write(str(disass))
##########################################################################################					Double-others		  ##########################################################################################
##########################################################################################					Double-others		##########################################################################################
##########################################################################################					Double-others		  ##########################################################################################
##########################################################################################					Double-others		  ##########################################################################################
	if model2=='ds' and model=='other':
		with open("result.txt","w") as resultout:
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[0]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['DA','DG','DC','DT']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									if (distance<= 10):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					#'NP':0,'CP':0,'CC':0,'NC':0,'OP':0	OO	CN	NN	SN	SC	SO	SP

					#print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('dna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			scmc=[]
			co_count=[]
			sn_count=[]
			inter_energy=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
					time.sleep(5)
					print(" ")
					sys.stdout.flush()
					print(" ")
					if len(df_inter.columns)>3:
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'_interaction_energy.csv')
						inter_energy.append(energy_df['total_energy'].sum())
						time.sleep(5)
						print(" ")
						sys.stdout.flush()
						print(" ")
						atoms_involve = inst2.f8_interaction_type(df_inter)
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						#energy_dict = inst2.f9_energy_div(energy_df)
						#df_en=pd.Series(energy_dict).to_frame()
						#df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
						for p0 in df_inter['prot_resno'].tolist():
							res_bind_count.append(p0)
						#scmc.append(df_en.loc[['scmc']].values)
						if 'SN' in df_at.index:
							sn_count.append(df_at.loc[['SN']].values)
					#cn_count.append(df_at.loc[['CN']].values/df_at['count'].sum())

			#print(res_bind_count)
			if len(sn_count)==0:
				sn1=0
			else:
				sn1=sum(sn_count)
			res_bind_uni=list(set(res_bind_count))
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb="+pdb_id_up+".pdb --complexWithDNA=true", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			time.sleep(5)
			print(" ")
			sys.stdout.flush()
			print(" ")
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				#print(f)
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					#print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				#print(fold_dict)
				shutil.copyfile("../3vvv.py","3vvv.py")
				shutil.copyfile("../vol_out.txt","vol_out.txt")
				p=subprocess.Popen("/opt/websites/anaconda/bin/python3.7 3vvv.py apo_"+pdb_id_up,shell=True)
				p.wait()
				with open("vol_out.txt") as file:
					g=file.readlines()
					vol1=g[0]
					surface_area=g[1]
					#for vol in file.readlines():
					#	vol1=vol					
					#	surface_area=vol[1]
			'''	
			shutil.copyfile("apo_"+pdb_id_up+".pdb","../apo_"+pdb_id_up+".pdb")
			os.system("chmod -R 777 ../apo_"+pdb_id_up+".pdb")
			os.chdir("../")
			p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".pdb > volume",shell=True)
			p4.wait()
			os.chdir(path)
			shutil.copyfile("../volume", "volume")
			with open("volume") as file:
				for vol in file.readlines():
					vol1=vol.split()[2]					
					surface_area=vol.split()[3]
			'''		
			#sn1=sum(sn_count)
			# apo_rsa_dna
			p3=subprocess.Popen("./naccess dna_"+pdb_id_up+".pdb ", shell=True)
			p3.wait()
			rsa={}
			i="dna_{}.pdb".format(pdb_id_up)
			rsa[i.split(".")[0].split("_")[1]]={}
			
			with open("dna_"+pdb_id_up+".rsa") as file:
				for rows in file.readlines():
					#print(rows)
					if rows[0:3]=='RES':
						if not rows.split(" ")[2] in rsa[i.split(".")[0].split("_")[1]]:
							rsa[i.split(".")[0].split("_")[1]][rows.split(" ")[2]]={}
						rows=[str for str in rows.split(" ") if str.strip()]
						#print(rows)
						if not rows[2].isalpha():
							new_rows=rows[0:2]
							new_rows.append(rows[2][0])
							new_rows.append(rows[2][1:])
							for kkkk in rows[3:]:
								new_rows.append(kkkk)
							rows=new_rows
						#print(rows)						
						if not rows[2] in rsa[i.split(".")[0].split("_")[1]]:
							rsa[i.split(".")[0].split("_")[1]][rows[2]]={}
						rsa[i.split(".")[0].split("_")[1]][rows[2]][rows[3]]=float(rows[4])
						print(rsa)		
			res={}
			print(rsa)
			new_rsa_1={}
			for i in glob.glob("*_atom_interaction.csv"):
				x=i.split(".")[0]
				print(x)
				if x.split("_")[0] in rsa:
					print(x)
					if len(rsa[x.split("_")[0]])>0:
						print("i")
						if not i.split(".")[0].split("_")[1] in res:
							res[i.split(".")[0].split("_")[1]]=[]
						if not i.split(".")[0].split("_")[2] in res:
							res[i.split(".")[0].split("_")[2]]=[]
						df=pd.read_csv(i)
						print(res)
						if 'na_resno' in df.columns:
						 	for kk in range(len(df['na_resno'])):
						 		if not str(int(df['na_resno'][kk])) in res[i.split(".")[0].split("_")[2]]:
						 			res[i.split(".")[0].split("_")[2]].append(str(int(df['na_resno'][kk])))
						print(res)
						if not x.split("_")[0] in new_rsa_1:
						 	new_rsa_1[x.split("_")[0]]=[]
						print(new_rsa_1)
						for j in res[i.split(".")[0].split("_")[2]]:
							print(j)
							print(new_rsa_1)
							if str(j) in rsa[i.split(".")[0].split("_")[0]][i.split(".")[0].split("_")[2]]:
								new_rsa_1[x.split("_")[0]].append(rsa[i.split(".")[0].split("_")[0]][i.split(".")[0].split("_")[2]][str(j)])
							else:
								new_rsa_1[x.split("_")[0]].append(0)
			apo_rsa_dna=sum(new_rsa_1[x.split("_")[0]])
			time.sleep(5)
			print(" ")
			sys.stdout.flush()
			print(" ")
			inter_clash_1=fold_dict['IntraclashesGroup1']
			inter_foldx=fold_dict['Interaction Energy']
			solv_hydro=fold_dict['Solvation Hydrophobic']
			#****surface_area***
			#***apo_rsa_dna***
			print(inter_clash_1,inter_foldx,solv_hydro,surface_area,sn1,apo_rsa_dna)
			predval=-9.98469859e-02*float(inter_clash_1)+1.22810352e-01*float(inter_foldx)+-1.62610279e-01*float(solv_hydro)-7.49730483e-05*float(surface_area)-3.61572487e-02*float(sn1)+-7.84279169e-05*float(apo_rsa_dna)-9.390454617677126
			predval="%.2f" % predval
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval)+" ± 0.83")
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))

'''
		# ######################################################		Residue Depth			  ######################################################
			def Convertst(string): 
				list1=[] 
				list1[:0]=string 
				return list1  
			from Bio.PDB.ResidueDepth import ResidueDepth
			from Bio.PDB.PDBParser import PDBParser
			from Bio.PDB.ResidueDepth import get_surface
			from Bio.PDB.ResidueDepth import min_dist
			from Bio.PDB.ResidueDepth import residue_depth
			binres=binres5
			chains=[i.split("_")[1] for i in binres]
			for rs in chains:
				break
			res_dep=0
			print (rs)
			print (Convertst(rs))
			for j in Convertst(rs):
				resnum=[int(i.split("_")[0][3:]) for i in binres if j==i.split("_")[1]]
				parser = PDBParser()
				structure = parser.get_structure("{0}", "{0}_atom.pdb".format(pdb_id_up))
				model = structure[0]
				rd = ResidueDepth(model)
				surface = get_surface(model)
				rschain = model["{}".format(j)]
				resdepth=0
				ax=[]
				for a in resnum:
					res = rschain[a]
					rd = residue_depth(res, surface)
					ax.append(rd)
				res_dep+=sum(ax)
			predval= 0.20665088*res_dep-0.00410831*num-0.22996925*sum(electro)-0.00685319*sum(energy_SolvP)-4.093052692485966
			predval="%.2f" % predval
			resultout.write("\n")
			resultout.write(str(predval)+" ± 0.529")
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
'''

########## common for all the classifications - start


timetaken=timeit.default_timer() - start

with open("result.py", "w") as polyout:
	print("ii")
	polyout.write("""#!/opt/websites/anaconda/envs/tf37/bin/python\nimport cgi\nimport cgitb; cgitb.enable()\nprint ('Content-Type: text/html\\r\\n')\n""")
	g=open("index.txt").readlines()
	for gg in g:
		gg=gg.rstrip()
		polyout.write("""print (\"\"\"{}\"\"\")""".format(gg))
		polyout.write("\n")
	f=open("result.txt").readlines()
	# polyout.write ("print ('<center>Time taken for the calculation: {:.2f} seconds</center>')".format(timetaken))
	# polyout.write("\n")
	# polyout.write("print ('<br/>')")
	# polyout.write("\n")
	# polyout.write ("print ('<center>Model selected: {}</center>')".format(modeltype))
	# polyout.write("\n")
	# polyout.write("print ('<br/>')")
	# polyout.write("\n")


	polyout.write("""print ('<font face="Times New Roman" ><table id="customers">')""")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td><b>Time taken</b></td><td>{:.2f} seconds</td></td>')".format(timetaken))
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td><b>Model selected</b></td><td>"+modeltype+"-"+modeltype1+"-"+modeltype2+"</td>')")
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('</table></font>')")
	polyout.write("\n")
	polyout.write("print ('<br/>')")
	polyout.write("\n")
	polyout.write("print ('<br/>')")
	polyout.write("\n")



	polyout.write("""print ('<font face="Times New Roman" ><table id="customers">')""")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td><b>PDB ID</b></td><td><b>protein Chain(s)</b></td><td><b>DNA chain (s) </b></td><td><b>Predicted &Delta;G (kcal/mol)</b></td><td><b>K<sub>d</sub> (M)</b></td>')")
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>')".format(f[0].rstrip(),f[1].rstrip(),f[2].rstrip(),f[3].rstrip(),f[4].rstrip()))
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('</table></font>')")
	polyout.write("\n")
	h=open("footer.txt").readlines()
	for hh in h:
		hh=hh.rstrip()
		polyout.write("""print (\"\"\"{}\"\"\")""".format(hh))
		polyout.write("\n")
os.system("chmod +x result.py")
os.system("ls > remfile")
'''
redirectURL = "/cgi-bin/%s/result.py" % randname
print ('Content-Type: text/html\r\n')
print ('<html>')
print ('  <head>')
print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
print ('	<title>You are going to be redirected</title>')
'''
'''
redirectURL = "%s/result.py" % randname
print ('Content-Type: text/html\r\n')
print ('<html>')
print ('  <head>')
print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
print('	<title>You are going to be redirected</title>')
'''

remfile=open("remfile").readlines()
for rf in remfile:
	rf=rf.rstrip()
	if rf =='molecules':
		os.rmdir('{}'.format(rf))
	elif rf =='autodock':
		shutil.rmtree('autodock')
	elif rf not in ['result.txt','result.py',"style4.css","vol_out.txt"]:
		os.remove('{}'.format(rf))

#print ('	Redirecting... <a href="%s">Click here if you are not redirected</a>' % redirectURL)
print ('</body>')
print ('</html>')

########## common for all the classifications - end