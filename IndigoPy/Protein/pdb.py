#!/usr/bin/python3
# -*- coding: UTF-8-*-
from urllib import request
from ..Math import coordinates as cdn
import copy
import os
t3t1={'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','PHE':'F','PRO':'P','SER':'S','THR':'T','HIS':'H','TRP':'W','CYS':'C','ASP':'D','GLU':'E','LYS':'K','TYR':'Y','MET':'M','ASN':'N','GLN':'Q','ARG':'R'}




#========================================================================================================================
#===========||====||||||==||||||==||||====||||||==||||||==||||||==||==||====||======||====||==||==||======||==||=========
#=========||||||==||||====||======||==||==||||||==||||||==||======||||||====||======||====||||====||======||||||=========
#=========||==||==||||||==||||||==||||====||||||==||======||||====||==||====||====||||====||==||==||||||==||||||=========
#========================================================================================================================
#===============================================Indigo Mad  vesion 0.01a1================================================
#========================================================================================================================
#=========||||||==||||||==||||======||||==||||======||||==||||||==||==||==||==||==||||||==||==||==||==||==||||===========
#=========||==||==||==||==||||======||||==||||======||======||====||==||==||==||==||||||====||======||======||===========
#=========||==||==||||||==||==========||==||==||==||||======||====||||||====||====||==||==||==||====||======||||=========
#========================================================================================================================




#========================================================================================================================
#=====================================================||||||||||||||=====================================================
#=====================================================||=================================================================
#=====================================================||=================================================================
#=====================================================||=================================================================
#=====================================================||=================================================================
#=====================================================||=================================================================
#=====================================================||||||||||||||=====================================================
#========================================================================================================================
#===================================================Classes definition===================================================

#========================================================================================================================
#=========================================================||||===========================================================
#=========================================================||||===========================================================
#=========================================================||=============================================================
#========================================================================================================================
#=======================================Basic position informations for a protein========================================
class Atom():
	def __init__(self,name,pos,ele):
		self.name=name
		self.pos=cdn.Position(pos)
		self.ele=ele
	def __str__(self):
		return self.name+" "+self.ele+str(self.pos)
	__repr__=__str__
	def __getitem__(self, key):
		if isinstance(key,str):
			if key=='x' or key=='X':
				return self.pos[0]
			elif key=='y' or key=='Y':
				return self.pos[1]
			elif key=='y' or key=='Y':
				return self.pos[2]
			else:
				raise KeyError("Expecting an intege or 'x','X','y','Y','Z','z' but a "+str(type(key))+"is given")
		elif isinstance(key,int):
			return self.pos[key]
		elif isinstance(key,slice):
			return self.pos[key]
		else:
			raise KeyError("Expecting an intege or 'x','X','y','Y','Z','z' but a "+str(type(key))+"is given")
	def __setitem__(self, key, value):
		if isinstance(key,str):
			if key=='x' or key=='X':
				self.pos[0]=value
			elif key=='y' or key=='Y':
				self.pos[1]=value
			elif key=='y' or key=='Y':
				self.pos[2]=value
			else:
				raise KeyError("Expecting an intege or 'x','X','y','Y','Z','z' but a "+str(type(key))+"is given")
		elif isinstance(key,int):
			self.pos[key]=value
		elif isinstance(key,slice):
			self.pos[key]=value
		else:
			raise KeyError("Expecting an intege or 'x','X','y','Y','Z','z' but a "+str(type(key))+"is given")






class Residue():
	def __init__(self,name,num,atomlist):
		self.name=name
		self.atoms=atomlist
		self.num=num
	def __getitem__(self, key):
		if isinstance(key,str):
			temlist=[]
			for i in self.atoms:
				if i.name==key:
					return i
					break
		elif isinstance(key,int):
			return self.atoms[key]
		elif isinstance(key,slice):
			return self.atoms[key]
		else:
			raise KeyError("Expecting an intege or a string but a "+str(type(key))+"is given")
	def __setitem__(self, key, value):
		if isinstance(key,str):
			temlist=[]
			for i in range(len(self.atoms)):
				if self.atoms[i].name==key:
					self.atoms[i]=value
					break
		elif isinstance(key,int):
			self.atoms[key]=value
		elif isinstance(key,slice):
			self.atoms[key]=value
		else:
			raise KeyError("Expecting an intege or a string but a "+str(type(key))+"is given")
	def __str__(self):
		return str(self.num)+" "+self.name
	__repr__=__str__
	def __len__(self):
		return len(self.atoms)
	def __iter__(self):
		return iter(self.atoms)










class Chain():
	def __init__(self,name,residueslist):
		self.name=name
		self.residues=residueslist
		self.seq=''
	def __getitem__(self, key):
		if isinstance(key,str):
			temlist=[]
			for i in self.residues:
				if i.name==key.upper():
					temlist.append(i)
				if i.name==t3t1.get(key):
					temlist.append(i)
			return temlist
		elif isinstance(key,int):
			for i in self.residues:
				if i.num==key:
					return i
					break
		elif isinstance(key,slice):
			if key.start!=None:
				if key.stop!=None:
					return self.residues[self.residues.index(self[key.start]):self.residues.index(self[key.stop]):key.step]
				else:
					return self.residues[self.residues.index(self[key.start])::key.step]
			else:
				if key.stop!=None:
					return self.residues[:self.residues.index(self[key.stop]):key.step]
				else:
					return self.residues[::key.step]

		else:
			raise KeyError("Expecting an intege or a string but a "+str(type(key))+"is given")

	def __setitem__(self, key, value):
		if type(key)==int:
			for i in range(len(self.residues)):
				if self.residues[i].num==key:
					self.residues[i].num=value
					break
			for i in self:
				try:
					self.seq=self.seq+t3t1[i.name]
				except:
					print('[Warning]:Unknow residues '+i.name+' included.')
		elif isinstance(key,slice):
			raise KeyError("Do not use slice to set values in a Chain")
		else:
			raise KeyError("Expecting an intege but a "+str(type(key))+"is given")
	def __str__(self):
		return self.name+"-"+self.seq
	__repr__=__str__
	def __iter__(self):
		return iter(self.residues)
	def __len__(self):
		return len(self.residues)
	def seqlen(self):
		return len(self.seq)
	def seqence(self, beg, end):
		seqtem=''
		for i in range(beg,end+1):
			try:
				seqtem=seqtem+t3t1[self[i].name]
			except:
				print('[Warming]:Unknow residue '+chaintem[i].name+' included.')
		return seqtem
	def __setattr__(self,name,value):
		if name == "residues":
			self.__dict__[name] = value
			temstr=''
			for i in self:
				try:
					temstr=temstr+t3t1[i.name]
				except:
					print('[Warning]:Unknow residues '+i.name+' included.')
			self.seq=temstr
		else:
			self.__dict__[name] = value










class Model():
	def __init__(self,chainlist):
		self.chains=chainlist
	def __getitem__(self, key):
		if isinstance(key,str):
			temlist=[]
			for i in self.chains:
				if i.name==key:
					return i
					break
		elif isinstance(key,int):
			return self.chains[key]
		elif isinstance(key,slice):
			return self.chains[key]
		else:
			raise KeyError("Expecting an intege or a string but a "+str(type(key))+"is given")
	def __setitem__(self, key, value):
		if isinstance(key,str):
			temlist=[]
			for i in range(len(self.chains)):
				if self.chains[i].name==key:
					self.chains[i]=value
					break
		elif isinstance(key,int):
			self.chains[key]=value
		elif isinstance(key,slice):
			self.chains[key]=value
		else:
			raise KeyError("Expecting an intege or a string but a "+str(type(key))+"is given")
	def __iter__(self):
		return iter(self.chains)
	def __len__(self):
		return len(self.chains)
	def seqlen(self):
		temnum=0
		for i in self:
			temnum+=len(i)
		return temnum







#========================================================================================================================
#===========================================================||||=========================================================
#===========================================================||===========================================================
#=========================================================||||===========================================================
#========================================================================================================================
#==================================================Secondary structures==================================================



"""
				TYPE OF HELIX          CLASS NUMBER 
									   (COLUMNS 39 - 40)
	  ----------------------------------------------
	  Right-handed alpha (default)       1
	  Right-handed omega                 2
	  Right-handed pi                    3
	  Right-handed gamma                 4
	  Right-handed 310                   5
	  Left-handed alpha                  6
	  Left-handed omega                  7
	  Left-handed gamma                  8
	  27 ribbon/helix                    9
	  Polyproline                       10
"""

class Helix():
	def __init__(self,pdbstr):
		self.ini=int(pdbstr[21:25])
		self.inic=pdbstr[19]
		self.end=int(pdbstr[33:37])
		self.endc=pdbstr[31]
		self.type=int(pdbstr[38:40])
	def __str__(self):
		return "Helix"+str(len(self))+" "+self.inic+str(self.ini)+"-"+self.endc+str(self.end)
	def __len__(self):
		return self.end-self.ini+1
	def seqlen(self):
		return self.end-self.ini+1


'''
Sense of strand with respect to previous strand in the sheet. 
0 if first strand, 
1 if parallel,
-1 if anti-parallel.
'''
class Strand():
	def __init__(self,pdbstr):
		self.ini=int(pdbstr[22:26])
		self.inic=pdbstr[21]
		self.end=int(pdbstr[33:37])
		self.endc=pdbstr[32]
		self.type=int(pdbstr[38:40])
	def __str__(self):
		return "Strand"+str(len(self))+" "+self.inic+str(self.ini)+"-"+self.endc+str(self.end)
	def __len__(self):
		return self.end-self.ini+1
	def seqlen(self):
		return self.end-self.ini+1

class Sheet():
	def __init__(self,Protein,pdblist):
		self.pro=Protein
		self.strands=pdblist
		self.link=[]
		for i in range(len(self.strands)):
			if i==0:
				self.strands[i]=Strand(self.strands[i])
			else:
				self.link.append(((self.strands[i][64],int(self.strands[i][65:69]),self.strands[i][56:60].replace(' ','')),(self.strands[i][49],int(self.strands[i][50:54]),self.strands[i][41:45].replace(' ',''))))
				self.strands[i]=Strand(self.strands[i])

	def __str__(self):
		tem=''
		temlist=[]
		temnum=len(self.strands)-1
		for i in range(temnum+1):
			if i==0:
				temlist.append([self.strands[i].type,None,self.link[i][0][1],self.strands[i].inic,self.strands[i].ini,self.pro.seq(self.strands[i].inic,self.strands[i].ini,self.strands[i].end)])
			elif i==temnum:
				temlist.append([self.strands[i].type,self.link[i-1][1][1],None,self.strands[i].inic,self.strands[i].ini,self.pro.seq(self.strands[i].inic,self.strands[i].ini,self.strands[i].end)])
			else:
				temlist.append([self.strands[i].type,self.link[i-1][1][1],self.link[i][0][1],self.strands[i].inic,self.strands[i].ini,self.pro.seq(self.strands[i].inic,self.strands[i].ini,self.strands[i].end)])
		for i in range(temnum+1):
			if i==0:
				temlist[i][0]='>'
			else:
				if temlist[i][0]==1:
					temlist[i][0]=temlist[i-1][0]
				else:
					if temlist[i-1][0]=='>':
						temlist[i][0]='<'
					else:
						temlist[i][0]='>'
		for i in range(temnum+1):
			if i==0:
				temlist[i].append(str(temlist[i][4])+temlist[i][3]+temlist[i][0])
				temlist[i].append(None)
				temlist[i].append(len(temlist[i][6])+temlist[i][2]-temlist[i][4])
				temlist[i][6]=temlist[i][6]+temlist[i][5]
			elif i==temnum:
				if temlist[i][0]=='>':
					temlist[i].append(str(temlist[i][4])+temlist[i][3]+temlist[i][0])
					temlist[i].append(len(temlist[i][6])+temlist[i][1]-temlist[i][4])
					temlist[i].append(None)
					temlist[i][6]=temlist[i][6]+temlist[i][5]
				else:
					temlist[i].append(temlist[i][0]+str(temlist[i][4])+temlist[i][3])
					temlist[i].append(len(temlist[i][6])+temlist[i][1]-temlist[i][4])
					temlist[i].append(None)
					temlist[i][6]=temlist[i][5][::-1]+temlist[i][6]
					temlist[i][7]=len(temlist[i][6])-temlist[i][7]-1
				gap=temlist[i-1][8]-temlist[i][7]
				if gap>0:
					temlist[i][6]=' '*gap+temlist[i][6]
					temlist[i][7]+=gap
				else:
					gap=-gap
					for j in range(i):
						temlist[j][6]=' '*gap+temlist[j][6]
			else:
				if temlist[i][0]=='>':
					temlist[i].append(str(temlist[i][4])+temlist[i][3]+temlist[i][0])
					temlist[i].append(len(temlist[i][6])+temlist[i][1]-temlist[i][4])
					temlist[i].append(len(temlist[i][6])+temlist[i][2]-temlist[i][4])
					temlist[i][6]=temlist[i][6]+temlist[i][5]
				else:
					temlist[i].append(temlist[i][0]+str(temlist[i][4])+temlist[i][3])
					temlist[i].append(len(temlist[i][6])+temlist[i][1]-temlist[i][4])
					temlist[i].append(len(temlist[i][6])+temlist[i][2]-temlist[i][4])
					temlist[i][6]=temlist[i][5][::-1]+temlist[i][6]
					temlist[i][7]=len(temlist[i][6])-temlist[i][7]-1
					temlist[i][8]=len(temlist[i][6])-temlist[i][8]-1
				gap=temlist[i-1][8]-temlist[i][7]
				if gap>0:
					temlist[i][6]=' '*gap+temlist[i][6]
					temlist[i][7]+=gap
					temlist[i][8]+=gap
				else:
					gap=-gap
					for j in range(i):
						temlist[j][6]=' '*gap+temlist[j][6]
		for i in temlist:
			tem=tem+i[6]+'\n'

		return tem
	def __len__(self):
		return len(self.strands)

	def __getitem__(self, key):
		return self.strands[key]










#========================================================================================================================
#=========================================================||==||=========================================================
#=========================================================||||||=========================================================
#=========================================================||||||=========================================================
#========================================================================================================================
#==============================================Main classes describing a protien=========================================

#========================================================================================================================
#=====================================================A Uniform Protein =================================================
#========================================================================================================================
class Protein():
	def __init__(self,pdbid=None,**kargs):
		if pdbid!=None:#It's a pdb ID, download it from the website.
			url='https://files.rcsb.org/download/'+pdbid+'.pdb'
			with request.urlopen(url) as f:
				if f.status== 200:
					f=f.read().decode()
				else:
					raise Exception("Failed to Download"+pdbid+'.pdb. Error '+f.status+' '+f.reason)
			try:
				self.name=kargs['name']
			except:
				self.name=pdbid
		else:
			try:
				f=kargs['file']
				try:
					self.name=kargs['name']
				except:
					self.name=None
					print('[Warning]Can not find the name of this protein, which can given by the argument: name=[the name of the protein]')
			except:
				raise Exception('pdb ID should be given or pass a string of pdb file with the argument: file=[the given string]')
		try:
			downloadpath=kargs['path']
		except:
			pass
		else:
			fso = open(downloadpath, 'w')
			fso.write(f)
			fso.close()
		lines=f.split("\n")
		self.models=[]
		self.helixes=[]
		self.sheets=[]
		modelendlist=[]
		#==================Atoms Building======================
		#preprocess: chunks incision
		#models incision
		for i in range(len(lines)):
			if lines[i].find('HELIX')==0:
				self.helixes.append(Helix(lines[i]))
			elif lines[i].find('SHEET')==0:
				self.sheets.append(lines[i])
			elif lines[i].find('MODEL')==0:
				self.models.append(i)
			elif lines[i].find('ENDMDL')==0:
				modelendlist.append(i)
		if len(self.models)==0:
			self.models.append(lines)
		else:
			for i in range(len(self.models)):
				self.models[i]=lines[self.models[i]+1:modelendlist[i]]
		#chains incision
		for i in range(len(self.models)):
			temlist=[]
			for j in range(len(self.models[i])):
				if self.models[i][j].find('TER')==0:
					temlist.append(j)
			temnum=len(temlist)-1
			for k in range(temnum+1):
				j=temnum-k
				if j==0:
					temlist[j]=self.models[i][0:temlist[j]]
				else:
					temlist[j]=self.models[i][temlist[j-1]:temlist[j]]
			self.models[i]=temlist
		#residues incision
		for i in self.models:
			for j in range(len(i)):
				chaintem=[]
				conuter=0
				residuestem=[]
				tematomlist=[]
				for k in i[j]:
					if k.find('ATOM')==0:
						tematomlist.append(k)
				i[j]=tematomlist
				for k in range(len(i[j])):
					if k==0:
						counter=i[j][k][22:26]
						residuestem.append(i[j][k])
					else:
						if i[j][k][22:26]==counter:
							residuestem.append(i[j][k])
						else:
							chaintem.append((counter,residuestem[:]))
							counter=i[j][k][22:26]
							residuestem=[]
							residuestem.append(i[j][k])
				i[j]=chaintem
		#Build!
		for m in range(len(self.models)):
			clist=[]
			for c in self.models[m]:
				cname=c[0][1][0][21]
				rlist=[]
				for r in c:
					rname=r[1][0][17:20].replace(' ','')
					rnum=int(r[1][0][22:26])
					for a in range(len(r[1])):
						r[1][a]=Atom(r[1][a][12:16].replace(' ',''),(float(r[1][a][30:38]),float(r[1][a][38:46]),float(r[1][a][46:54])),r[1][a][76:78].replace(' ',''))
					rlist.append(Residue(rname,rnum,r[1]))
				clist.append(Chain(cname,rlist))
			self.models[m]=Model(clist)
		#=======================Secondary Structure Building=================
		#preprocess: chunks incision
		#sheet groups incision
		temlastlist=[]
		temsheet=''
		if len(self.sheets)!=0:
			for i in range(len(self.sheets)):
				if i==0:
					temsheet=self.sheets[i][11:14]
					temlastlist.append(i)
				else:
					if temsheet!=self.sheets[i][11:14]:
						temlastlist.append(i)
						temsheet=self.sheets[i][11:14]
			temnum=len(temlastlist)-1
			for i in range(temnum+1):
				if i==temnum:
					temlastlist[i]=self.sheets[temlastlist[i]:]
				else:
					temlastlist[i]=self.sheets[temlastlist[i]:temlastlist[i+1]]
			self.sheets=temlastlist
			for i in range(temnum+1):
				self.sheets[i]=Sheet(self,self.sheets[i])


	def __getitem__(self, key):
		if isinstance(key,int):
			return self.models[key]
		elif isinstance(key,slice):
			return self.models[key]
		else:
			raise KeyError("Expecting an intege but a "+str(type(key))+"is given")
	def __setitem__(self, key, value):
		if isinstance(key,int):
			self.models[key]=value
		elif isinstance(key,slice):
			self.models[key]=value
		else:
			raise KeyError("Expecting an intege but a "+str(type(key))+"is given")

	def __str__(self):
		return self.name

	def __len__(self):
		return len(self.models)
	def seqlen(self):
		return self[0].seqlen()
	def __iter__(self):
		return iter(self.models)



	def copy(self,newname=None, selectmodel=None, selectchian=None, selectresidue=None):
		tem=copy.deepcopy(self)
		if newname!=None:
			tem.name=newname
		if selectmodel!=None:
			try:
				temobj=selectmodel[0]
			except:
				tem.models=[tem[selectmodel]]
			else:
				temlist=[]
				for i in selectmodel:
					temlist.append(tem[i])
				tem.models=temlist
		if selectchian!=None:
			print('[warning]Secondary structures will not update')
			try:
				temobj=selectchian[0]
			except:
				for m in tem:
					m.chains=[m[selectchian]]
			else:
				for m in tem:
					temlist=[]
					for i in selectchian:
						temlist.append(m[i])
					m.chains=temlist
		if selectresidue!=None:
			print('[warning]Secondary structures will not update')
			try:
				temobj=selectresidue[0]
			except:
				for m in tem:
					for c in m:
						c.residues=[c[selectresidue]]
			else:
				for m in tem:
					for c in m:
						temlist=[]
						for i in selectresidue:
							temlist.append(c[i])
						c.residues=temlist
		return tem
	def exportpdb(self,path=None):
		if path==None:
			if self.name==None:
				exist=True
				index=1
				path=os.getcwd()+os.sep+'export'+str(index)+'.pdb'
				while os.path.isfile(path)==True:
					index+=1
					path=os.getcwd()+os.sep+'export'+str(index)+'.pdb'
			else:
				path=os.getcwd()+os.sep+self.name+'.pdb'
		with open(path,'a') as f:
			if len(self)==1:
				counter=1
				for c in self[0]:
					for r in c:
						for a in r:
							if len(a.name)<=3:
								f.write('ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f                      %2s  \n'% (counter, a.name, r.name, c.name, r.num, a.pos[0], a.pos[1], a.pos[2], a.ele))
							else:
								f.write('ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f                      %2s  \n'% (counter, a.name, r.name, c.name, r.num, a.pos[0], a.pos[1], a.pos[2], a.ele))
							counter+=1
					f.write('TER   %5d      %3s %1s% 4d                                                      \n'% (counter, r.name, c.name, r.num))
				f.write('END                                                                             ')
			else:
				for m in self:
					f.write('MODEL     % 4d                                                                  \n')
					counter=1
					for c in m:
						for r in c:
							for a in r:
								if len(a.name)<=3:
									f.write('ATOM  %5d  %-3s %3s %1s% 4d    %8.3f%8.3f%8.3f                      %2s  \n'% (counter, a.name, r.name, c.name, r.num, a.pos[0], a.pos[1], a.pos[2], a.ele))
								else:
									f.write('ATOM  %5d %4s %3s %1s% 4d    %8.3f%8.3f%8.3f                      %2s  \n'% (counter, a.name, r.name, c.name, r.num, a.pos[0], a.pos[1], a.pos[2], a.ele))
								counter+=1
						f.write('TER   %5d      %3s %1s%4d                                                      \n'% (counter, r.name, c.name, r.num))
					f.write('ENDMDL                                                                          \n')
				f.write('END                                                                             \n')

	#Two ends of the insert will align to the residues 'beg' and 'end' and then be removed.
	def couple(self,chain,beg,end,insert, name=None, score=False, debug=False):
		if name==None:
			tem=self.copy(newname=self.name+'ins',selectmodel=0)
		else:
			tem=self.copy(newname=name,selectmodel=0)
		newchain=tem[0][chain]
		begr=newchain[beg]
		endr=newchain[end]
		listA=newchain[:beg+1]
		listB=copy.deepcopy(insert)
		listC=newchain[end:]
		#renumber
		counter=beg
		for i in listB:
			i.num=counter
			counter+=1
		counter=counter-1-end
		for i in listC:
			i.num=i.num+counter


		#align insert to target
		p1=begr['CA'].pos
		p2=begr['C'].pos
		p3=endr['CA'].pos
		p4=endr['N'].pos

		q1=listB[0]['CA'].pos
		q2=listB[0]['C'].pos
		q3=listB[-1]['CA'].pos
		q4=listB[-1]['N'].pos

		ccp=(p1+p2+p3+p4)/4
		ccq=(q1+q2+q3+q4)/4
		cp24=(p2+p4)/2
		cq24=(q2+q4)/2
		mat=(cq24-ccq).rot_matrix(cp24-ccp)
		for i in listB:
			for j in i:
				j.pos=mat*(j.pos-ccq)+ccp
		q1=listB[0]['CA'].pos
		q2=listB[0]['C'].pos
		q3=listB[-1]['CA'].pos
		q4=listB[-1]['N'].pos
		ccq=(q1+q2+q3+q4)/4
		cq24=(q2+q4)/2
		cp34=(p3+p4)/2
		cq34=(q3+q4)/2
		cqplane=cdn.Plane((cq24-ccq,ccq))
		projcp34=cp34+cqplane
		projcq34=cq34+cqplane
		vp=projcp34-ccp
		vq=projcq34-ccq
		cosa=vp*vq/abs(vp)/abs(vq)
		rotafunc=cdn.rotation(ccq,ccq-cq24,cosa)

		for i in listB:
			for j in i:
				j.pos=rotafunc(j.pos)

		q1=listB[0]['CA'].pos
		q2=listB[0]['C'].pos
		q3=listB[-1]['CA'].pos
		q4=listB[-1]['N'].pos

		ascore=((p1-q1)**2)+((p2-q2)**2)+((p3-q3)**2)+((p4-q4)**2)
		if debug==True:
			print(ascore)
			#plot.pointviewer([[p1,q1],[p2,q2],[p3,q3],[p4,q4]])

		listBcopy=copy.deepcopy(listB)
		temascore=1000000
		temascorepre=10000000
		while temascore<temascorepre:
			q1=listBcopy[0]['CA'].pos
			q2=listBcopy[0]['C'].pos
			q3=listBcopy[-1]['CA'].pos
			q4=listBcopy[-1]['N'].pos
			ccq=(q1+q2+q3+q4)/4


			d1=(p1-q1)
			d3=(p3-q3)
			d=(d1-d3)/6
			t=q1+d
			m=(q1-ccq).rot_matrix(t-ccp)
			for i in listBcopy:
				for j in i:
					j.pos=(m*(j.pos-ccq))+ccp

			q1=listBcopy[0]['CA'].pos
			q2=listBcopy[0]['C'].pos
			q3=listBcopy[-1]['CA'].pos
			q4=listBcopy[-1]['N'].pos
			ccq=(q1+q2+q3+q4)/4

			d2=(p2-q2)
			d4=(p4-q4)
			d=(d2-d4)/6
			t=q2+d
			m=(q2-ccq).rot_matrix(t-ccp)
			for i in listBcopy:
				for j in i:
					j.pos=(m*(j.pos-ccq))+ccp


			
			temascorepre=temascore
			temascore=((p1-q1)**2)+((p2-q2)**2)+((p3-q3)**2)+((p4-q4)**2)

			if debug==True:
				print(temascore)
				#plot.pointviewer([[p1,q1],[p2,q2],[p3,q3],[p4,q4]])

		if temascore<ascore:
			ascore=temascore
			listB=listBcopy


		#build a new chain
		if debug==True:
			debugP=copy.deepcopy(tem)
			temlistA=listA[:]
			listA.extend(listC)
			newchain.residues=listA
			debugP[0][chain]=newchain
			debugP[0].chains.append(Chain('W',listB[:]))
			debugP.exportpdb(os.getcwd()+os.sep+debugP.name+'_debug.pdb')

		listB=listB[1:-1]
		listA.extend(listB)
		listA.extend(listC)
		newchain.residues=listA
		tem[0][chain]=newchain
		
		if score==True:
			with open(os.getcwd()+os.sep+'couple_scoring','a') as f:
				f.write(tem.name+' '+str(ascore)+'\n')
		elif isinstance(score,str):
			with open(score,'a') as f:
				f.write(tem.name+' '+str(ascore))
		return tem









#========================================================================================================================
#=================================================A Protein with screen tasks============================================
#========================================================================================================================


#A class for screening
class ProteinS(Protein):
	def __init__(self,pdbid=None,**kargs):
		debug=kargs.get('debug',False)
		if debug==True:
			print('#### Debug Mode Start')
		if pdbid!=None:#It's a pdb ID, download it from the website.
			url='https://files.rcsb.org/download/'+pdbid+'.pdb'
			with request.urlopen(url) as f:
				if f.status== 200:
					f=f.read().decode()
				else:
					raise Exception("Failed to Download"+pdbid+'.pdb. Error '+f.status+' '+f.reason)
			try:
				self.name=kargs['name']
			except:
				self.name=pdbid
			if debug==True:
				print('#### Start processing '+pdbid+' - '+self.name)
		else:
			try:
				f=kargs['file']
				try:
					self.name=kargs['name']
				except:
					self.name=None
					print('[Warning]Can not find the name of this protein, which can given by the argument: name=[the name of the protein]')
			except:
				raise Exception('pdb ID should be given or pass a string of pdb file with the argument: file=[the given string]')
			if debug==True:
				print('#### Start processing '+self.name)
		self.J=True
		self.lines=f.split("\n")
		self.models=[]
		self.helixes=[]
		self.sheets=[]
		modelendlist=[]
		self.score=None
		#===================Preprocess Conditions=====================
		if debug==True:
			print('#### Preprocess Conditions')
		try:
			conditions=kargs['precon']
		except:
			pass
		else:
			try:
				temcon=conditions[0]
			except:
				self.J=conditions(self,*(kargs.get('preconargs',())),**(kargs.get('preconkargs',{})))
			else:
				for i in range(len(conditions)):
					if conditions[i](self,*(kargs.get('preconargs',())[i]),**(kargs.get('preconkargs',{})[i]))==False:
						self.J=False
						break
		#===================Preprocess: Chunks Incision=====================
		#models incision
		if self.J==True:
			if debug==True:
				print('#### Preprocessing')
			for i in range(len(self.lines)):
				if self.lines[i].find('HELIX')==0:
					self.helixes.append(Helix(self.lines[i]))
				elif self.lines[i].find('SHEET')==0:
					self.sheets.append(self.lines[i])
				elif self.lines[i].find('MODEL')==0:
					self.models.append(i)
				elif self.lines[i].find('ENDMDL')==0:
					modelendlist.append(i)
				else:
					try:
						conditions=kargs['linetrue']
					except:
						pass
					else:
						try:
							temcon=conditions[0]#A list or not?
						except:
							self.J=(self.lines[i].find(conditions)!=-1)
						else:#a list
							for j in conditions:
								if self.lines[i].find(j)==-1:
									self.J=False
									break
					if self.J==False:
						break

					try:
						conditions=kargs['linefalse']
					except:
						pass
					else:
						try:
							temcon=conditions[0]#A list or not?
						except:
							self.J=(self.lines[i].find(conditions)==-1)
						else:#a list
							for j in conditions:
								if self.lines[i].find(j)!=-1:
									self.J=False
									break
					if self.J==False:
						break


					try:
						conditions=kargs['linecon']
					except:
						pass
					else:
						try:
							temcon=conditions[0]#A list or not?
						except:#just one function
							self.J=conditions(self,self.lines[i],*(kargs.get('lineconargs',())),**(kargs.get('lineconkargs',{})))
						else:#a list
							for j in range(len(conditions)):
								if conditions[j](self,self.lines[i],*(kargs.get('lineconargs',())[j]),**(kargs.get('lineconkargs',{})[j]))==False:
									self.J=False
									break
					if self.J==False:
						break

			if self.J==True:
				if len(self.models)==0:
					self.models.append(self.lines)
				else:
					for i in range(len(self.models)):
						self.models[i]=self.lines[self.models[i]+1:modelendlist[i]]
				del self.lines
				#chains incision
				for i in range(len(self.models)):
					temlist=[]
					for j in range(len(self.models[i])):
						if self.models[i][j].find('TER')==0:
							temlist.append(j)
					temnum=len(temlist)-1
					for k in range(temnum+1):
						j=temnum-k
						if j==0:
							temlist[j]=self.models[i][0:temlist[j]]
						else:
							temlist[j]=self.models[i][temlist[j-1]:temlist[j]]
					self.models[i]=temlist
				#residues incision
				for i in self.models:
					for j in range(len(i)):
						chaintem=[]
						conuter=0
						residuestem=[]
						tematomlist=[]
						for k in i[j]:
							if k.find('ATOM')==0:
								tematomlist.append(k)
						i[j]=tematomlist
						for k in range(len(i[j])):
							if k==0:
								counter=i[j][k][22:26]
								residuestem.append(i[j][k])
							else:
								if i[j][k][22:26]==counter:
									residuestem.append(i[j][k])
								else:
									chaintem.append((counter,residuestem[:]))
									counter=i[j][k][22:26]
									residuestem=[]
									residuestem.append(i[j][k])
						i[j]=chaintem


				#==================Atoms Building Conditions======================
				if debug==True:
					print('#### Atoms Building Conditions')
				try:
					conditions=kargs['atmcon']
				except:
					pass
				else:
					try:
						temcon=conditions[0]
					except:
						self.J=conditions(self,*(kargs.get('atmconargs',())),**(kargs.get('atmconkargs',{})))
					else:
						for i in range(len(conditions)):
							if conditions[i](self,*(kargs.get('atmconargs',())[i]),**(kargs.get('atmconkargs',{})[i]))==False:
								self.J=False
								break
				#==================Atoms Building======================

				if self.J==True:
					if debug==True:
						print('#### Building Atoms')
					for m in range(len(self.models)):
						clist=[]
						for c in self.models[m]:
							cname=c[0][1][0][21]
							rlist=[]
							for r in c:
								rname=r[1][0][17:20].replace(' ','')
								rnum=int(r[1][0][22:26])
								for a in range(len(r[1])):
									r[1][a]=Atom(r[1][a][12:16].replace(' ',''),(float(r[1][a][30:38]),float(r[1][a][38:46]),float(r[1][a][46:54])),r[1][a][76:78].replace(' ',''))
								rlist.append(Residue(rname,rnum,r[1]))
							clist.append(Chain(cname,rlist))
						self.models[m]=Model(clist)
					#=======================Secondary Structure Building=================
					if debug==True:
						print('#### Building Secondary Structures')
					#preprocess: chunks incision
					#sheet groups incision
					temlastlist=[]
					temsheet=''
					if len(self.sheets)!=0:
						for i in range(len(self.sheets)):
							if i==0:
								temsheet=self.sheets[i][11:14]
								temlastlist.append(i)
							else:
								if temsheet!=self.sheets[i][11:14]:
									temlastlist.append(i)
									temsheet=self.sheets[i][11:14]
						temnum=len(temlastlist)-1
						for i in range(temnum+1):
							if i==temnum:
								temlastlist[i]=self.sheets[temlastlist[i]:]
							else:
								temlastlist[i]=self.sheets[temlastlist[i]:temlastlist[i+1]]
						self.sheets=temlastlist
						for i in range(temnum+1):
							self.sheets[i]=Sheet(self,self.sheets[i])



					
					#==================Scoring======================	
					if debug==True:
						print('#### Scoring')
					try:
						scoring=kargs['score']
					except:
						pass
					else:
						try:
							temcon=scoring[0]
						except:
							self.J=scoring(self,*(kargs.get('scoreargs',())),**(kargs.get('scorekargs',{})))
						else:
							for i in range(len(scoring)):
								if scoring[i](self,*(kargs.get('scoreargs',())[i]),**(kargs.get('scorekargs',{})[i]))==False:
									self.J=False
									break

					#==================Download======================	
					if self.J==True:
						if debug==True:
							print('#### Downloading')
						try:
							downloadpath=kargs['path']
						except:
							pass
						else:
							fso = open(downloadpath, 'w')
							fso.write(f)
							fso.close()
					else:
						if debug==True:
							print('#### Conditions are not met, the downloading is aborted')
				else:
					if debug==True:
						print('#### Conditions are not met, the atoms building is aborted')
			else:
				if debug==True:
					print('#### Conditions are not met, the preprocessing is breaked off')
		else:
			if debug==True:
				print('#### Conditions are not met, the preprocessing is aborted')







#========================================================================================================================
#=====================================================||||||||||||||=====================================================
#=====================================================||==========||=====================================================
#=====================================================||=================================================================
#=====================================================||||||||||||||=====================================================
#=================================================================||=====================================================
#=====================================================||==========||=====================================================
#=====================================================||||||||||||||=====================================================
#========================================================================================================================
#===================================================Classes definition===================================================



#========================================================================================================================
#=========================================================||||===========================================================
#=========================================================||||===========================================================
#=========================================================||=============================================================
#========================================================================================================================
#================================================Preprocessing Conditions================================================






#==[Warning] Preprocessing Conditions below might be slow, use arguments 'linetrue' and 'linefalse' in ProteinS instead.=
#==[Warning] Preprocessing Conditions below might be slow, use arguments 'linetrue' and 'linefalse' in ProteinS instead.=
#==[Warning] Preprocessing Conditions below might be slow, use arguments 'linetrue' and 'linefalse' in ProteinS instead.=
#Return whether the protein is expressed by the specific systems, e.g. 'ESCHERICHIA COLI' or ['ESCHERICHIA COLI','HOMO SAPIENS']
def screxp(protein,expressionsystem):
	print("[Warning] function 'screxp' might be slow, use argument 'linetrue' in ProteinS instead.")
	if type(expressionsystem)==list:
		for i in range(len(protein.lines)):
			for j in expressionsystem:
				if protein.lines[i].find('EXPRESSION_SYSTEM: '+j+';')!=-1:
					return True
		return False
	else:
		for i in range(len(protein.lines)):
			if protein.lines[i].find('EXPRESSION_SYSTEM: '+expressionsystem+';')!=-1:
				return True
		return False


#Return whether the protein is not expressed by the specific systems, e.g. 'ESCHERICHIA COLI' or ['ESCHERICHIA COLI','HOMO SAPIENS']
def scrnotexp(protein,expressionsystem):
	print("[Warning] function 'scrnotexp' might be slow, use argument 'linefalse' in ProteinS instead.")
	if type(expressionsystem)==list:
		for i in range(len(protein.lines)):
			for j in expressionsystem:
				if protein.lines[i].find('EXPRESSION_SYSTEM: '+j+';')!=-1:
					return False
		return True
	else:
		for i in range(len(protein.lines)):
			if protein.lines[i].find('EXPRESSION_SYSTEM: '+expressionsystem+';')!=-1:
				return False
		return True


#Return whether the protein comes from the specific organism, e.g. 'ESCHERICHIA COLI' or ['ESCHERICHIA COLI','HOMO SAPIENS']
def scrorg(protein,expressionsystem):
	print("[Warning] function 'scrorg' might be slow, use argument 'linetrue' in ProteinS instead.")
	if type(expressionsystem)==list:
		for i in range(len(protein.lines)):
			for j in expressionsystem:
				if protein.lines[i].find('ORGANISM_SCIENTIFIC: '+j+';')!=-1:
					return True
		return False
	else:
		for i in range(len(protein.lines)):
			if protein.lines[i].find('ORGANISM_SCIENTIFIC: '+expressionsystem+';')!=-1:
				return True
		return False


#Return whether the protein comes from organisms except the specific organism, e.g. 'ESCHERICHIA COLI' or ['ESCHERICHIA COLI','PYROCOCCUS HORIKOSHII']
def scrnotorg(protein,expressionsystem):
	print("[Warning] function 'scrnotorg' might be slow, use argument 'linefalse' in ProteinS instead.")
	if type(expressionsystem)==list:
		for i in range(len(protein.lines)):
			for j in expressionsystem:
				if protein.lines[i].find('ORGANISM_SCIENTIFIC: '+j+';')!=-1:
					return False
		return True
	else:
		for i in range(len(protein.lines)):
			if protein.lines[i].find('ORGANISM_SCIENTIFIC: '+expressionsystem+';')!=-1:
				return False
		return True


#========================================================================================================================
#=========================================================||=============================================================
#=========================================================||=============================================================
#=========================================================||||||=========================================================
#========================================================================================================================
#================================================Lines Reading Conditions================================================

# Record the ORGANISM_SCIENTIFIC information
def recorg(protein, line):
	if line.find('ORGANISM_SCIENTIFIC: ')!=-1:
		protein.organism= line[line.find('ORGANISM_SCIENTIFIC: ')+21:line.find(';')]
	return True

# Record the EXPRESSION_SYSTEM information
def recexp(protein, line):
	if line.find('EXPRESSION_SYSTEM: ')!=-1:
		protein.expression= line[line.find('EXPRESSION_SYSTEM: ')+19:line.find(';')]	
	return True



#========================================================================================================================
#=========================================================||||||=========================================================
#=========================================================||||||=========================================================
#=========================================================||||||=========================================================
#========================================================================================================================
#=================================================Atoms Building Conditions==============================================



#Return whether the protein have a sheet/ sheets
def scrsheet(protein):
	if len(protein.sheets)==0:
		return False
	else:
		return True


#Return whether the protein is without any sheet
def scrnosheet(protein):
	if len(protein.sheets)==0:
		return True
	else:
		return False

#Return whether the protein have a helix/ helixes
def scrhelix(protein):
	if len(protein.helixes)==0:
		return False
	else:
		return True


#Return whether the protein is without any helix
def scrnohelix(protein):
	if len(protein.helixes)==0:
		return True
	else:
		return False


#========================================================================================================================
#===========================================================||||=========================================================
#===========================================================||===========================================================
#=========================================================||||===========================================================
#========================================================================================================================
#========================================================Scoring=========================================================

#score acorrding to the full length of the protein 
def scrfulllength(protein, scorefunction, threshold=None):
	protein.flscore=scorefunction(protein.seqlen())
	if threshold==None:
		return True
	else:
		if protein.flscore>threshold:
			return True
		else:
			return False

#score acorrding to the ratio of loops length to full length of the protein 
def scrloopratio(protein, scorefunction=lambda x:1-x, threshold=None):
	lenh=0
	for i in protein.helixes:
		lenh+=len(i)
	lens=0
	for i in protein.sheets:
		for j in i:
			lens+=len(j)
	lenf=protein.seqlen()
	lenl=lenf-lenh-lens
	protein.lrscore=scorefunction(lenl/lenf)
	if threshold==None:
		return True
	else:
		if protein.lrscore>threshold:
			return True
		else:
			return False

#score acorrding to the full length of the protein 
def scrsheetratio(protein, scorefunction=lambda x:x, threshold=None):
	lenh=0
	for i in protein.helixes:
		lenh+=len(i)
	lens=0
	for i in protein.sheets:
		for j in i:
			lens+=len(j)
	protein.srscore=scorefunction(lens/(lens+lenh))
	if threshold==None:
		return True
	else:
		if protein.srscore>threshold:
			return True
		else:
			return False


#score acorrding to the full length of the protein 
def scrhelixratio(protein, scorefunction=lambda x:x, threshold=None):
	lenh=0
	for i in protein.helixes:
		lenh+=len(i)
	lens=0
	for i in protein.sheets:
		for j in i:
			lens+=len(j)
	protein.hrscore=scorefunction(lenh/(lens+lenh))
	if threshold==None:
		return True
	else:
		if protein.hrscore>threshold:
			return True
		else:
			return False

#score acorrding to a position function. 
def scrposition(protein,positionfunction, scorefunction, threshold=None):
	protein.pscore=scorefunction(positionfunction(protein))
	if threshold==None:
		return True
	else:
		if protein.pscore>threshold:
			return True
		else:
			return False