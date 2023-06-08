#!/usr/bin/python

import sys
import os
import re
import gzip
import numpy as np
from collections import defaultdict 

class BlastReport():

	class __BlastQuery():
		
		def __init__(self,query_name,query_length,query_hits):
			self.name = query_name
			self.length = query_length
			self.__query_hits = query_hits						
			
		def hits(self):
			return self.__query_hits

	class __SingleHit(): #__init__(self,subject_ref_data,hit_data): ## class for a signle hit
		
		def __init__(self,subject_ref_name,subject_ref_length,hit_data):
						
			self.reference_name = subject_ref_name
			self.reference_length = subject_ref_length
			self.__hit_data = hit_data
			self.__parse()

		def __parse(self):
			
			self.score = 0
			self.evalue = np.inf
			self.identities = (0,0,0)
			self.gaps = (0,0,0)
			self.strand = './.'			
			self.start_pos = 0
			self.end_pos = 0
			self.seq = ''	
			self.alignment_report = ''	
			
			if self.__hit_data == []:
				return		
			
			
			self.score = eval(re.findall('[0-9]+ bits',self.__hit_data[0])[0].replace(' bits',''))
			self.evalue = eval(self.__hit_data[0].split('Expect = ')[1])
			
			identities_tmp,gaps_tmp = self.__hit_data[1].split(', ')
			self.identities = tuple(eval_list(re.findall('[0-9]+',identities_tmp)))
			self.gaps = tuple(eval_list(re.findall('[0-9]+',gaps_tmp)))
			
			self.strands = self.__hit_data[2].replace('Strand=','').replace('Plus','+').replace('Minus','-')
			
			#--alignment parsing----
			
			blast_aln = re.compile('[-,A,T,G,C,a,t,g,c,0-9]+')
			alignment_data = self.__hit_data[3:]			
			alignment_report_tmp = []						
			query_seq = ''
			query_pos = []
			subject_seq = ''
			subject_pos = []
			
			for i in range(len(alignment_data)):
				line = alignment_data[i]
				if re.search('Lambda',line):
					break
				alignment_report_tmp.append(line)
				if i%3 == 2:
					alignment_report_tmp.append('\n')			

				if re.search('Query',line): # query line
					q_pos1,q_seq,q_pos2 = eval_list(re.findall(blast_aln,line.replace('Query','').strip()))
					query_seq += q_seq
					query_pos += [q_pos1,q_pos2]
				elif re.search('Sbjct',line): # subject line
					s_pos1,s_seq,s_pos2 = eval_list(re.findall(blast_aln,line.replace('Sbjct','').strip()))
					subject_seq += s_seq
					subject_pos += [s_pos1,s_pos2]
				else:
					continue
					
			self.start_pos = subject_pos[0]
			self.end_pos = subject_pos[-1]
			self.seq = subject_seq.replace('-','')	
			self.alignment_report = '\n'.join(alignment_report_tmp).strip()

	
	def __init__(self,report_file):
		self.__fetch_out = []
		self.__report_file = report_file
		self.__load_report()
		

	def fetch(self):
		return self.__fetch_out
		

	def __load_report(self):
		
		
		## hit set save : 1. >chr1 = new subject / 2. Score = sample subject, new alignment , 3. Lambda = end alignment
		
		## flag ; 0 = off / 1 = on
				
		query_flag = 0
		hit_flag = 0
		subject_ref_flag = 0		
		query_name =  ''
		subject_ref_name = ''
		query_length = 0
		nohit_flag = 1

		hit_set = []
		
		finp = open(self.__report_file)
		for line in finp:
			line = line.rstrip()
			if line == '':
				continue
			
			if subject_ref_flag == 1:
				subject_ref_length = eval(re.findall('[0-9]+',line)[0])
				subject_ref_flag = 0			
			#----------------
			
			if query_flag == 1:
				query_length = eval(re.findall('[0-9]+',line)[0])
				query_flag = 0			
			#----------------			
			if re.search('Query=',line): ## query start
				query_name = line.strip().split('= ')[1]
				query_hits = []
				query_flag = 1
				nohit_flag = 1
							
			elif re.search('^>',line): # hit_start
				if not hit_set == []:				
					query_hits.append(self.__SingleHit(subject_ref_name,subject_ref_length,hit_set))
					subject_ref_flag,hit_flag,hit_set = (0,0,[]) # flush
				
				subject_ref_name = line.strip().replace('>','')
				subject_ref_flag = 1
				hit_flag = 0
				
			elif re.search('Score = ',line): # hit start & previous hit end

				if not hit_set == []:				
					query_hits.append(self.__SingleHit(subject_ref_name,subject_ref_length,hit_set))
					hit_set  = [] # hit_set flush
													
				#-----hit save ---
				hit_set.append(line)
				nohit_flag = 0
				hit_flag = 1									
								
			elif re.search('Effective search space used',line): ## hit all end
				
				#----hit save----
				if nohit_flag == 1: # no hit found
					query_hits.append(self.__SingleHit('*',0,hit_set))				
				else:				
					
					if not hit_set == []:
						query_hits.append(self.__SingleHit(subject_ref_name,subject_ref_length,hit_set))
						subject_ref_flag,hit_flag,hit_set = (0,0,[]) # flush
				self.__fetch_out.append(self.__BlastQuery(query_name,query_length,query_hits))				
												
			else:								
				if hit_flag == 1:
					hit_set.append(line)															
						
			#-------------------------			
			
		finp.close()

if __name__=="__main__":
	
	header = '''
< Class Usages>

1. BlastReport
blast = BlastReport(report_file = 'blast.out')

for query in blast.fetch():

	query.name
	query.length
		
	for hit in query.hits():
	
		hit.reference_name
		hit.reference_length
		hit.score
		hit.evalue
		hit.identities
		hit.gaps
		hit.strand
		hit.start_pos
		hit.end_pos
		hit.seq
		hit.alignment_report
		
	'''	
	print header;sys.exit(2)
