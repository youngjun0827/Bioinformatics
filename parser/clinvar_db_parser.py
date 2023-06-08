#!/data/youngjun/local/bin/python	
#-*-coding:utf-8-*-

from unidecode import unidecode
import sys,os,pickle,string,unicodedata,gzip,time,getopt,datetime,re,pymysql
import pyfaidx
import numpy as np
import xml.etree.cElementTree as ET
reload(sys)
sys.setdefaultencoding('utf-8')

#######################################
### Manual patch : Go to line 553   ###
#######################################

def usage(errmsg):
	title		=	'''
#####################################################
## Clinvar DB Parser v2.0 (by youngjun 2018.12.20) ##
#####################################################
'''
	sys.stderr.write(title)
	if not errmsg == '':
		sys.stderr.write('\n'+errmsg+'\n')
	errtext	=	'''
usage: %s [ -T input type ( vcf / xml ) ] [ Options ]

<Input type & options>

 <vcf>
	-R [ reference fasta ] (required)
	-v [ clinvar vcf file ] (required)
	-o [ output npz name ] (required)  
	
 <xml>
	-x [ clinvar xml file ] (required)
	-o [ output npz name ] (required)
'''%sys.argv[0]
	sys.stderr.write(errtext+'\n')

def required_option_check(required_options,parsed_options):
	for required_opt in required_options:
		if not required_opt in parsed_options:
			errmsg                    =       '#ERROR: option %s is required'%required_opt
			usage(errmsg);sys.exit(2)
			
def option_file_check(parsed_options,parsed_args,input_file_list):
	for check in input_file_list:
		if check	==	'Null':
			continue
		if not os.path.isfile(check):
			errmsg	=	'#ERROR: file %s was not found'%check
			usage(errmsg);sys.exit(2)

	check_list  =       list(set(parsed_args))
	for check in check_list:
		if parsed_args.count(check)>1:
			errmsg    =       '#ERROR:  same argument was found ( filename : %s )'%check
			usage(errmsg);sys.exit(2)

def load_gene_list(gene_list_file):
	gene_list = list()
	finp = open(gene_list_file)
	for line in finp:
		gene = line.strip()
		gene_list.append(gene)
	finp.close()	
	return gene_list

def hgvs_p_check(hgvs_p):	
	hgvs_p = hgvs_p.replace('Ter','*')
	if not hgvs_p.startswith('p.'):
		return hgvs_p
	pos_set = re.findall('[0-9]+',hgvs_p)
	aa_pos	= int(pos_set[0])
	if not len(pos_set) == 1:
		return hgvs_p
	tToken	=	hgvs_p.split('.')[1].split(str(aa_pos))
	if tToken[0]==tToken[1]:
		converted_hgvs_p	='p.%s%s='%(tToken[0],aa_pos)
	else:
		f_aaset = re.findall('[A-Z]{1}[a-z]{2}',tToken[0])
		r_aaset = re.findall('[A-Z]{1}[a-z]{2}',tToken[1])
		if not len(f_aaset) == len(r_aaset):
			return hgvs_p
		conv_f_aaset = f_aaset
		conv_r_aaset = r_aaset
		for i in range(len(f_aaset)):
			f_aa = f_aaset[i]
			r_aa = r_aaset[i]
			if f_aa == r_aa:
				aa_pos += 1
				conv_f_aaset = conv_f_aaset[1:]
				conv_r_aaset = conv_r_aaset[1:]
			else:
				break
		converted_hgvs_p = 'p.%s%s%s'%(''.join(conv_f_aaset),aa_pos,''.join(conv_r_aaset))				
	return converted_hgvs_p
		
def load_transcript_list():
	transcript = dict()
	try:
		finp = gzip.open('/data/DB_Analy/Annotation/refFlat.txt.gz')
	except:
		sys.stderr.write('file "/data/DB_Anlay/Annotation/refFlat.txt.gz" was not found!');sys.exit(2)
	for line in finp:
		token = line.strip().split('\t')
		gene_name,transcript_id = token[:2]
		if not re.search('NM\_[0-9]+',transcript_id):
			continue
		transcript[transcript_id] = gene_name
	finp.close()
	return transcript

def parse_var_set(ref_assertion):
	condition_set = list()	
	def parse_var_summary(var_summary):
		try:
			nm_raw_id = re.findall('NM\_[0-9]+\.[0-9]+',var_summary)[0]
			nm_id = nm_raw_id.split('.')[0]
		except:
			nm_raw_id,nm_id = '-','-'
		if re.search('c\.',var_summary):
			hgvs_c = 'c.%s'%(var_summary.split('c.')[1].split('(')[0].strip())
		elif re.search('r\.',var_summary):
			hgvs_c = 'r.%s'%(var_summary.split('r.')[1].split('(')[0].strip())
		else:
			hgvs_c = '-'			
		if re.search('p\.',var_summary):
			hgvs_p = 'p.%s'%(hgvs_p_check(var_summary.split('p.')[1].replace(')','').strip()))
		else:
			hgvs_p = '-'
		return nm_id,hgvs_c,hgvs_p
		
	genotype_set = ref_assertion.find('GenotypeSet')
	try:
		measure_set_group = genotype_set.findall('MeasureSet')
		var_type = genotype_set.attrib['Type']
		var_id = int(genotype_set.attrib['ID'].replace('"',''))
		var_name = genotype_set.find('Name').findtext('ElementValue') 
	except:
		measure_set_group = ref_assertion.findall('MeasureSet')
		var_type = measure_set_group[0].attrib['Type']
		var_id = int(measure_set_group[0].attrib['ID'].replace('"',''))
		var_name = measure_set_group[0].find('Name').findtext('ElementValue')	
	var_set = dict()		
#-----NMid / hgvs c / hgvs p parsing
	for measure_set in measure_set_group:
 		for measure in measure_set.findall('Measure'):
 			allele_type = measure.attrib['Type']
 			allele_id = int(measure.attrib['ID'].replace('"',''))
 			values = measure.findall('Name')
 			for value in values:
 				for ele_value in value.findall('ElementValue'):
 					if not ele_value.attrib['Type'] == 'Preferred':
 						continue
 					var_summary = ele_value.text
 			try:
 		 		nm_id,hgvs_c,hgvs_p = parse_var_summary(var_summary)
 		 	except:
 		 		nm_id,hgvs_c,hgvs_p = '-','-','-' 					
 		 	var_set[allele_id] = (nm_id,hgvs_c,hgvs_p) 			
 		for measure_name in measure_set.findall('Name'):
 			for ele_name in measure_name.findall('ElementValue'):
 				if not ele_name.attrib['Type'] == 'Preferred':
 					continue
 				tmp_var_sumary = ele_name.text
 				tmp_hgvs_set = parse_var_summary(tmp_var_sumary)
 	if len(var_set.keys()) == 1:
 		allele_id = var_set.keys()[0]
 		if not re.search('NM\_',var_set[allele_id][0]) and re.search('NM\_',tmp_hgvs_set[0]): 
 			var_set[allele_id] = tmp_hgvs_set 					
#----condition parsing------------------------- 					
	trait_set = ref_assertion.find('TraitSet')
	for trait in trait_set.findall('Trait'):
		condition_id = int(trait.attrib['ID'])
		for name in trait.findall('Name'):
			type = name.find('ElementValue').attrib['Type']
			condition = name.findtext('ElementValue').replace('"','')
			if not type == 'Preferred':
				continue
			condition_set.append((condition_id,condition))
 	return var_id,var_type,var_name,var_set,condition_set
 			
def parse_clinical_assertion(clinvar_assertion):		
	submit_acc = clinvar_assertion.find('ClinVarAccession').attrib['Acc'][0:]
	submitter = clinvar_assertion.find('ClinVarSubmissionID').attrib['submitter'][0:]
	submit_lastdate = clinvar_assertion.find('ClinVarAccession').attrib['DateUpdated'][0:]
	record_status = clinvar_assertion.findtext('RecordStatus')[0:]
	
	try:
		significance = clinvar_assertion.find('ClinicalSignificance').findtext('Description')[0:]
	except:
		significance = 'Not provided'	
	inheritance_mode = 'Not Provided'	
	attributeset = clinvar_assertion.find('AttributeSet')	
	try:
		attributeset[0]
		for attribute in attributeset.findall('Attribute'):
			if attribute.attrib['Type'] == 'ModeOfInheritance':
				inheritance_mode = attribute.text
	except:
		inheritance_mode = 'Not Provided'
	review_status = clinvar_assertion.find('ClinicalSignificance').findtext('ReviewStatus')[0:]
	sample_origin = clinvar_assertion.find('ObservedIn').find('Sample').findtext('Origin').lower().capitalize()
	collection_method = clinvar_assertion.find('ObservedIn').find('Method').findtext('MethodType')[0:]	
	return submit_acc,record_status,significance,submitter,submit_lastdate,sample_origin,review_status,collection_method,inheritance_mode

def vcf_parser(input_vcf):

	allowed_chroms = []
	for i in range(1,23):
		allowed_chroms.append('chr%s'%i)
	allowed_chroms += ['chrX','chrY','chrM']
	required_header_list = ['CHROM','POS','REF','ALT']
	desc_header = dict()
	dh_idx = 1
	raw_vcf = dict()
	vcf_idx = 1
	variant = dict()
	v_idx = 1
	header = 'null'
	header_idx = dict()
	if '.gz' in input_vcf:
		finp = gzip.open(input_vcf)
	else:
		finp = open(input_vcf)	
	for line in finp:
		if re.search('^##',line):
			if re.search('^##INFO',line.strip()):
				desc_header[dh_idx] = line.strip()
				dh_idx += 1
			continue
		if re.search('^#',line):
			header = line.strip()
			header_token = line.strip()[1:].split('\t')
			for req_header in required_header_list:
				try:
					header_idx[req_header]=header_token.index(req_header)
				except:
					errmsg = '\n#ERROR: required vcf header was not found! : %s\n'%req_header
					sys.stderr.write(errmsg+'\n');sys.exit(2)
			if 'INFO' in header_token:
				header_idx['INFO']=header_token.index('INFO')
			else:
				header_idx['INFO']=len(header_token)
				header = '%s\tINFO'%header					
			continue
		if line.strip() == '': # skipping empty lines
			continue	
		token = line.strip().split('\t')
		chrom,pos,ref,alt_tmp = token[header_idx['CHROM']],token[header_idx['POS']],token[header_idx['REF']],token[header_idx['ALT']]
						
		if not re.search('chr',chrom):
			chrom = 'chr'+chrom
		if not chrom in allowed_chroms:
			continue
		alt_idx = header_idx['ALT']
		allele_index = 0
		for alt in alt_tmp.split(','):
			allele_index += 1
			try:
				info = token[header_idx['INFO']]
				split_line_token = token[:alt_idx]+[alt]+token[alt_idx+1:]				
			except:
				info = ''
				split_line_token = token[:alt_idx]+[alt]+token[alt_idx+1:]+['']
			vcf_line = '\t'.join(split_line_token)
			raw_vcf[vcf_idx] = vcf_line
			variant[vcf_idx] = (chrom,int(pos),ref,alt,allele_index)
			vcf_idx += 1
	finp.close()
	if header == 'null':
		errmsg = '\n#ERROR: vcf header was not defined!\n'
		sys.stderr.write(errmsg+'\n');sys.exit(2)
		
	return desc_header,header,header_idx,raw_vcf,variant		

def seq_reverse(seq):
	tmp = []
	for nt in seq:
		tmp.append(nt)
	return ''.join(reversed(tmp))


def summary_var_info(var_info,chrom):
	pos_list = sorted(var_info.keys())
	start = pos_list[0]
	end = pos_list[-1]
	change = var_info[start]
	var_summary = '%s[%s:%s-%s]'%(change,chrom,start,end)
	return var_summary
	
		
class IndelPositioning():
	
	def __init__(self,reference_file):		
		self.prev_chrom = 'null'
		self.ref_fasta = pyfaidx.Fasta(reference_file) 
		
	def __del__(self):
		sys.stderr.write('Done.\n')
				
	def define_variant(self,variant_tmp):
		
		variant = variant_tmp[:4]
		self.chrom,self.pos,self.ref,self.alt = variant
		if self.chrom == 'chrMT':
			self.chrom = 'chrM'
		self.variant = '\t'.join(['%s']*len(variant))%variant
		self.refend = self.pos+len(self.ref)-1
		if not self.chrom == self.prev_chrom:
			
			if not self.prev_chrom == 'null':
				sys.stderr.write('Done.\n')
			sys.stderr.write('Processing in %s ... '%self.chrom)
			self.chrom_seq = str(self.ref_fasta[self.chrom]).upper()
			self.prev_chrom = self.chrom
		
		
	def positioning_call(self):
		
		if (self.ref == self.alt) or (self.alt == '.'):
			var_summary = '%s>%s[%s:%s-%s]'%(self.ref,self.alt,self.chrom,self.pos,self.pos)
			return var_summary
		
		ref_allele = self.chrom_seq[self.pos-1:self.refend]
		
		if not ref_allele == self.ref:
			errmsg = '#ERROR: vcf file and referenece sequence are not matched! (%s:%s-%s)'%(self.chrom,self.pos,self.alt)
			sys.stderr.write(errmsg+'\n')
			var_summary = 'ref not matched'
			return var_summary
		
		
		self.var_trimming()		
		self.define_vartype()						
		
		if self.var_type == 'INS':
			
			var_info = self.search_insertion()
			try:			
				var_summary = summary_var_info(var_info,self.chrom)
			except:
				return 'null'
		elif self.var_type == 'DEL':
			var_info = self.search_deletion()
			try:
				var_summary = summary_var_info(var_info,self.chrom)
			except:
				return 'null'
		else:
			var_info = 'SNP'
			var_summary = '%s>%s[%s:%s-%s]'%(self.ref,self.alt,self.chrom,self.pos,self.pos)
						
		return var_summary
				
	def var_trimming(self): # left alignment
#---reverse trimming-------		
		tmp_rev_ref = seq_reverse(self.ref)
		tmp_rev_alt = seq_reverse(self.alt)
		max_rev_complen = min([len(self.ref),len(self.alt)])
		for i in range(max_rev_complen-1):
			ref_rev_nt = seq_reverse(self.ref)[i]
			alt_rev_nt = seq_reverse(self.alt)[i]
			if ref_rev_nt == alt_rev_nt:
				tmp_rev_ref = tmp_rev_ref[1:]
				tmp_rev_alt = tmp_rev_alt[1:]
			else:
				break
		tmp_ref = seq_reverse(tmp_rev_ref)
		tmp_alt = seq_reverse(tmp_rev_alt)
		
#---forward trimming--------
		conv_ref = tmp_ref
		conv_alt = tmp_alt
		max_forw_complen = min([len(tmp_ref),len(tmp_alt)])
		for j in range(max_forw_complen-1):
			ref_forw_nt = tmp_ref[j]
			alt_forw_nt = tmp_alt[j]
			if ref_forw_nt == alt_forw_nt:
				conv_ref = conv_ref[1:]
				conv_alt = conv_alt[1:]
				self.pos += 1
			else:
				break
		self.ref = conv_ref
		self.alt = conv_alt

	def define_vartype(self):
		
		if self.ref[0] == self.alt[0]: ## indel << after trimming
			if len(self.ref) > len(self.alt):
				self.var_type = 'DEL'
			else:
				self.var_type = 'INS'
		else:			
			if len(self.ref)*len(self.alt) == 1:
				self.var_type = 'SNV'
			else:
				self.var_type = 'COMPLEX'

	def search_insertion(self): 
		normalized_pos = 4000
		shift_pos_raw = self.pos+1
		insert_seq = self.alt[1:]
		full_strand = self.chrom_seq
		
# seq trimming and repositioning ( +- 10000 bp >> 20001 bp ) 

		if self.pos > normalized_pos:
			shift_pos = normalized_pos+1 ## normalized
		else:
			shift_pos = shift_pos_raw
		left_pos = shift_pos
		right_pos = shift_pos
		cor_pos = shift_pos_raw-shift_pos
		raw_strand = full_strand[(shift_pos_raw-(shift_pos-1)-1):(shift_pos_raw+normalized_pos)]
		
		control_strand = raw_strand[:left_pos-1]+insert_seq+raw_strand[left_pos-1:]

		insert_info = dict()
		insert_info[shift_pos+cor_pos] = '+%s'%insert_seq
#left search-------------

								
		while left_pos > 0:
			left_pos -= 1
			if raw_strand[left_pos-2] == 'N': # python index - 1
				break
			test_l_strand = raw_strand[:left_pos-1]+'X'*len(insert_seq)+raw_strand[left_pos-1:]
			blank_index = test_l_strand.index('X')
			comp_l_idx = blank_index-1
			comp_r_idx = blank_index+len(insert_seq)
			if (control_strand[comp_l_idx] == test_l_strand[comp_l_idx]) and (control_strand[comp_r_idx] == test_l_strand[comp_r_idx]): 
				insert_frag = control_strand[comp_l_idx+1:comp_r_idx]
				test_l_strand = raw_strand[:left_pos-1]+insert_frag+raw_strand[left_pos-1:]
				insert_info[left_pos+cor_pos] = '+%s'%insert_frag
				control_strand = test_l_strand
			else:
				break
		control_strand = raw_strand[:right_pos-1]+insert_seq+raw_strand[right_pos-1:]
		
#right search --------------
		
		while right_pos < len(raw_strand):
			right_pos += 1
			try:
				raw_strand[right_pos+len(insert_seq)-2]
			except:
				return 'null'
#				print raw_strand
#				print self.pos
#				print right_pos+len(insert_seq)-2
#				sys.exit()
			if raw_strand[right_pos+len(insert_seq)-2] == 'N':
				break
			test_r_strand = raw_strand[:right_pos-1]+'X'*len(insert_seq)+raw_strand[right_pos-1:]
			blank_index = test_r_strand.index('X')
			comp_l_idx = blank_index-1
			comp_r_idx = blank_index+len(insert_seq)

			if (control_strand[comp_l_idx] == test_r_strand[comp_l_idx]) and (control_strand[comp_r_idx] == test_r_strand[comp_r_idx]):				
				insert_frag = control_strand[comp_l_idx+1:comp_r_idx]
				test_r_strand = raw_strand[:right_pos-1]+insert_frag+raw_strand[right_pos-1:]
				insert_info[right_pos+cor_pos] = '+%s'%insert_frag
				control_strand = test_r_strand
			else:
				break
			
		return insert_info
							
	def search_deletion(self):
		shift_pos = self.pos+1
		delete_seq = self.ref[1:]
		raw_strand = self.chrom_seq
		left_pos = shift_pos
		right_pos = shift_pos
		delete_info = dict()
		delete_info[shift_pos] = '-%s'%delete_seq
#left search------------- test : left , compare : previous variant						
		while left_pos > 0:
			left_pos -= 1
			if raw_strand[left_pos-1] == 'N': # python index - 1
				break
			test_frag = raw_strand[left_pos-1:left_pos+len(delete_seq)][0]
			compare_frag = raw_strand[left_pos-1:left_pos+len(delete_seq)][-1]
			delete_frag = raw_strand[left_pos-1:left_pos+len(delete_seq)-1]					
			if test_frag == compare_frag:
				delete_info[left_pos] = '-%s'%delete_frag
			else:
				break
			
#right search -------------- test : previous variant compare : right >> stack position+1
		while right_pos < len(raw_strand):			
			if raw_strand[right_pos+len(delete_seq)-2] == 'N':
				break
			test_frag = raw_strand[right_pos-1:right_pos+len(delete_seq)][-1]
			compare_frag = raw_strand[right_pos-1:right_pos+len(delete_seq)][0]
			delete_frag = raw_strand[right_pos:right_pos+len(delete_seq)]
			if test_frag == compare_frag:
				delete_info[right_pos+1] = '-%s'%delete_frag
			else:
				break
			right_pos += 1
			
		return delete_info
											
class XmlParser():
	
	def __init__(self):
		self.variant_to_allele = dict()
		self.allele_to_variant = dict()
		self.condition_table = dict()
		self.variation_info = dict()
		self.allele_info = dict()
						
	def collect_clinvar_set(self,clinvarset):		
		def summary_var_set(var_id,var_set):
			allele_ids = []
			hgvs_forms = []
			for allele_id in sorted(var_set.keys()):
				allele_ids.append(str(allele_id))
				nm_id,hgvs_c,hgvs_p = var_set[allele_id]
				try:
					gene_name = transcript_list[nm_id]
				except:
					gene_name = '-'				
				hgvs_form = '%s|%s|%s|%s'%(nm_id,gene_name,hgvs_c,hgvs_p)
				hgvs_forms.append(hgvs_form)				
			return ','.join(allele_ids),','.join(hgvs_forms)
				
		ref_assertion = clinvarset.find('ReferenceClinVarAssertion')
		clinical_assertion_set = clinvarset.findall('ClinVarAssertion')			
		var_id,var_type,var_name,var_set,condition_set = parse_var_set(ref_assertion)		
		try:
			nm_id = re.findall('NM\_[0-9]+',var_name)[0]
			gene_name = transcript_list[nm_id]
		except:
			nm_id,gene_name = '-','-'
		allele_ids,hgvs_forms = summary_var_set(var_id,var_set)		
		hgvs_forms_set = hgvs_forms.split(',')
		allele_id_set = allele_ids.split(',')
		for i in range(len(allele_id_set)):
			allele_id = int(allele_id_set[i])
			hgvs_form = hgvs_forms_set[i]
			if not allele_id in self.allele_info:
				self.allele_info[allele_id] = hgvs_form
		allele_num = len(var_set.keys())		
		for allele_id in var_set.keys():
			if not allele_id in self.allele_to_variant:
				self.allele_to_variant[allele_id] = dict()
			if not var_id in self.allele_to_variant[allele_id]:
				self.allele_to_variant[allele_id][var_id] = var_type	
		if not var_id in self.variant_to_allele:
			self.variant_to_allele[var_id] = np.array(sorted(var_set.keys()))			
		for clinvar_assertion in clinical_assertion_set:
			submit_acc,record_status,significance,submitter,submit_lastdate,sample_origin,review_status,collection_method,inheritance_mode = parse_clinical_assertion(clinvar_assertion)
			
## Manual patch -----------------------------------------			

# significance : Pathogenic / Likely pathogenic / Uncertain significance / Likely benign / Benign
			
			if submit_acc == 'SCV000264384':
				significance = 'Benign'
			elif submit_acc == 'SCV000024411':
				significance = 'Benign'
			elif submit_acc == 'SCV000268745':
				significance = 'Likely benign'
			elif submit_acc == 'SCV000800200':
				significance = 'Uncertain significance'
			elif submit_acc == 'SCV000915634':
				significance = 'Uncertain significance'
## -------------------------------------------------------			
			submitter = submitter.strip()
			for condition_id,condition in condition_set:
				uniq_key = '%s:%s'%(submit_acc,condition_id)
# clinvar info stack ----------------				
				info = (var_type,var_name,allele_num,significance,inheritance_mode,condition_id,sample_origin,submitter,review_status,collection_method,submit_lastdate)
#------------------------------------												
				if not var_id in self.variation_info:
					self.variation_info[var_id] = dict()
				self.variation_info[var_id][uniq_key] = info
				if not condition_id in self.condition_table:
					self.condition_table[condition_id] = condition
	def save(self,outfile_name):
		np.savez_compressed(outfile_name,
			variant_to_allele = self.variant_to_allele,
			allele_to_variant = self.allele_to_variant,
			condition_table = self.condition_table,
			variation_info = self.variation_info,
			allele_info = self.allele_info)

if __name__=="__main__":

	if len(sys.argv)==1:
		usage('');sys.exit(2)

	input_type_names	=	['xml','vcf']
	try:
		input_type_opts,input_type_args	=	getopt.getopt(sys.argv[1:3],'T:')
	except getopt.GetoptError, error:
		errmsg	=	'#ERROR: Input type was not defined'
		usage(errmsg);sys.exit(2)
	try:
		input_type	=	input_type_opts[0][1]
		if not input_type in input_type_names:
			errmsg	=	'#ERROR: Input type "%s" was not recognized'%input_type
			usage(errmsg);sys.exit(2)
	except IndexError:
		errmsg	=	'#ERROR: Input type was not defined'
		usage(errmsg);sys.exit(2)
	rest_argvs		=	sys.argv[3:]


	start_page = '''		
#####################################################
## Clinvar DB Parser v2.0 (by youngjun 2018.12.18) ##
#####################################################
'''
	
#-------Xml parse-----------	
	if input_type == 'xml':
		s_time = time.time()
		try:
			opts,args	=	getopt.getopt(rest_argvs,'x:o:') # -x Clinvar Xml full release -o output name
		except getopt.GetoptError,err:
			errmsg	=	'#ERROR : %s'%str(err)
			usage(errmsg);sys.exit(2)
		xml_file,output_name = ['Null']*2
		parsed_options		=	[]
		parsed_args		=	[]
		required_list		=	['-x','-o']
		for opt,arg in opts:
			parsed_options.append(opt)
			parsed_args.append(arg)
			if opt	== '-x':
				xml_file = arg
			elif opt ==	'-o':
				output_name = arg
		required_option_check(required_list,parsed_options)
		input_file_list			=	[xml_file]
		option_file_check(parsed_options,parsed_args,input_file_list)
		
		sys.stderr.write(start_page)
		sys.stderr.write('''
##               Xml Parsing Mode                  ##
#####################################################
'''.strip()+'\n\n')

		sys.stderr.write('Xml file loading ...')	 		
		if re.search('\.gz',xml_file):
			fxml = gzip.open(xml_file)
		else:
			fxml = open(xml_file)
		sys.stderr.write(' Done.\n')
		sys.stderr.write('Transcript info loading ...')
		transcript_list = load_transcript_list()
		sys.stderr.write(' Done.\n')
		sys.stderr.write('Xml tree parsing ...')
		tree	=	ET.parse(fxml)
		root	=	tree.getroot()
		fxml.close()
		sys.stderr.write(' Done.\n')
		gene_compare = dict()
		data = dict()
		XP = XmlParser()
		sys.stderr.write('Data collecting ...')
		for clinvarset in root.findall('ClinVarSet'):  ## clinvar set -> variant for 1 condition including submittions
			XP.collect_clinvar_set(clinvarset)
		sys.stderr.write(' Done.\n')
		sys.stderr.write('Data saving ...')
		XP.save(output_name)
		del XP
		sys.stderr.write(' Done.\n')
		end_text = '''
Elapsed time : %s mins                     
#####################################################
'''%(format(((time.time()-s_time)/60),'0.2f'))		
		sys.stderr.write('\n'+end_text.strip()+'\n')		
	
	elif input_type == 'vcf':
		try:
			opts,args	=	getopt.getopt(rest_argvs,'R:v:o:') # -x Clinvar Xml full release -o output name
		except getopt.GetoptError,err:
			errmsg	=	'#ERROR : %s'%str(err)
			usage(errmsg);sys.exit(2)
		ref_fasta,clinvar_vcf,output_name = ['Null']*3
		parsed_options		=	[]
		parsed_args		=	[]
		required_list		=	['-R','-v','-o']
		for opt,arg in opts:
			parsed_options.append(opt)
			parsed_args.append(arg)
			if opt	== '-R':
				ref_fasta = arg
			elif opt ==	'-v':
				clinvar_vcf = arg
			elif opt == '-o':
				output_name = arg
		required_option_check(required_list,parsed_options)
		input_file_list			=	[ref_fasta,clinvar_vcf]
		option_file_check(parsed_options,parsed_args,input_file_list)
		
		sys.stderr.write(start_page)
		sys.stderr.write('''
##               VCF Parsing Mode                  ##
#####################################################
'''.strip()+'\n\n')
		
		total_process_time = time.time()
		sys.stderr.write('clinvar vcf file loading ...')
		desc_header,header,header_idx,raw_vcf,variant = vcf_parser(input_vcf=clinvar_vcf)
		desc_header_idx_max = max(desc_header.keys())

		extra_desc_header_idx = desc_header_idx_max+1
		desc_header[extra_desc_header_idx] = '##INFO=<ID=VARID,Number=1,Type=Integer,Description="ClinVar Variation ID">'
		extra_desc_header_idx += 1
		desc_header[extra_desc_header_idx] = '##INFO=<ID=CLNSIGREFVAR,Number=1,Type=String,Description="\
Clinical significance for a reference allele type at this variant position. Reported as pairs of VariationID:clinical significance.\
">'
		
		sys.stderr.write(' Done.\n')
		info_idx = header_idx['INFO']
		varid_idx = 2
		upv_to_info = dict()
		upv_to_alleleid = dict()
		ref_allele_info = dict()
		pre_chrom = 'null'
		sys.stderr.write('UPV calculation started.\n\n')
		IP = IndelPositioning(reference_file=ref_fasta)

		for vcf_idx in sorted(variant.keys()):
			chrom,pos,ref,alt,allele_index = variant[vcf_idx] 
			raw_token = raw_vcf[vcf_idx].strip().split('\t')
			var_id = raw_token[varid_idx]			
			info = '%s;VARID=%s'%(raw_token[info_idx],var_id)
			IP.define_variant(variant[vcf_idx])			
			upv = IP.positioning_call()
						
			if upv == 'null':
				continue
			if upv == 'RefNotMatched':
				continue			
			allele_id = int(re.findall('ALLELEID=[0-9]+',info)[0].split('=')[1])
			try:
				clnsig = re.findall('CLNSIG=[A-z,a-z,/,_]+',info)[0].replace(',_',',').split('=')[1]
			except:
				clnsig = 'not_provided'			
			
			upv_to_info[upv] = info
			upv_to_alleleid[upv] = allele_id
			
			if (alt == ref) or (alt == '.'):
				ref_allele_info[(chrom,pos,ref)] = 'CLNSIGREFVAR=%s:%s'%(var_id,clnsig)
		del IP
#------ data save---------		

		sys.stderr.write('\nSaving data ... ')
		np.savez_compressed(output_name,
			desc_header = desc_header,
			upv_to_info = upv_to_info,
			upv_to_alleleid = upv_to_alleleid,
			ref_allele_info = ref_allele_info
		)
		
		sys.stderr.write('Done.\n')
		total_elapsed_hours = format((time.time()-total_process_time),'0.2f')					
		end_text = '''
Total elapsed time : %s secs                    
#####################################################
'''%total_elapsed_hours		
		sys.stderr.write('\n'+end_text.strip()+'\n')		
		
