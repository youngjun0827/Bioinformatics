#!/data/youngjun/local/bin/python
#-*-coding:utf-8-*-

# This code was written with reference to a research paper "UPS-indel: a Universal Positioning System for Indels." (https://www.nature.com/articles/s41598-017-14400-1)


from unidecode import unidecode
import sys,os,pickle,string,unicodedata,gzip,time,getopt,datetime,re,pymysql
import pyfaidx
import numpy as np
import xml.etree.cElementTree as ET
reload(sys)
sys.setdefaultencoding('utf-8')

def usage(errmsg):
	title		=	'''
#################################################
## UPV annotator v1.0 (by youngjun 2018.12.20) ##
#################################################
'''
	sys.stderr.write(title)	
	if not errmsg == '':
		sys.stderr.write('\n'+errmsg+'\n')
	errtext	=	'''
usage: %s -R [reference fasta ] -v [input vcf] -c [parsed clinvar vcf npz (optional)]
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

##------------------------------------

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
				
	def define_variant(self,variant_tmp):
		
		variant = variant_tmp[:4]
		self.chrom,self.pos,self.ref,self.alt = variant
		if self.chrom == 'chrMT':
			self.chrom = 'chrM'
		self.variant = '\t'.join(['%s']*len(variant))%variant
		self.refend = self.pos+len(self.ref)-1
		if not self.chrom == self.prev_chrom:
			self.chrom_seq = str(self.ref_fasta[self.chrom]).upper()
			self.prev_chrom = self.chrom
		
		
	def positioning_call(self):
		
		if (self.ref == self.alt) or (self.alt == '.'):
			var_summary = '%s>%s[%s:%s-%s]'%(self.ref,self.ref,self.chrom,self.pos,self.pos)
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
			
			var_summary = summary_var_info(var_info,self.chrom)
		elif self.var_type == 'DEL':
			var_info = self.search_deletion()
			var_summary = summary_var_info(var_info,self.chrom)
		else:
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
		shift_pos_raw = self.pos+1
		insert_seq = self.alt[1:]
		full_strand = self.chrom_seq
		
# seq trimming and repositioning ( +- 10000 bp >> 20001 bp ) 

		if self.pos > 20000:
			shift_pos = 20001 ## normalized
		else:
			shift_pos = shift_pos_raw
		left_pos = shift_pos
		right_pos = shift_pos
		cor_pos = shift_pos_raw-shift_pos
		raw_strand = full_strand[(shift_pos_raw-(shift_pos-1)-1):(shift_pos_raw+20000)]
		
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
					

def parse_info(raw_info):
	info_keys = []
	line_info = dict()
	for info in raw_info.split(';'):
		if re.search('=',info):
			info_name,info_text = info.split('=')
		else:
			info_name = info
			info_text = ''
		info_keys.append(info_name)
		line_info[info_name] = info_text
	return info_keys,line_info


if __name__=="__main__":
	if len(sys.argv) == 1:
		usage('');sys.exit(2)
	try:
		opts,args	=	getopt.getopt(sys.argv[1:],'R:v:c:') # -R : reference fasta -v input vcf -c upv-added clinvar vcf npz ( optional)
	except getopt.GetoptError,err:
		errmsg	=	'#ERROR : %s'%str(err)
		usage(errmsg);sys.exit(2)
	ref_fasta,input_vcf,clinvar_npz = ['Null']*3
	parsed_options		=	[]
	parsed_args		=	[]
	required_list		=	['-R','-v']
	for opt,arg in opts:
		parsed_options.append(opt)
		parsed_args.append(arg)
		if opt	== '-R':
			ref_fasta = arg
		elif opt ==	'-v':
			input_vcf = arg
		elif opt == '-c':
			clinvar_npz = arg
	required_option_check(required_list,parsed_options)
	input_file_list			=	[ref_fasta,input_vcf,clinvar_npz]
	option_file_check(parsed_options,parsed_args,input_file_list)
	desc_header,header,header_idx,raw_vcf,variant = vcf_parser(input_vcf)
	dh_max_idx = max(desc_header.keys())
	
	add_idx = dh_max_idx + 1
	desc_header[add_idx] = '##INFO=<ID=UPV,Number=.,Type=String,Description="Universal Position of Variant">'
	add_idx += 1
	desc_header[add_idx] = '##INFO=<ID=ALLELEINDEX,Number=1,Type=Integer,Description="\
Allele index from original variant calling. e.g. ref - 0, 1st alt - 1, 2nd alt - 2, etc.\
">' 		
	
	info_idx = header_idx['INFO']
	clinvar_flag = 0
	try:
		clinvar_flag = 1
		clinv_npz = np.load(clinvar_npz)		
		clinv_header = clinv_npz['desc_header'].item()
		upv_to_info = clinv_npz['upv_to_info'].item()
		upv_to_alleleid = clinv_npz['upv_to_alleleid'].item()
		ref_allele_info = clinv_npz['ref_allele_info'].item()
		add_idx += 1
		desc_header[add_idx] = '##INFO=<ID=VARID,Number=1,Type=Integer,Description="ClinVar Variation ID">'
		add_idx += 1
		desc_header[add_idx] = '##INFO=<ID=CLNSIGREFVAR,Number=1,Type=String,Description="\
Clinical significance for a reference allele type at this variant position. Reported as pairs of VariationID:clinical significance.\
">'
				
	except:
		clinv_npz = 'null' 
		clinv_header = 'null'
		upv_to_info = 'null'
		upv_to_alleleid = 'null'
		ref_allele_info = {}
		
	for header_i in sorted(desc_header.keys()):
		print desc_header[header_i]
	try:
		for clinv_header_i in sorted(clinv_header.keys()):
			header_text = clinv_header[clinv_header_i]
			if (re.search('##INFO',header_text)) and (not header_text in desc_header.values()):
				print header_text
		print header
	except:
		print header
	
	IP = IndelPositioning(reference_file=ref_fasta)
	for i in sorted(raw_vcf.keys()):
		raw_vcf_line = raw_vcf[i]
		token = raw_vcf_line.strip().split('\t')
		
		try:
			raw_info = token[info_idx]
		except:
			raw_info = ''

		info_keys,line_info = parse_info(raw_info)

		site_variant_set = variant[i]
		chrom,pos,ref,alt,allele_index = site_variant_set
		IP.define_variant(site_variant_set)
		upv = IP.positioning_call()
		added_info = '%s;UPV=%s;ALLELEINDEX=%s'%(raw_info,upv,allele_index)
		if (clinvar_flag == 1) and (upv in upv_to_info):
			added_info = '%s;%s'%(added_info,upv_to_info[upv])		
		ref_varset = (chrom,pos,ref)
		if ref_varset in ref_allele_info:
			added_info = '%s;%s'%(added_info,ref_allele_info[ref_varset])		
		added_info_keys,added_line_info = parse_info(added_info)
		for added_info_key in added_info_keys:
			if not added_info_key in info_keys:
				info_keys.append(added_info_key)
				line_info[added_info_key] = added_line_info[added_info_key]
		final_info_tmp = []
		for info_key in info_keys:
			info_value = line_info[info_key]
			if not info_value == '':
				final_info_tmp.append('%s=%s'%(info_key,info_value))
			else:
				final_info_tmp.append(info_key)
		final_info = ';'.join(final_info_tmp)
		token[info_idx] = final_info
		
		print '\t'.join(token)
			
	del IP	
	sys.exit()

	
