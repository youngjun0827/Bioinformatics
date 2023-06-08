#!/data/youngjun/local/bin/python

import sys
import os
import time
import getopt
import re
import numpy as np
import pysam
import pickle
import subprocess
import numpy

sliding_size = 10000
window_size = 3000000
sliding_size_name = '%skb'%int(sliding_size/1000)
window_size_name = '%sMb'%int(window_size/1000000)


def convert_numlist(input_list):
	output_list = list()
	for comp in input_list:
		try:
			converted_comp = eval(comp)
		except NameError:
			converted_comp = comp
		output_list.append(converted_comp)
	return output_list

def gc_convert(raw_gc):			
	tmp_gc = float(format(raw_gc,'0.2f'))
	if raw_gc > tmp_gc:
		conv_gc = round(tmp_gc+0.005,4)
	else:
		conv_gc = round(tmp_gc,4)
	return conv_gc

def MAD(values):
	median = np.median(values)
	ads = list()
	for value in values:
		ads.append(abs(value-median))
	return np.median(ads)
	
def convert_bam(bam_file_name, min_shift=4, threshold=4, mapq=20, isize_min_filt=30, isize_max_filt=600):
	
	total_bin_length = 0

	# Flush the current stack of reads and return filtered count
	def flush(read_buff,counts):
		reads_passed = 0
		reads_filtout = 0
		stair_size = len(read_buff)
		if stair_size <= threshold or threshold < 0:
			reads_passed += stair_size
			for read in read_buff:
				location = read.pos / sliding_size
				counts[int(location)] += 1
		else:
			reads_filtout += len(read_buff)
		return reads_filtout,reads_passed

	# Calculation of fragment shifting position
	def cal_shifted_pos(read,prev_read):
		if read.isize>0:
			read_start_pos = read.pos
			read_end_pos = read.pos+read.qlen - 1
		else:
			read_start_pos = read.pos-read.qlen + 1
			read_end_pos = read.pos				
		if prev_read.isize > 0:
			prev_read_start_pos = prev_read.pos
			prev_read_end_pos= prev_read.pos+prev_read.qlen - 1
		else:
			prev_read_start_pos = prev_read.pos-prev_read.qlen + 1
			prev_read_end_pos = prev_read.pos										
		start_shift = read_start_pos - prev_read_start_pos
		end_shift = read_end_pos - prev_read_end_pos
		if start_shift*end_shift < 0: # different direction
			shifted_pos = abs(start_shift)+abs(end_shift)
		else:
			shifted_pos = max([abs(start_shift),abs(end_shift)])
		return shifted_pos					
		
	sam_file = pysam.AlignmentFile(bam_file_name, "rb")

	reads_total = 0
	reads_dup = 0
	reads_mapqf = 0
	reads_proper_pairf = 0
	reads_unmapped = 0
	reads_qcpassed = 0
	reads_samecluster = 0
	reads_analyzed = 0
	
	last_readisize = -1
	last_readpos = -1

	prev_shift = 0
	new_shift = 0
	
	chrom_window_bin_counts = dict()
	chrom_bin_maxindex = dict()
	chrom_length = dict()
	for i in range(1,23):
		chrom_window_bin_counts['chr%s'%i] = None
	for sex_chrom in ['chrX','chrY']:
		chrom_window_bin_counts[sex_chrom] = None
	for index,chrom in enumerate(sam_file.references):

		chrom_name = chrom
		if not re.search('chr',chrom_name):
			chrom_name = 'chr%s'%chrom_name
		if chrom_name not in chrom_window_bin_counts:
			continue


		max_index = (sam_file.lengths[index]/sliding_size)+1 # 1-based
		
		chrom_length[chrom_name] = sam_file.lengths[index]
		chrom_bin_maxindex[chrom_name] = max_index
		
		counts = np.zeros(int(sam_file.lengths[index] / float(sliding_size) + 1))  # count list with zero value , 32bit->64fit (faster)
				
		read_buff = []
		sam_iter = sam_file.fetch(chrom)
		try:
			prev_read = sam_iter.next()
		except StopIteration:
			continue

		for read in sam_iter:
			fliter_flag = 0			
			reads_total += 1
			#>>> process in common
			if (last_readpos == read.pos) and (last_readisize == read.isize): # dup failed
				reads_dup += 1
				fliter_flag	+=	1
			else: # dup-passed
				if read.mapping_quality < mapq: # mapq failed
					reads_mapqf	+=	1
					fliter_flag	+=	1
				if read.is_paired: # paired-end
					if (not read.reference_name == read.next_reference_name) or (abs(read.isize)<isize_min_filt) or (abs(read.isize)>isize_max_filt):
						reads_proper_pairf	+=	1
						fliter_flag	+=	1
				if read.is_unmapped: # unmapped
					reads_unmapped	+=	1
					fliter_flag	+=	1

			last_readpos = read.pos
			last_readisize = read.isize
			#<<<
									
			if fliter_flag > 0:				
				continue			
			# QC passed read only
			reads_qcpassed	+=	1
						
											
			shifted_pos = cal_shifted_pos(read,prev_read)										
			if shifted_pos > min_shift:
				reads_filtout,reads_passed = flush(read_buff, counts)			
				reads_samecluster += reads_filtout
				reads_analyzed += reads_passed				
				read_buff = []

			read_buff.append(read)
			prev_read = read

		reads_filtout,reads_passed = flush(read_buff, counts)		
		reads_samecluster += reads_filtout
		reads_analyzed += reads_passed
		
		chrom_window_bin_counts[chrom_name] = counts
		
	qual_info = {'total_reads':reads_total,
				 'mapped':sam_file.mapped,
				 'unmapped':reads_unmapped,
				 'filter_mapq':reads_mapqf,
				 'filter_rmdup':reads_dup,
				 'pair_fail':reads_proper_pairf,
				 'qc_passed':reads_qcpassed,
				 'cluster_filtered':reads_samecluster,
				 'analyzed_read':reads_analyzed,
				 }

	return chrom_window_bin_counts,chrom_bin_maxindex,chrom_length,qual_info	
	
class BamProcessing():
	
	def __init__(self,bamfilename,isize_min_filt=30,isize_max_filt=600):
		self.chrom_bin_raw_counts,self.chrom_bin_maxindex,self.chrom_length,self.bam_qual_info = convert_bam(bamfilename,isize_min_filt=isize_min_filt,isize_max_filt=isize_max_filt)
		predata = np.load('predata.data')
		self.filter_region = list()
		for filt_region in predata['filter_region']:
			chrom_name,bin_start,bin_end = convert_numlist(filt_region)
			location = (chrom_name,bin_start,bin_end)
			self.filter_region.append(location)

		self.bin_feature_matrix = predata['bin_feature_matrix']
		
 				
		self.bin_location_index = dict()
		self.bin_location = dict()

		self.chrom_names = []
		for i in range(1,23):
			self.chrom_names.append('chr%s'%i)
		self.chrom_names += ['chrX','chrY']
		
							
	def raw_bin_counting(self):
		self.full_bin_raw_counts = np.zeros(sum(self.chrom_bin_maxindex.values()))
		self.full_50kb_bin_raw_counts = dict()
		bin_index = 0
		for chrom_name in self.chrom_names:
			self.full_50kb_bin_raw_counts[chrom_name] = dict()
			bin_50kb_index = -1
			chrom_len = self.chrom_length[chrom_name]
			for i in range(self.chrom_bin_maxindex[chrom_name]):
				bin_start = i*sliding_size
				bin_end = bin_start+sliding_size
				if bin_end > chrom_len:
					bin_end = chrom_len
				location = (chrom_name,bin_start,bin_end) 	
				bin_50kb_start = (bin_start/50000)*50000
				if bin_start == bin_50kb_start:
					bin_50kb_index += 1
				try:
					raw_count = self.chrom_bin_raw_counts[chrom_name][i]
				except:
					raw_count = 0.0
				if not bin_50kb_index in self.full_50kb_bin_raw_counts[chrom_name]:
					self.full_50kb_bin_raw_counts[chrom_name][bin_50kb_index] = 0.0
				self.full_50kb_bin_raw_counts[chrom_name][bin_50kb_index] += raw_count
					
				if location in self.filter_region:
					raw_count *= 0				
				self.bin_location_index[location] = bin_index
				self.bin_location[bin_index] = location
				self.full_bin_raw_counts[bin_index] += raw_count
				bin_index += 1
								
		self.bin_max_index = bin_index
				
	def count_outlier_filtering(self):		
		self.filt_bin_counts = np.zeros(self.bin_max_index)		
		non_zero_counts = sorted(self.full_bin_raw_counts[np.nonzero(self.full_bin_raw_counts)])
		
		avr = np.average(non_zero_counts)
		std = np.std(non_zero_counts)
	
		
		count_high_index = int(len(non_zero_counts)*(1-(0.1/100)))
		if count_high_index >= len(non_zero_counts):
			count_high_index = len(non_zero_counts)-1
				
		high_count = non_zero_counts[count_high_index]
		
				
		for bin_index in range(self.bin_max_index):
			location = self.bin_location[bin_index]
			raw_count = self.full_bin_raw_counts[bin_index]
			if (raw_count > high_count):
				continue
			self.filt_bin_counts[bin_index] += raw_count
			

	def construct_window_bin(self):
		self.window_feature = dict()
		self.window_location = dict()
		self.full_chrom_window_index = dict()
					
		window_index = 0
		start = time.time()		
		for chrom_name in self.chrom_names:

			#---- 3mb window-----
			self.full_chrom_window_index[chrom_name] = dict()
			end_flag = 0
			chrom_length = self.chrom_length[chrom_name]
			max_index = self.chrom_bin_maxindex[chrom_name]
			chrom_window_start = -1
			chrom_window_end = chrom_length
			for i in range(max_index):
				if end_flag > 0:
					continue				
				window_start = i*sliding_size
				window_end = window_start+window_size
				
				if window_end > chrom_length:
					window_end = chrom_length
					end_flag += 1
					if window_size == sliding_size:
						continue
					
				window_location = (chrom_name,window_start,window_end)
								
				self.window_location[window_index] = window_location
				
				bin_start_location = (chrom_name,window_start,window_start+sliding_size)
				tmp_bin_start_pos = (window_end/sliding_size)*sliding_size
				if window_end == tmp_bin_start_pos:
					bin_end_location = (chrom_name,tmp_bin_start_pos-sliding_size,window_end)					
				else:
					bin_end_location = (chrom_name,tmp_bin_start_pos,window_end)

				bin_start_index = self.bin_location_index[bin_start_location]
				bin_end_index = self.bin_location_index[bin_end_location]			
				window_len = window_end-window_start
				window_feature = self.bin_feature_matrix[bin_start_index:bin_end_index+1].sum(axis=0)[0]
				window_num_gc,window_num_nt,window_num_n,mappable_cnt = window_feature
				window_mappability = float(mappable_cnt)/window_len
				if window_num_nt == 0:
					window_pct_gc = 0.0
				else:
					window_pct_gc = float(window_num_gc)/window_num_nt
				window_pct_n = float(window_num_n)/window_len
				self.window_feature[window_index] = (window_pct_gc, window_pct_n, window_mappability,(bin_start_index,bin_end_index+1))
				window_index += 1

		self.window_max_index = window_index

	def window_mappability_correction(self): # <0.5 filter & correction & calculation gc correction weight
		self.window_mappability_corrected_counts = np.zeros(self.window_max_index)
		window_counts_collect_per_gc = dict()
		autosomal_global = dict()
		for window_index in range(self.window_max_index):
			window_feature = self.window_feature[window_index]
			window_pct_gc, window_pct_n, window_mappability,window_bin_index_set = window_feature
			window_bin_f_index,window_bin_r_index = window_bin_index_set
			window_counts = sum(self.full_bin_raw_counts[window_bin_f_index:window_bin_r_index])
			chrom_name,window_start,window_end = self.window_location[window_index]

			if (window_pct_n > 0.1) or (window_mappability < 0.5):
				mappability_corrected_window_count = 0.0
			else:
				mappability_corrected_window_count = window_counts/window_mappability			
			self.window_mappability_corrected_counts[window_index] += mappability_corrected_window_count
			if (not re.search('chr[X,Y,M]',chrom_name)) and ( mappability_corrected_window_count > 0): # autosome & non-zero
				ref_gc = gc_convert(window_pct_gc)
				if not ref_gc in window_counts_collect_per_gc:
					window_counts_collect_per_gc[ref_gc] = dict()
				window_counts_collect_per_gc[ref_gc][window_index] = mappability_corrected_window_count
				autosomal_global[window_index] = mappability_corrected_window_count
				
		self.gc_weight = dict()
		global_avr = np.average(autosomal_global.values())
		for ref_gc in sorted(window_counts_collect_per_gc.keys()):
			weight = global_avr/np.average(window_counts_collect_per_gc[ref_gc].values())
			self.gc_weight[ref_gc] = weight

	def window_gc_correction(self):
		self.normalized_chrom_counts = dict()
		self.normalized_window_counts = np.zeros(self.window_max_index)
		for window_index in range(self.window_max_index):
			window_feature = self.window_feature[window_index]
			window_pct_gc, window_pct_n, window_mappability,window_bin_index_set = window_feature
			ref_gc = gc_convert(window_pct_gc)
			chrom_name,window_start,window_end = self.window_location[window_index]
			mb_corrected_window_counts = self.window_mappability_corrected_counts[window_index]
			try:
				full_normalized_window_counts = self.gc_weight[ref_gc]*mb_corrected_window_counts
			except KeyError:
				full_normalized_window_counts = 0.0
			self.normalized_window_counts[window_index] += full_normalized_window_counts
			if not chrom_name in self.normalized_chrom_counts:
				self.normalized_chrom_counts[chrom_name] = 0.0
			if window_start%window_size == 0:
				self.normalized_chrom_counts[chrom_name] += full_normalized_window_counts

					
	def save_count(self,outfile_name):
		np.savez_compressed(outfile_name,
							window_feature = self.window_feature,
							window_location = self.window_location,
							raw_bin_counts = self.chrom_bin_raw_counts,
							bin_50kb_counts = self.full_50kb_bin_raw_counts,
							normalized_window_counts = self.normalized_window_counts,
							normalized_chrom_counts = self.normalized_chrom_counts,
							qual_info = self.bam_qual_info							
							)
							
	
	def exe(self,outfile_name):
		self.raw_bin_counting()
		self.count_outlier_filtering()
		self.construct_window_bin()
		self.window_mappability_correction()
		self.window_gc_correction()
		self.save_count(outfile_name)

if __name__ == '__main__':
	bamfilename	=	sys.argv[1]
	bamloc_header = bamfilename.split('.')[0]
	outfile_name = '%s.readcount'%bamloc_header
	processed_bam =	BamProcessing(bamfilename)
	processed_bam.exe(outfile_name)
	
