import sys,os,random,re


seqtk_path = '/data/Tools/seqtk/seqtk'
sickle_path = '/data/Tools/sickle-1.2/sickle'
aligner = '/data/Tools/align/bwa-0.6.2/bwa'
refSeq = '/data//Ref/bwa062/hg19.fasta'
samtools_path = '/data/Tools/samtools-1.2/samtools'
picard_path = '/data/Tools/picard-tools-1.102'
bedtools_path = '/data/Tools/bedtools-2.17.0/bin/bedtools'
	
platfm ='Illumina'
library = 'HG19'
qual_cut = 20
tagID = 'TEST'

mdset = ['1p36del','criduchat','digeorge','jacobsen','praderwilli_angelman']

test_number = 9
total_read_set = [10000000,15000000,20000000,25000000,30000000]
mother_samples = ['M1','M2']
fetus_samples = ['F2','F3']
ff_set = [4,6,8,10,12,14,16,18,20]

def calculate_read_size(ff,n_fetus,total_reads):
	if not (n_fetus == 1 or n_fetus == 2):
		print 'Wrong fetus number input: %s'%n_fetus;sys.exit()
	ff_rate = float(ff)/100
	maternal_reads = int(total_reads*(1-ff_rate)/2
	fetal_reads = int(total_reads*(ff_rate))/n_fetus/2
	return maternal_reads,fetal_reads

def run_seqtk(params):
	iter_num,ff,mother,fetus,total_reads = params
	test_tag ='rep%s'%iter_num
	run_id = '%s-%s-%sper-%s'%(mother,fetus,ff,test_tag)
	maternal_reads,fetal_reads	=	CALCULATE_READ_SIZE(ff,1,total_reads)
	fetal_alleleA_reads = fetal_reads/2
	fetal_alleleB_reads = fetal_reads/2
									
	mother_raw_fastq_read1	=	'RAW_FASTQ/%s_1.fastq.gz'%mother
	mother_raw_fastq_read2	=	'RAW_FASTQ/%s_2.fastq.gz'%mother
	fetus_raw_fastq_read1	=	'RAW_FASTQ/%s_1.fastq.gz'%fetus
	fetus_raw_fastq_read2	=	'RAW_FASTQ/%s_2.fastq.gz'%fetus

	maternal_seed = random.randint(10000000,99999999)
	fetal_alleleA_seed =	random.randint(10000000,99999999)
	fetal_alleleB_seed =	random.randint(10000000,99999999)
					
	mother_indiv_run_id = '%s-%s'%(run_id,mother)
	fetus_alleleA_indiv_run_id = '%s-%sA'%(run_id,fetus)
	fetus_alleleB_indiv_run_id = '%s-%sB'%(run_id,fetus)
					
	mother_outfastq_read1 = 'FASTQ/%s_1.fastq.gz'%mother_indiv_run_id
	mother_outfastq_read2 = 'FASTQ/%s_2.fastq.gz'%mother_indiv_run_id
	fetus_alleleA_outfastq_read1 = 'FASTQ/%s_1.fastq.gz'%fetus_alleleA_indiv_run_id
	fetus_alleleA_outfastq_read2 = 'FASTQ/%s_2.fastq.gz'%fetus_alleleA_indiv_run_id
	fetus_alleleB_outfastq_read1 = 'FASTQ/%s_1.fastq.gz'%fetus_alleleB_indiv_run_id
	fetus_alleleB_outfastq_read2 = 'FASTQ/%s_2.fastq.gz'%fetus_alleleB_indiv_run_id

					
	mother_fastq_read1_command = '%s sample -s %s %s %s | gzip > %s'%(seqtk_path,maternal_seed,mother_raw_fastq_read1,maternal_reads,mother_outfastq_read1) 
	mother_fastq_read2_command = '%s sample -s %s %s %s | gzip > %s'%(seqtk_path,maternal_seed,mother_raw_fastq_read2,maternal_reads,mother_outfastq_read2)
	fetus_alleleA_fastq_read1_command = '%s sample -s %s %s %s | gzip > %s'%(seqtk_path,fetal_alleleA_seed,fetus_raw_fastq_read1,fetal_alleleA_reads,fetus_alleleA_outfastq_read1)
	fetus_alleleA_fastq_read2_command = '%s sample -s %s %s %s | gzip > %s'%(seqtk_path,fetal_alleleA_seed,fetus_raw_fastq_read2,fetal_alleleA_reads,fetus_alleleA_outfastq_read2)
	fetus_alleleB_fastq_read1_command = '%s sample -s %s %s %s | gzip > %s'%(seqtk_path,fetal_alleleB_seed,fetus_raw_fastq_read1,fetal_alleleB_reads,fetus_alleleB_outfastq_read1)
	fetus_alleleB_fastq_read2_command = '%s sample -s %s %s %s | gzip > %s'%(seqtk_path,fetal_alleleB_seed,fetus_raw_fastq_read2,fetal_alleleB_reads,fetus_alleleB_outfastq_read2)
					
	print mother_fastq_read1_command
	print mother_fastq_read2_command
	print fetus_alleleA_fastq_read1_command
	print fetus_alleleA_fastq_read2_command
	print fetus_alleleB_fastq_read1_command
	print fetus_alleleB_fastq_read2_command
	
	out_set = [mother_indiv_run_id,fetus_alleleA_indiv_run_id,fetus_alleleB_indiv_run_id]
	return out_set

def fastq_trimming(out_set):
	for indiv_run_id in out_set:
		raw_fastq_read_1 = 'FASTQ/%s_1.fastq.gz'%indiv_run_id
		raw_fastq_read_2 = 'FASTQ/%s_2.fastq.gz'%indiv_run_id

		cut_fastq_read_1 = 'FASTQ/%s-cut_1.fastq.gz'%indiv_run_id
		cut_fastq_read_2 = 'FASTQ/%s-cut_2.fastq.gz'%indiv_run_id
		
		cut_command_read1 = 'zcat %s | cut -c -75 | gzip > %s'%(raw_fastq_read_1,cut_fastq_read_1)
		cut_command_read2 = 'zcat %s | cut -c -75 | gzip > %s'%(raw_fastq_read_2,cut_fastq_read_2)

		sickle_fastq_read1 = 'FASTQ/%s-cleaned_1.fastq'%indiv_run_id
		sickle_fastq_read2 = 'FASTQ/%s-cleaned_2.fastq'%indiv_run_id
		
		sickle_command = '%s pe -t sanger -n -q 20 -x -f %s -r %s -o %s -p %s -s FASTQ/%s_single.fastq'%(sickle_path,cut_fastq_read_1,cut_fastq_read_2,sickle_fastq_read1,sickle_fastq_read2,indiv_run_id)
		
		gzip_command_read1 = 'gzip %s'%sickle_fastq_read1
		gzip_command_read2 = 'gzip %s'%sickle_fastq_read2
		
		del_command = 'rm FASTQ/%s-cut_*.fastq.gz FASTQ/%s_single.fastq'%(indiv_run_id,indiv_run_id)
		
		print cut_command_read1
		print cut_command_read2
		print sickle_command
		print gzip_command_read1
		print gzip_command_read2
		print del_command
		
		
def make_bam(run_id):
	
	run_id_set = []
		
	align_read1_command = '%s aln -t 8 -n 2 -q %s %s FASTQ/%s-cleaned_1.fastq.gz > Indiv_bams/%s_1.sai'%(aligner,qual_cut,refSeq,run_id,run_id)
 	align_read2_command = '%s aln -t 8 -n 2 -q %s %s FASTQ/%s-cleaned_2.fastq.gz > Indiv_bams/%s_2.sai'%(aligner,qual_cut,refSeq,run_id,run_id) 	
	sampe_command = "%s sampe -r '@RG\tPL:%s\tID:%s\tSM:%s\tLB:%s' %s Indiv_bams/%s_1.sai Indiv_bams/%s_2.sai FASTQ/%s-cleaned_1.fastq.gz FASTQ/%s-cleaned_2.fastq.gz | %s  view -q %s -bS - | %s sort -m 16000000000 - Indiv_bams/%s.sorted.%s "%(aligner,platfm,tagID,run_id,library,refSeq,run_id,run_id,run_id,run_id,samtools_path,qual_cut,samtools_path,run_id,library)

	print align_read1_command
	print align_read2_command
 	print sampe_command

	
	if re.search('F[1-3]B',run_id):
		rmdup_command = 'java -jar %s/MarkDuplicates.jar  I=Indiv_bams/%s.sorted.%s.bam  O=Indiv_bams/%s.rmdup.%s.tmp.bam  M=Indiv_bams/%s.picard_metrics.txt  AS=true  REMOVE_DUPLICATES=true  VALIDATION_STRINGENCY=LENIENT TMP_DIR=./tmp_rmdup >  Indiv_bams/%s.rmdup.log '%(picard_path,run_id,library,run_id,library,run_id,run_id)
		bam_index_command = '%s index Indiv_bams/%s.rmdup.%s.tmp.bam'%(samtools_path,run_id,library)
		print rmdup_command
		print bam_index_command		
		for mdtarget in mdset:
			filt_bed = '%s.bed'%mdtarget
			deletion_command = '%s intersect -abam Indiv_bams/%s.rmdup.%s.tmp.bam -b %s -v > Indiv_bams/%s-%s.rmdup.%s.bam'%(bedtools_path,run_id,library,filt_bed,run_id,mdtarget,library)
			delbam_index_command = '%s index Indiv_bams/%s-%s.rmdup.%s.bam'%(samtools_path,run_id,mdtarget,library)
			print deletion_command  
			print delbam_index_command
			
			run_id_set.append('Indiv_bams/%s-%s.rmdup.HG19.bam'%(run_id,mdtarget))		
	else:
		rmdup_command = 'java -jar %s/MarkDuplicates.jar  I=Indiv_bams/%s.sorted.%s.bam  O=Indiv_bams/%s.rmdup.%s.bam  M=Indiv_bams/%s.picard_metrics.txt  AS=true  REMOVE_DUPLICATES=true  VALIDATION_STRINGENCY=LENIENT TMP_DIR=./tmp_rmdup >  Indiv_bams/%s.rmdup.log '%(picard_path,run_id,library,run_id,library,run_id,run_id)	
		bam_index_command = '%s index Indiv_bams/%s.rmdup.%s.bam'%(samtools_path,run_id,library)
		print rmdup_command
		print bam_index_command
 	
 		run_id_set.append('Indiv_bams/%s.rmdup.HG19.bam'%run_id)
 		
 	return run_id_set
 	

if __name__=="__main__":
	for iter_num in range(1,test_number+1):
		test_tag	=	'rep%s'%iter_num		
		for ff in ff_set:
			for mother in mother_samples:
				for fetus in fetus_samples:
					for total_reads in total_read_set:
						base_bam_set = []
						md_allele_set = []
						params = [iter_num,ff,mother,fetus,total_reads]
						out_set = run_seqtk(params)
						fastq_trimming(out_set)
						
						for indiv_run_id in out_set:
							out_bam_set = make_bam(indiv_run_id)
							if len(out_bam_set) == 1:
								base_bam_set.append(out_bam_set[0])
							else:
								md_allele_set = out_bam_set
						for md_del_bam in md_allele_set:
							md_name = md_del_bam.split('-')[5].split('.')[0]
							merge_bam_list = base_bam_set+[md_del_bam]
							merged_bam_name = 'Merged_bams/%s-%s-%s-%sper-%s.rmdup.HG19.bam'%(mother,fetus,md_name,ff,test_tag)
							merged_bam_head = merged_bam_name.split('.')[0]
							merge_command = '%s merge %s %s'%(samtools_path,merged_bam_name,' '.join(merge_bam_list))
							merge_index_command = '%s index %s'%(samtools_path,merged_bam_name)
	
							#---exe ---						
							print merge_command
							print merge_index_command
							
							#---------
						
						
											

						



					
