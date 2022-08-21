# import packages
import pysam
import logging
import argparse
import os
from weblogo import *

# parse argument
parser = argparse.ArgumentParser(description = "")
parser.add_argument('--input_bamfile', required = True, type = str, help = 'Input bam file.')
parser.add_argument('--index', type = str, required = False, help = 'Input index file of the bam file.')
parser.add_argument('--ref', required = True, type = str, help = 'Input the reference file.')
parser.add_argument('--output', required = False, help = 'Output fastq filename')
args = parser.parse_args()

# get path
outputPath = os.getcwd() + "/"

filename = args.input_bamfile.split(".")[:-1]
filename = "".join(filename)
filename = filename.split("/")[-1]
filename = outputPath + "/" + filename

# build error logger
errorLogger = logging.getLogger('errorLogger')
errorLogger.setLevel(logging.ERROR)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
errorLogger.addHandler(sh)

# build info logger to print out log information
infoLogger = logging.getLogger('infoLogger')
infoLogger.setLevel(logging.INFO)
fh = logging.FileHandler(f'{filename}.log')
fh.setFormatter(logging.Formatter('%(message)s'))
infoLogger.addHandler(fh)
infoLogger.info(f'The input Bam file is {args.input_bamfile}\nThe index file is {args.index}\nThe reference file is {args.ref}\nThe output File is {args.output}\n')

# sort bam file
sorted_bamfile = ''.join(args.input_bamfile.split('.')[:-1]) + ".sorted.bam"
os.system(f'samtools sort -n {args.input_bamfile} -o {sorted_bamfile}')

# find index file and reference file
if args.index is None:
	index_file = f"{args.input_bamfile}.bai"
else:
	index_file = args.index

try:
	open(index_file,'r')
except FileNotFoundError:
	errorLogger.error('The index file is not exist. Please double check.')
	raise SystemExit(1)

try:
	with open(args.ref, 'r') as reffile:
		rawfile = reffile.readlines()
	rawfile.pop(0)
	count = 0
	while count < len(rawfile):
		rawfile[count] = rawfile[count].strip()
		count += 1
	refseq = ''.join(rawfile)
	refseq.replace("\n","")
except FileNotFoundError:
	errorLogger.error('File {} cannot be opened for reading. Please double check.'.format(args.ref))
	raise SystemExit(1)

#open bam file
infoLogger.info(f'Opening bam file {sorted_bamfile} with the index {index_file} and the reference {args.ref}\n')
try:
	samfile = pysam.AlignmentFile(sorted_bamfile,"rb", reference_filename = args.ref, index_filename = index_file)
except FileNotFoundError:
	errorLogger.error('File {} cannot be opened for reading. Please double check.'.format(sorted_bamfile))
	raise SystemExit(1)

# start reading reads and find the mapped paired reads
mapped_count = 0
unmapped_count = 0
secondary_count = 0
supplementary_count = 0
paired_count = 0
primary_count = 0
mapped_paired_list = []
for r in samfile.fetch(until_eof = True):
	if r.is_secondary:
		secondary_count += 1
	elif r.is_supplementary:
		supplementary_count += 1
	elif r.is_unmapped:
		unmapped_count += 1
	else:
		primary_count += 1
		if r.is_paired:
                        paired_count += 1
                        if r.mate_is_unmapped is False:
                                mapped_count += 1
                                mapped_paired_list.append(r)
			

infoLogger.info(f'Count of primary alignment: {primary_count}\nCount of paired and mapped reads: {mapped_count}\nCount of secondary alignment: {secondary_count}\nCount of supplementary alignment: {supplementary_count}\nCount of unmapped reads: {unmapped_count}\n')

# function to identify if the paired reads are overlaped
def if_overlap(read1, read2):
	if read1.is_reverse and read2.is_reverse:
		return False
	elif not read1.is_reverse and read2.is_reverse:
		return False
	elif not read1.is_reverse and read2.is_reverse:
		if read2.reference_start <= (read1.reference_end-1) and read2.reference_end >= read1.reference_end:
			overlap_len = read1.reference_end - read2.reference_start
			return f"+{overlap_len}"
		else:
			return False
	elif not read2.is_reverse and read1.is_reverse:
		if read1.reference_start <= (read2.reference_end-1) and read1.reference_end >= read2.reference_end:
			overlap_len = read2.reference_end - read1.reference_start
			return f"-{overlap_len}"
		else:
			return False
	else:
		return False

# function to reverse the reverse reads
def resverse_seq(seq):
	seq_list = list(seq)
	seq_list = list(reversed(seq_list))
	re_seq = "".join(seq_list)
	new_seq = ""
	for nucleotide in re_seq:
		if nucleotide.upper() == "A":
			new_seq += "T"
		elif nucleotide.upper() == "T":
			new_seq += "A"
		elif nucleotide.upper() == "G":
			new_seq += "C"
		elif nucleotide.upper() == "C":
			new_seq += "G"
		elif nucleotide.upper() == "N":
			new_seq += "-"
		elif nucleotide == "-":
			new_seq += "-"
	return new_seq

def resverse(seq, qual):
	re_qual = list(reversed(qual))
	new_seq = resverse_seq(seq)
	return (new_seq, re_qual)


# function to align the paired reads with reference sequence
def read_aligner(read, ref_seq, strand):
	cigar = read.cigartuples
	ref_len = len(ref_seq)
	seq = read.get_forward_sequence()
	quality = read.get_forward_qualities()
	if strand == "-":
		seq = resverse(seq, quality)[0]
		quality = resverse(seq, quality)[1]
	pos = 0
	new_seq = ""
	new_quality = ""	
	for item in cigar:
		if item[0] == 4 or item[0] == 6:
			pos += item[1]
		elif item[0] == 0 or item[0] == 3 or item[0] == '=' or item[0] == 'X':
			new_seq += seq[pos:pos+item[1]]
			quality_list = quality[pos:pos+item[1]]
			for q in quality_list:
				qual = chr(int(q) + 33)
				new_quality += qual
			pos += item[1]
		elif item[0] == 1:
			pos += item[1]
			continue
		elif item[0] == 2:
			new_seq += "-" * item[1]
			new_quality += "*" * item[1]
	if strand == "+":
		add_len = ref_len - len(new_seq)
		new_seq = new_seq + " " * add_len
		new_quality = new_quality + " " * add_len
	elif strand == "-":
		add_len = ref_len - len(new_seq)
		new_seq =  " " * add_len + new_seq
		new_quality = " " * add_len + new_quality
	return (new_seq,new_quality)

# class to record each error
class error_record:
	def __init__(self, error_pos_read, error_upstream, error_type, ref_base, ref_pos):
		if len(error_upstream) < 5:
			error_upstream = "-" * (5-len(error_upstream)) + error_upstream
		self.error_pos = ref_pos
		self.error_upstream = error_upstream
		self.error_type = error_type
		self.ref_base = ref_base
		self.error_read_pos = error_pos_read

#function to record the insertion information
class ins_info:
	def __init__(self,read):
		ins = {}
		ins_pos = 0
		real_pos = 0
		self.read = read
		for i in read.cigartuples:
			if i[0] == 1:
				ins[ins_pos] = read.get_forward_sequence()[real_pos:(real_pos+i[1])]
				real_pos += i[1]
			elif i[0] == 4:
				real_pos += i[1]
			elif i[0] != 2:
				ins_pos += i[1]
				real_pos += i[1]
		self.ins = ins
	def len_before_position(self,position):
		length = 0
		for key in self.ins.keys():
			if key < position:
				length += len(self.ins[key])
		pos = 0
		for i in self.read.cigartuples:
			if i[0] == 4 and pos <= position:
				length += i[1]
			elif i[0] in [1,5,6]:
				continue
			else:
				pos += i[1]
		return length


# class to process pairs
class ref_pair_align:
	def __init__(self, ref, read1, read2):
		self.ref_seq = ref[int(read1.reference_start): int(read2.reference_end)]
		self.forward = read1
		self.reverse = read2
		self.aligned_forward = read_aligner(self.forward, self.ref_seq, "+")
		self.aligned_reverse = read_aligner(self.reverse, self.ref_seq, "-")

		self.alignment = f">Ref\n{self.ref_seq}\n>{self.forward.query_name}_1\n{self.aligned_forward[0]}\n>{self.reverse.query_name}_2\n{self.aligned_reverse[0]}"
	def merge(self):
		merged_seq = ""
		merged_qual = ""
		pos = 0 
		while pos < len(self.ref_seq):
			if self.aligned_forward[0][pos] == " ":
				if self.aligned_reverse[0][pos] == " ":
					pos += 1
					continue
				elif self.aligned_reverse[0][pos] == "-":
					merged_seq += "-"
					merged_qual += "*"
				else:
					merged_seq += self.aligned_reverse[0][pos]
					merged_qual += self.aligned_reverse[1][pos]
			elif self.aligned_forward[0][pos] == "-":
				if self.aligned_reverse[0][pos] == " " or self.aligned_reverse[0][pos] == "-":
					merged_seq += "-"
					merged_qual += "*"
				else:
					merged_seq += self.aligned_reverse[0][pos]
					merged_qual += self.aligned_reverse[1][pos]
			else:
				if self.aligned_reverse[0][pos] == " " or self.aligned_reverse[0][pos] == "-":
					merged_seq += self.aligned_forward[0][pos]
					merged_qual += self.aligned_forward[1][pos]
				elif self.aligned_reverse[0][pos] == self.aligned_forward[0][pos]:
					merged_seq += self.aligned_forward[0][pos]
					if self.aligned_forward[1][pos] > self.aligned_reverse[1][pos]:
						merged_qual += self.aligned_forward[1][pos]
					else:
						merged_qual += self.aligned_reverse[1][pos]
				else:
					if self.aligned_forward[1][pos] > self.aligned_reverse[1][pos]:
						merged_seq += self.aligned_forward[0][pos]
						merged_qual += self.aligned_forward[1][pos]
					else:
						merged_seq += self.aligned_reverse[0][pos]
						merged_qual += self.aligned_reverse[1][pos]
			pos += 1
		return (merged_seq,merged_qual)
	def if_error(self):
		pos = 0 
		self.error_statistics = []
		ins_forward = ins_info(self.forward)
		ins_reverse = ins_info(self.reverse)
		forward_error = 0
		reverse_error = 0
		while pos < len(self.ref_seq):
			add_len_forward = ins_forward.len_before_position(pos)
			add_len_reverse = ins_reverse.len_before_position(len(self.ref_seq)-1-pos)
			if self.aligned_forward[0][pos] != " " and self.aligned_reverse[0][pos] != " " and self.aligned_forward[0][pos] != self.aligned_reverse[0][pos]:
				if self.aligned_reverse[0][pos] == self.ref_seq[pos]:
					forward_error += 1
					self.error_statistics.append(error_record(f"{pos+add_len_forward}+", self.aligned_forward[0][(pos-5):(pos)], f"{self.aligned_reverse[0][pos]}-{self.aligned_forward[0][pos]}", self.ref_seq[pos], self.forward.reference_start+pos))
				elif self.aligned_forward[0][pos] == self.ref_seq[pos]:
					reverse_error += 1
					self.error_statistics.append(error_record(f"{len(self.ref_seq)-1-pos+add_len_reverse}-", resverse_seq(self.aligned_reverse[0][(pos+1):(pos+6)]), f"{resverse_seq(self.aligned_forward[0][pos])}-{resverse_seq(self.aligned_reverse[0][pos])}", self.ref_seq[pos], self.forward.reference_start+pos))
				else:
					forward_error += 1
					reverse_error += 1
					self.error_statistics.append(error_record(f"{pos+add_len_forward}+", self.aligned_forward[0][(pos-5):(pos)], f"{self.ref_seq[pos]}-{self.aligned_forward[0][pos]}", self.ref_seq[pos], self.forward.reference_start+pos))
					self.error_statistics.append(error_record(f"{len(self.ref_seq)-1-pos+add_len_reverse}-", resverse_seq(self.aligned_reverse[0][(pos+1):(pos+6)]), f"{resverse_seq(self.ref_seq[pos])}-{resverse_seq(self.aligned_reverse[0][pos])}", self.ref_seq[pos], self.forward.reference_start+pos))
			pos += 1
		self.forward_error_count = forward_error
		self.reverse_error_count = reverse_error
		if len(self.error_statistics) == 0:
			return False
		else:
			return True


# process the mapped paired reads and write the merged output and find error
if args.output is None:
	output_file = f"{filename}.merged.fastq"
else:
	output_file = args.output
output_write = open(output_file, 'wt',newline='')
read_count = 0
pair_alignment = []
all_error = []
overlap_size = 0
overlap_count = 0
overlap_list = []
error_read_count = 0
position_overlap = {}

while read_count < len(mapped_paired_list) - 1:
	if mapped_paired_list[read_count].query_name ==  mapped_paired_list[read_count + 1].query_name:
		if_overlap_result = if_overlap(mapped_paired_list[read_count], mapped_paired_list[read_count + 1])
		if if_overlap_result is False:
			read_count += 2
			continue
		elif if_overlap_result[0] == "+":
			overlap_len = int(if_overlap_result[1:])
			overlap_list.append(overlap_len)
			overlap_size += overlap_len
			overlap_count += 1
			for overlap_pos in range(mapped_paired_list[read_count + 1].reference_start,(mapped_paired_list[read_count].reference_end + 1)):
				position_overlap[overlap_pos] = position_overlap.get(overlap_pos,0) + 1
			pair = ref_pair_align(refseq, mapped_paired_list[read_count], mapped_paired_list[read_count + 1])
			pair_alignment.append(pair)
			output_write.write(f"@{mapped_paired_list[read_count].query_name}\n{pair.merge()[0]}\n+\n{pair.merge()[1]}\n")	
			if pair.if_error() is True:
				#print(f">Ref\n{pair.ref_seq}\n>{pair.forward.query_name}_1\n{pair.forward.get_forward_sequence()}\n>{pair.reverse.query_name}_2\n{pair.reverse.get_forward_sequence()}\n#aligned\t{pair.forward.cigartuples}\t{pair.reverse.cigartuples}\n{pair.alignment}\n\n")
				all_error += pair.error_statistics
				if pair.forward_error_count != 0 and pair.reverse_error_count != 0:
					error_read_count += 2
				else:
					error_read_count += 1
			read_count += 2
		elif if_overlap_result[0] == "-":
			overlap_len = int(if_overlap_result[1:])
			overlap_list.append(overlap_len)
			overlap_size += overlap_len
			overlap_count += 1
			for overlap_pos in range(mapped_paired_list[read_count].reference_start,(mapped_paired_list[read_count + 1].reference_end + 1)):
				position_overlap[overlap_pos] = position_overlap.get(overlap_pos,0) + 1
			pair = ref_pair_align(refseq, mapped_paired_list[read_count + 1], mapped_paired_list[read_count])
			pair_alignment.append(pair)
			output_write.write(f"@{mapped_paired_list[read_count].query_name}\n{pair.merge()[0]}\n+\n{pair.merge()[1]}\n")
			if pair.if_error() is True:
				#print(f">Ref\n{pair.ref_seq}\n>{pair.forward.query_name}_1\n{pair.forward.get_forward_sequence()}\n>{pair.reverse.query_name}_2\n{pair.reverse.get_forward_sequence()}\n#aligned\t{pair.forward.cigartuples}\t{pair.reverse.cigartuples}\n{pair.alignment}\n\n")
				all_error += pair.error_statistics
				if pair.forward_error_count != 0 and pair.reverse_error_count != 0:
					error_read_count += 2
				else:
					error_read_count += 1
			read_count += 2
		else:
			errorLogger.error("Read processing error")
			raise SystemExit(1)
	else:
		read_count += 1

output_write.close()
infoLogger.info(f'The output fastq file {output_file} is done')

# statistical analysis of error
ave_overlap_size = overlap_size/overlap_count

# overlap statistic
overlap_write = open(f"{filename}_overlap_statistic.tsv","wt", newline='')
overlap_write.write("Overlap_length\tCount\n")
overlap_dict = {}
for key in overlap_list:
	overlap_dict[key] = overlap_dict.get(key,0) + 1
sorted_keys = sorted(overlap_dict.keys())
for key in sorted_keys:
	overlap_write.write(f"{key}\t{overlap_dict[key]}\n")
overlap_write.close()
	

error_table = open(f"{filename}_error_info.tsv", 'wt',newline='')
upstream_fatsa = open(f"{filename}_upstream.fasta","w", newline='')

total_error_count = len(all_error)
forward_error_count = 0
reverse_error_count = 0
error_type_count = {}

for item in all_error:
	if item.error_read_pos[-1] == "+":
		forward_error_count += 1
	elif item.error_read_pos[-1] == "-":
		reverse_error_count += 1
	if item.error_type in error_type_count.keys():
		error_type_count[item.error_type] += 1
	else:
		error_type_count[item.error_type] = 1

error_table.write(f"#System error information of {args.input_bamfile}\n#Average overlap size: {ave_overlap_size}\n#Error read count: {error_read_count}\n#Total error count: {total_error_count}\n#Forward error: {forward_error_count}\n#Reverse error: {reverse_error_count}\n")
for key in error_type_count.keys():
	error_table.write(f"#{key} mutation count: {error_type_count[key]}\n")

error_table.write("Position_in_reference\tRead_position\tReference_nucleotide_base\tDirection\tMutation\tUpstream\n")

upstream_seq_count = 1
for item in all_error:
	error_table.write(f"{item.error_pos}\t{item.error_read_pos[:-1]}\t{item.ref_base}\t{item.error_read_pos[-1]}\t{item.error_type}\t{item.error_upstream}\n")
	upstream_fatsa.write(f">{upstream_seq_count}\n{item.error_upstream}\n")
	upstream_seq_count += 1

error_table.close()
upstream_fatsa.close()

#draw weblogo plot of upstream
try:
	upstream_data = open(f"{filename}_upstream.fasta","r")
	upstream_seqs = read_seq_data(upstream_data)
	upstream_seqs_data = LogoData.from_seqs(upstream_seqs)
	upstream_options = LogoOptions()
	upstream_options.fineprint = "Upstream of error"
	upstream_format = LogoFormat(upstream_seqs_data, upstream_options)
	eps = eps_formatter(upstream_seqs_data, upstream_format)
	output_upstream_plot = open(f"{filename}_upstream_plot.eps",'wb')
	output_upstream_plot.write(eps)
	output_upstream_plot.close()
	os.system(f"rm {filename}_upstream.fasta")
except FileNotFoundError:
	errorLogger.error('Fail to draw Upstream Weblogo plot. Please double check.')


# funtion to calculate error count in a position
def pos_error(pos, all_error, coverage, pos_overlap_count):
	total_count = 0
	A_mutation = 0
	T_mutation = 0
	G_mutation = 0
	C_mutation = 0
	deletion = 0
	ref_base = refseq[pos]
	for item in all_error:
		if item.error_pos == pos:
			total_count += 1
			ref_base = item.ref_base
			if item.error_type[2] == "A":
				A_mutation += 1
			elif item.error_type[2] == "T":
				T_mutation += 1
			elif item.error_type[2] == "G":
				G_mutation += 1
			elif item.error_type[2] == "C":
				C_mutation += 1
			elif item.error_type[2] == "-":
				deletion += 1
			else:
				errorLogger.error("Error statistical calculation error")
				raise SystemExit(1)
	if int(coverage) != 0:
		frequency = float(total_count/int(coverage))
	else:
		frequency = 0
	return f"{pos}\t{ref_base}\t{A_mutation}\t{T_mutation}\t{G_mutation}\t{C_mutation}\t{deletion}\t{total_count}\t{coverage}\t{frequency}\t{pos_overlap_count}\n"
os.system(f"samtools depth -aa -d0 {args.input_bamfile} > depth.txt")
with open("depth.txt", "r") as depth:
	depth_table = depth.readlines()

error_stat = open(f"{filename}_error_statistic.tsv", 'wt', newline="")
error_stat.write("Reference_position\tReference_base\tA\tT\tG\tC\tDeletion\tTotal_mutation\tCoverage\tMutation_frequency\tOverlap_count\n")

position = 0
for base in refseq:
	#print(depth_table[position].split("\t")[1],position+1)
	if int(depth_table[position].split("\t")[1]) == position+1:
		coverage = depth_table[position].split("\t")[2].strip()
	else:
		errorLogger.error("Error statistical position infomation error")
		raise SystemExit(1)
	try:
		pos_overlap_count = position_overlap[position+1]
	except KeyError:
		pos_overlap_count = 0
	base_info = pos_error(position, all_error, coverage, pos_overlap_count)
	#print(base_info,base)
	if base_info.split("\t")[1] != base:
		errorLogger.error("Error statistical base infomation error")
		raise SystemExit(1)
	error_stat.write(base_info)
	position += 1
error_stat.close()
os.system("rm depth.txt")

#Write plotting R script
r_write = open(f"{filename}_plotting.r", 'wt', newline="")
r_write.write(f"library(ggplot2)\nlibrary(RColorBrewer)\nlibrary(reshape2)\ndatafile = read.delim('{filename}_error_info.tsv', comment.char = '#')\ndatafile = data.frame(as.matrix(datafile))\ndatafile$Position_in_reference = as.numeric(datafile$Position_in_reference)\ndatafile$Read_position = as.numeric(datafile$Read_position)\ndensity_plot = ggplot(data = datafile,aes(x=Read_position)) + geom_histogram(bins = 200, color='black', fill='lightblue') + geom_vline(xintercept= 0, linetype='dashed', color = 'black', size=0.6) + geom_hline(yintercept= 0, linetype='dashed', color = 'black', size=0.6) + theme(legend.position = 'none') + labs(y='Error count')\nupstream = as.matrix(datafile$Upstream)\nout = strsplit(as.character(upstream),'')\nupstream_table = do.call(rbind,out)\nupstream_table = data.frame(upstream_table)\nnames(upstream_table) = c('5','4','3','2','1')\nupstream_table.m = melt(upstream_table,measure.vars=c('5','4','3','2','1'), variable.names='position', value.name='nucl')\nupstream_table.m = data.frame(upstream_table.m)\nnames(upstream_table.m) = c('position','nucleotide')\nupstream_bar_plot = ggplot(upstream_table.m, aes(x=position,fill=nucleotide)) + geom_bar(position = 'stack') + scale_fill_brewer(palette = 'Blues') + labs(x='Base position before error')\nupstream_table$`2_bases_before` = paste(upstream_table$`2`, upstream_table$`1`, sep='')\nupstream_table_sub = subset(upstream_table, grepl('-', `2_bases_before`)!=TRUE)\nupstream_2 = ggplot(data=upstream_table_sub, aes(x=`2_bases_before`)) + geom_bar(color='black', fill='lightblue') + labs(x='2 bases before error', y='Error count')\nmutation_type = datafile$Mutation\nmutation_table = data.frame(as.matrix(table(mutation_type)))\nmutation_table$type = row.names(mutation_table)\nmutation_out = strsplit(as.character(mutation_table$type),'')\nmutation_type_table = do.call(rbind, mutation_out)\nmutation_type_table = data.frame(mutation_type_table[,-2])\nnames(mutation_table) = c('count','type')\nmutation_type_table$count = mutation_table$count\nmutation_heatmap = ggplot(mutation_type_table, aes(x=X1, y=X2)) + geom_tile(alpha=0.5, aes(fill=count)) + labs(x='Reference base', y='Mutated base') + theme(panel.grid = element_blank()) + scale_fill_gradient(low = 'yellow', high = 'red')\noverlap_file = read.table('{filename}_overlap_statistic.tsv',header=1)\noverlap_file$Count = as.numeric(overlap_file$Count)\noverlap_file$Frequency = (overlap_file$Count)/sum(overlap_file$Count)\noverlap_dot_plot = ggplot(overlap_file, aes(x=Overlap_length, y=Frequency)) + geom_point(shape = 21,size=2)\npng(filename = '{filename}_position_error_count.png', width=800,height=800)\ndensity_plot\ndev.off()\npng(filename = '{filename}_upstream_barplot.png', width=800,height=800)\nupstream_bar_plot\ndev.off()\npng(filename = '{filename}_upstream_2.png', width=800,height=800)\nupstream_2\ndev.off()\npng(filename = '{filename}_mutation_heatmap.png', width=800,height=800)\nmutation_heatmap\ndev.off()\npng(filename = '{filename}_overlap_dotplot.png', width=800,height=800)\noverlap_dot_plot\ndev.off()\n")
r_write.close()
try:
	os.system(f"R --vanilla --slave < {filename}_plotting.r")
except FileNotFoundError:
	errorLogger.error('Fail to draw plots by R. Please double check.')


