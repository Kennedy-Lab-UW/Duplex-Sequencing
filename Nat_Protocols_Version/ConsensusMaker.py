#!/usr/bin/env python
"""
Consensus Maker
Version 2.0
By Brendan Kohrn and Scott Kennedy(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195
Based on work by Scott Kennedy
January 21, 2014

Written for Python 2.7.3
Required modules: Pysam, Samtools

Inputs: 
	A position-sorted paired-end BAM file containing reads with a duplex tag in the header.

Outputs:
	1: A paired-end BAM file containing SSCSs
	2: A single-end BAM file containing unpaired SSCSs (if --read_type is 'd')
	3: A single-end BAM file containing reads with less common cigar strings
	4: A single-end BAM file containing reads not in --read_type
	5: A tagcounts file

	Note that quality scores in outputs 1, 2, and 3 are just space fillers and do not signify anything about the
	quality of the sequence.

The program starts at the position of the first good read, determined by the type of read specified on startup.  
It then goes through the file until it finds a new position, saving all reads as it goes.  When it finds a new position, 
it sends the saved reads to the consensus maker, one tag at a time, untill it runs out of tags.  Consensus sequences are 
saved until their mates come up, at which point both are written to the output BAM file, read 1 first.  After making 
consensuses with the reads from the first position, it continues on through the origional file until it finds another 
new position, sends those reads to the consensus maker, and so on until the end of the file.  At the end of the file, 
any remaining reads are sent through the consensus maker, and any unpaired consensuses are written to a file ending in 
UP.bam.  

In the future, this program may be able to autodetect read length.  

usage: consensus_maker.py [-h] [--infile INFILE] [--tag_file tag_file] [--tag_stats tag_stats]
						[--outfile OUTFILE] [--rep_filt REP_FILT]
						[--minmem MINMEM] [--maxmem MAXMEM] [--cut_off cut_off]
						[--Ncut_off Ncut_off] [--read_length READ_LENGTH]
						[--read_type READ_TYPE] [--isize ISIZE]
						[--read_out ROUT] [--filt FILT] [--sam_tag SAM_TAG]

optional arguments:
	-h, --help            show this help message and exit
	--infile INFILE       input BAM file
	--tag_file tag_file     output tagcounts file
	--tag_stats tag_stats output tag_stats file
	--outfile OUTFILE     output BAM file
	--rep_filt REP_FILT   Remove tags with homomeric runs of nucleotides of
						length x. [9]
	--minmem MINMEM     Minimum number of reads allowed to comprise a
						consensus. [3]
	--maxmem MAXMEM     Maximum number of reads allowed to comprise a
						consensus. [1000]
	--cut_off cut_off   Percentage of nucleotides at a given position in a
						read that must be identical in order for a consensus
						to be called at that position. [0.7]
	--Ncut_off Ncut_off With --filt 'n', maximum fraction of Ns allowed in a
						consensus [1.0]
	--read_length READ_LENGTH
						Length of the input read that is being used. [84]
	--read_type READ_TYPE
						A string specifying which types of read to consider.
						Read types: n: Neither read 1 or read 2 mapped. m:
						Either read 1 or read 2 mapped, but not both. p: Both
						read 1 and read 2 mapped, not a proper pair. d: Both
						read 1 and read 2 mapped, proper pair. s: Single
						ended reads ['dpm']
	--isize ISIZE       maximum distance between read pairs
	--read_out ROUT     How often you want to be told what the program is
						doing. [1000000]
	--filt FILT         A string indicating which filters should be
						implemented. Filters: s: soft_clipping filter. o:
						Overlap filter. n: N filter. ['osn']
	--sam_tag SAM_TAG     The SAM tag that store the duplex tag sequence (can
						be set one more times).  Otherwise use the sequence
						in the read name."

Details of different arguments:
	--minmem and --maxmem set the range of family sizes (constrained by cigar score) that can be used to make a
	consensus sequence.  Examples use --minmem of 3 and --maxmem of 1000
		Example 1:
			Ten reads (read_length = 80) have a particular barcode.  Of these ten, nine of them have a cigar string of
			80M, while one has a cigar string of 39M1I40M.  Only the nine with a cigar string of 80M are sent on to be
			made into a SSCS.
		Example 2:
			Three reads (read_length 80) have a particular barcode.  Of these, two have a cigar string of 80M, and one
			has a cigar string of 20M1D60M.  No SSCS results.
		Example 3:
			A family with over 1000 members exists.  A random sample of 1000 reads from that family is used to make a
			SSCS.
	--cut_off sets the strictness of the consensus making.
		Example (--cut_off = 0.7):
			Four reads (read_length = 10) are as follows:
				Read 1: ACTGATACTT
				Read 2: ACTGAAACCT
				Read 3: ACTGATACCT
				Read 4: ACTGATACTT
			The resulting SSCS is:
				ACTGATACNT
	--Ncut_off, with --filt n enabled, sets the maximum percentage of Ns allowed in a SSCS.
		Example (--Ncut_off = .1, --read_length = 20):
			Two SSCSs are generated as follows:
				SSCS 1: ACGTGANCTAGTNCTNTACC
				SSCS 2: GATCTAGTNCATGACCGATA
			SSCS 2 passes the n filter (10%) with 1/20 = 5% Ns, while SSCS 1 does not with 3/20 = 15% Ns.
	--read_length sets the length of the reads imputed.  If this value is set incorrectly, the program will often crash
	with an error message about sequence length not matching quality score length, or will output an empty SSCS bam
	file.
	--read_type sets which reads are considered to have 'good' flags.  Options are:
		d:  Paired-end reads where both reads in the pair map, and where the two are properly paired (read 2 maps in
		the opposite direction and on the opposite strand from read 1).  Flags are 99, 83, 163, and 147  .
		p: Paired-end reads where both reads in the pair map, but the two are not properly paired.  Flags are 97, 81, 
		161, 145, 129, 65, 177, and 113.
		m: Paired-end reads where only one read in the pair maps.  Flags are 181, 117, 137, 133, 73, 89, 69, and 153.
		n: Paired-end reads where neither read in the pair maps, and single end unmapped reads.  Flags are 141, 77, 
		and 4.  
		s: Single end mapped reads.  Flags are 0 and 16.  
	--filt sets which filters are used.  Options are:
		o: Overlap filter. Filters out any read pairs which overlap.  Only works on  reads of type d (see above).
		s: soft_clipping filter.  Filters out any reads which have been soft-clipped in alignment.  This avoids later 
		problems with hard-clipping.  
		n: N filter. Filters out consensus sequences with a higher percentage of Ns than the threshold imposed by 
		--Ncut_off.  Without this option, --Ncut_off doesn't do anything.  
	--isize
		If not -1, sets the maximum distance between read 1 and read 2 for the two to not be considered unpaired.  Only 
		works if --read_type is 'd'
"""

import sys
import pysam
import random
from collections import defaultdict
from argparse import ArgumentParser


def print_read(read_in):
	sys.stderr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read_in.qname, read_in.flag, read_in.tid, 
																			read_in.pos, read_in.mapq, read_in.cigar, 
																			read_in.mrnm, read_in.mpos, read_in.isize, 
																			read_in.seq, read_in.qual, read_in.tags))


def consensus_maker(grouped_reads_list,  cut_off,  read_length):
	# The consensus maker uses a simple "majority rules" algorithm to qmake a consensus at each base position.  If no 
	# nucleotide majority reaches above the minimum theshold (--cut_off), the position is considered undefined and an 'N' 
	# is placed at that position in the read.'''
	nuc_identity_list = [0, 0, 0, 0, 0, 0]  # In the order of T, C, G, A, N, Total
	nuc_key_dict = {0: 'T', 1: 'C', 2: 'G', 3: 'A', 4: 'N'}
	consensus_read = ''

	for i in xrange(read_length):  # Count the types of nucleotides at a position in a read. i is the nucleotide index 
		# within a read in grouped_reads_list
		for j in xrange(len(grouped_reads_list)):  # Do this for every read that comprises a SMI group. j is the read 
			# index within grouped_reads_list
			try:
				if grouped_reads_list[j][i] == 'T':
					nuc_identity_list[0] += 1
				elif grouped_reads_list[j][i] == 'C':
					nuc_identity_list[1] += 1
				elif grouped_reads_list[j][i] == 'G':
					nuc_identity_list[2] += 1
				elif grouped_reads_list[j][i] == 'A':
					nuc_identity_list[3] += 1
				elif grouped_reads_list[j][i] == 'N':
					nuc_identity_list[4] += 1
				else:
					nuc_identity_list[4] += 1
				nuc_identity_list[5] += 1
			except:
				break
		try:
			for j in [0, 1, 2, 3, 4]:
				if float(nuc_identity_list[j])/float(nuc_identity_list[5]) > cut_off:
					consensus_read += nuc_key_dict[j]
					break
				elif j == 4:
					consensus_read += 'N'
		except:
			consensus_read += 'N'
		nuc_identity_list = [0, 0, 0, 0, 0, 0]  # Reset for the next nucleotide position
	return consensus_read, len(grouped_reads_list)


def tag_stats(tag_counts_file, tag_stats_file):
	fam_size_counts = defaultdict(lambda: 0)

	in_file = open(tag_counts_file, 'r')
	out_file = open(tag_stats_file, 'w')
	for line in in_file:
		fam_size_counts[int(line.strip().split()[1].split(":")[0])] += 1
	in_file.close()

	totals = 0
	for size in fam_size_counts.keys():
		fam_size_counts[size] *= int(size)
		totals += int(fam_size_counts[size])

	for size in sorted(fam_size_counts.keys()):
		out_file.write("%s\t%s\n" % (size, float(fam_size_counts[size])/float(totals)))

	out_file.close()
	return True


def main():
	# Parameters to be input.
	parser = ArgumentParser()
	parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
	parser.add_argument("--tag_file",  action="store",  dest="tag_file", help="output tagcounts file",  
						default='sys.stdout', required=True)
	parser.add_argument("--tag_stats",  action="store",  dest="tag_stats", help="output tag_stats file",  
						default='tag_stats.txt', required=True)
	parser.add_argument("--outfile",  action="store", dest="outfile", help="output BAM file", required=True)
	parser.add_argument("--rep_filt", action="store",  type=int, dest='rep_filt', 
						help="Remove tags with homomeric runs of nucleotides of length x. [9]", default=9 )
	parser.add_argument('--minmem', type=int, default=3, dest='minmem', 
						help="Minimum number of reads allowed to comprise a consensus. [3] ")
	parser.add_argument('--maxmem', type=int, default=1000, dest='maxmem', 
						help="Maximum number of reads allowed to comprise a consensus. [1000]")
	parser.add_argument('--cut_off', type=float, default=.7, dest='cut_off', 
						help="Percentage of nucleotides at a given position in a read that must be identical in order for a consensus to be called at that position. [0.7]")
	parser.add_argument('--Ncut_off', type=float, default=1, dest='Ncut_off', 
						help="With --filt 'n', maximum fraction of Ns allowed in a consensus [1.0]")
	parser.add_argument('--read_length', type=int, default=84, dest='read_length',
						help="Length of the input read that is being used. [80]")
	parser.add_argument('--read_type', type=str, action="store", dest='read_type', default="dpm", 
						help="A string specifying which types of read to consider.  Read types: n: Neither read 1 or "
							"read 2 mapped.  m: Either read 1 or read 2 mapped, but not both.  p: Both read 1 and "
							"read 2 mapped, not a propper pair.  d: Both read 1 and read 2 mapped, propper pair.  "
							"s: Single ended reads\n\t\t['dpm']")
	parser.add_argument('--isize', type=int, default=-1, dest='isize', help="maximum distance between read pairs")
	parser.add_argument('--read_out', type=int, default=1000000, dest='rOut',
						help='How often you want to be told what the program is doing. [1000000]')
	parser.add_argument('--filt', type=str, default='osn', dest='filt', 
						help="A string indicating which filters should be implemented.  Filters: s: soft_clipping filter.  o: Overlap filter.  n: N filter.  ['osn']")
	parser.add_argument('--sam_tag', action='append', type=str, dest='samtags', 
						help="The SAM tag that store the duplex tag sequence (can be set one more times). "
							" Otherwise use the sequence in the read name.", default=list())
	o = parser.parse_args()

	# Initialization of all global variables, main input/output files, and main iterator and dictionaries.
	good_flag = []
	if 'd' in o.read_type:
		good_flag.extend((99, 83, 163, 147))
	if 'm' in o.read_type:
		good_flag.extend((181, 117, 137, 133, 73, 89, 69, 153))
	if 'p' in o.read_type:
		good_flag.extend((97, 81, 161, 145, 129, 65, 177, 113))
	if 'n' in o.read_type:
		good_flag.extend((141, 77, 4))
	if 's' in o.read_type:
		good_flag.extend((0, 16))
	if 'u' in o.read_type:
		good_flag.extend((103, 167))

	in_bam_file = pysam.Samfile(o.infile, "rb")  # Open the input BAM file
	out_bam_file = pysam.Samfile(o.outfile, "wb", template=in_bam_file)  # Open the output BAM file
	outNC1 = pysam.Samfile(o.outfile.replace(".bam", "_LCC.bam"), "wb", template=in_bam_file)
	nonmapped_file = pysam.Samfile(o.outfile.replace(".bam", "_NM.bam"), "wb", template=in_bam_file)  # File for reads
	# with strange flags
	if o.read_type == 'd':
		extraneous_read_bam = pysam.Samfile(o.outfile.replace(".bam", "_UP.bam"), "wb", template=in_bam_file)
	samtags = o.samtags

	read_number_count = 0
	nM = 0
	bF = 0
	oL = 0
	sC = 0
	rT = 0
	nC = 0

	LCC = 0
	consenuses_made = 0
	if o.read_type == 'd':
		UP = 0

	file_done = False  # Initialize end of file bool
	read_one = False

	quality_score = 'J' * o.read_length  # Set a dummy quality score

	bam_entry = in_bam_file.fetch(until_eof=True)  # Initialize the iterator
	read_window = [bam_entry.next(), '']  # Get the first read
	window_position = 0

	read_dict = {}  # Initialize the read dictionary
	tag_dict = defaultdict(lambda: 0)  # Initialize the tag dictionary

	consensus_dict = {}

# Start going through the input BAM file, one position at a time.
	for line in bam_entry:
		window_position += 1
		read_window[window_position % 2] = line
		# Reinitialize first line
		if read_one is True:
			window_position -= 1
		while (read_window[window_position % 2].pos == read_window[(window_position-1) % 2].pos and
							file_done is False and read_one is False) or read_one is True:
			if read_number_count % o.rOut == 0:
				sys.stderr.write("Reads processed:" + str(read_number_count) + "\n")
			try:
				if 0 < len(samtags):
					tag = "".join([tag_tuple[1] for tag_tuple in read_window[window_position%2].tags if tag_tuple[0] in samtags])
				else:
					tag = read_window[window_position % 2].qname.split('|')[1].split('/')[0]
				tag += (":1" if read_window[window_position % 2].is_read1 is True
							else (":2" if read_window[window_position % 2].is_read2 is True else ":se"))
				tag_dict[tag] += 1
			except:
				print read_number_count
				raise

			# Overlap filter: filters out overlapping reads (with --filt o)
			overlap = False
			if 'o' in o.filt:
				if read_window[window_position % 2].pos < read_window[window_position % 2].mpos and \
								read_window[window_position % 2].mpos < read_window[window_position % 2].pos + \
								o.read_length and int(read_window[window_position % 2].flag) in (83, 99, 147, 163):
					overlap = True
				elif read_window[window_position % 2].pos > read_window[window_position % 2].mpos \
						and read_window[window_position % 2].pos < read_window[window_position % 2].mpos + o.read_length \
						and int(read_window[window_position % 2].flag) in (83, 99, 147, 163):
					overlap = True
				elif read_window[window_position % 2].pos == read_window[window_position % 2].mpos \
						and int(read_window[window_position % 2].flag) in (83, 99, 147, 163):
					overlap = True
					
			read_number_count += 1

			# soft_clip filter: filters out soft_clipped reads (with --filt s)
			soft_clip = False
			if 's' in o.filt:
				if read_window[window_position % 2].cigar is not None:
					for tupple in read_window[window_position % 2].cigar:
						if tupple[0] == 4:
							soft_clip = True

			# Check if the given read is good data
			if int(read_window[window_position % 2].flag) in good_flag and overlap is False and soft_clip is False:
				if ('A' * o.rep_filt in tag) or ('C' * o.rep_filt in tag) or ('G' * o.rep_filt in tag) \
						or ('T' * o.rep_filt in tag):
					# Check for bad barcodes
					nM += 1
					nonmapped_file.write(read_window[window_position % 2])
					rT += 1
				else:
					# Add the sequence to the read dictionary
					if tag not in read_dict:
						read_dict[tag] = [read_window[window_position % 2].flag, read_window[window_position % 2].rname,
										read_window[window_position % 2].pos, read_window[window_position % 2].mrnm,
										read_window[window_position % 2].mpos, read_window[window_position % 2].isize,
										{str(read_window[window_position % 2].cigar):[0, read_window[window_position%2].cigar]}]

					if str(read_window[window_position % 2].cigar) not in read_dict[tag][6]:
						read_dict[tag][6][str(read_window[window_position % 2].cigar)] = \
							[0, read_window[window_position % 2].cigar]

					read_dict[tag][6][str(read_window[window_position % 2].cigar)].append(read_window[window_position%2].seq)
					read_dict[tag][6][str(read_window[window_position % 2].cigar)][0] += 1
			else:
				nM += 1
				nonmapped_file.write(read_window[window_position % 2])
				if int(read_window[window_position % 2].flag) not in good_flag:
					bF += 1
				elif overlap is True:
					oL += 1
				elif soft_clip is True:
					sC += 1

			window_position += 1
			if read_one is False:
				try:  # Keep StopIteration error from happening at the end of a file
					read_window[window_position % 2] = bam_entry.next()  # Iterate the line
				except:
					file_done = True  # Tell the program that it has reached the end of the file
			else:
				read_one = False
		else:

			# Send reads to consensus_maker
			read_one = True
			for dict_tag in read_dict.keys():  # Extract sequences to send to the consensus maker
				# Cigar string filtering
				cigar_string_set = {}

				for cigar_string in read_dict[dict_tag][6].keys():  # Determine the most common cigar string
					cigar_string_set[cigar_string] = read_dict[dict_tag][6][cigar_string][0]

				max_cigar = max(cigar_string_set)

				if cigar_string_set[max_cigar] >= o.minmem:
					if cigar_string_set[max_cigar] <= o.maxmem:
						consenuses_made += 1
						consensus, fam_size = consensus_maker(read_dict[dict_tag][6][max_cigar][2:], o.cut_off,
																o.read_length)
					else:
						consenuses_made += 1
						consensus, fam_size = consensus_maker(random.sample(read_dict[dict_tag][6][max_cigar][2:],
																			o.maxmem), o.cut_off, o.read_length)

					for cigar_string in read_dict[dict_tag][6].keys():
						if cigar_string != max_cigar:
							for n in xrange(2, len(read_dict[dict_tag][6][cigar_string][2:])):
								a = pysam.AlignedRead()
								a.qname = dict_tag + ':' + str(fam_size)
								a.flag = read_dict[dict_tag][0]
								a.seq = read_dict[dict_tag][6][cigar_string][n]
								a.rname = read_dict[dict_tag][1]
								a.pos = read_dict[dict_tag][2]
								a.mapq = 255
								a.cigar = read_dict[dict_tag][6][cigar_string][1]
								a.mrnm = read_dict[dict_tag][3]
								a.mpos = read_dict[dict_tag][4]
								a.isize = read_dict[dict_tag][5]
								a.qual = quality_score  
								outNC1.write(a)
								LCC += 1

					# Filter out consensuses with too many Ns in them
					if (consensus.count("N")/ float(len(consensus)) <= o.Ncut_off and 'n' in o.filt) \
							or ('n' not in o.filt):
						# Write a line to the consensus_dictionary
						a = pysam.AlignedRead()
						a.qname = dict_tag + ":" + str(fam_size)
						a.flag = read_dict[dict_tag][0]
						a.seq = consensus
						a.rname = read_dict[dict_tag][1]
						a.pos = read_dict[dict_tag][2]
						a.mapq = 255
						a.cigar = read_dict[dict_tag][6][max_cigar][1]
						a.mrnm = read_dict[dict_tag][3]
						a.mpos = read_dict[dict_tag][4]
						a.isize = read_dict[dict_tag][5]
						a.qual = quality_score

						# Write SSCSs to output BAM file in read pairs.
						altTag = dict_tag.replace(("1" if "1" in dict_tag else "2"), ("2" if "1" in dict_tag else "1"))

						if altTag in consensus_dict:
							if a.is_read1 is True:
								out_bam_file.write(a)
								out_bam_file.write(consensus_dict.pop(altTag))
							else:
								out_bam_file.write(consensus_dict.pop(altTag))
								out_bam_file.write(a)
						else:
							consensus_dict[dict_tag] = a
					else:
						nC += 1
		read_dict = {}  # Reset the read dictionary
		if o.read_type == 'd':
			if o.isize != -1:
				for consensus_tag in consensus_dict.keys():
					if consensus_dict[consensus_tag].pos + o.isize < read_window[window_position % 2].pos:
						extraneous_read_bam.write(consensus_dict.pop(consensus_tag))
						UP += 1

	# Write unpaired SSCSs
	for consensus_tag in consensus_dict.keys():
		if o.read_type == 'd':
			extraneous_read_bam.write(consensus_dict.pop(consensus_tag))
			UP += 1
		else:
			out_bam_file.write(consensus_dict.pop(consensus_tag))

	# Close BAM files
	in_bam_file.close()
	out_bam_file.close()
	nonmapped_file.close()
	outNC1.close()

	if o.read_type == 'd':
		extraneous_read_bam.close()

	# Write summary statistics
	sys.stderr.write("Summary Statistics: \n")
	sys.stderr.write("Reads processed:" + str(read_number_count) + "\n")
	sys.stderr.write("Bad reads: %s\n" % nM)
	sys.stderr.write("\tReads with Bad Flags: %s\n" % bF)
	sys.stderr.write("\tOverlapping Reads: %s\n" % oL)
	sys.stderr.write("\tsoft_clipped Reads: %s\n" % sC)
	sys.stderr.write("\tRepetitive Duplex Tag: %s\n" % rT)
	sys.stderr.write("Reads with Less Common Cigar Strings: %s\n" % LCC)
	sys.stderr.write("Consensuses Made: %s\n" % consenuses_made)
	sys.stderr.write("Consensuses with Too Many Ns: %s\n\n" % nC)

	# Write the tag counts file.
	tag_file = open( o.tag_file, "w" )
	tag_file.write ( "\n".join(["%s\t%d" % (SMI, tag_dict[SMI]) 
								for SMI in sorted(tag_dict.keys(), key=lambda x: tag_dict[x], reverse=True ) ] ))
	tag_file.close()
	tag_stats(o.tag_file, o.tag_stats)

if __name__ == "__main__":
	main()
