#!/usr/bin/env python
'''
DCS Filter
Version 2.0
By Brendan Kohrn, Scott Kennedy(1), and Mike Schmitt(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
Based on work by Scott Kennedy, Mike Schmitt
December 17, 2013

Written for Python 2.7.3
Required modules: Pysam, Samtools, BioPython

Inputs:
    A position-sorted paired-end BAM file containing SSCSs
    
Outputs: 
    1: A paired-end BAM file containing DCSs
    2: A single-end BAM file containing unpaired DCSs
    3: A pair of fastq files containing DCSs for use in realigning.
    
    Note: Quality scores and cigar strings in these files are meaningless. 

This program goes through the input file by position, making DCSs as it goes and writing them to file.  At the end of
the run, any unpaired DCSs are written to a file ending in _UP.bam.

usage: DuplexMaker.py [-h] [--infile INFILE] [--outfile OUTFILE]
                      [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]
                      [--barcode_length BLENGTH] [--read_out ROUT]
                      [--gzip-fqs]

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --outfile OUTFILE     output BAM file
  --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1.0]
  --readlength READ_LENGTH
                        Length of the input read that is being used. [84]
  --barcode_length BLENGTH
                        Length of the duplex tag sequence. Should match the value in tag_to_header.  [12]
  --read_out ROUT       How often you want to be told what the program is
                        doing. [1000000]
  --gzip-fqs            Output gzipped fastqs [False]
'''

import sys
import pysam
import re
import gzip
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import defaultdict
from argparse import ArgumentParser


def print_read(read_in):
	sys.stderr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read_in.qname, read_in.flag, read_in.tid,
																		read_in.pos, read_in.mapq, read_in.cigar,
																		read_in.mrnm, read_in.mpos, read_in.isize,
																		read_in.seq, read_in.qual, read_in.tags))
	return

def dcs_maker(grouped_reads_list,  read_length):
	# The Duplex maker substitutes an N if the two input sequences are not identical at a position.
	consensus_read = ''
	for i in xrange(read_length):  # rebuild consensus read taking into account the cutoff percentage
		if grouped_reads_list[0][i] == grouped_reads_list[1][i]:
			consensus_read += grouped_reads_list[0][i]
		else:
			consensus_read += "N"
	return consensus_read


def fastq_open(outfile, gzip_fastq, end):
	fn = outfile.replace('.bam', '') + "." + end + ".fq"
	if gzip_fastq:
		fn += ".gz"
		return gzip.open(fn, 'wb')
	else:
		return open(fn, 'w')


def main():
	# Parameters to be input.
	parser = ArgumentParser()
	parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
	parser.add_argument("--outfile",  action="store", dest="outfile", help="output BAM file",  required=True)
	parser.add_argument('--Ncutoff', type=float, default=1.0, dest='Ncutoff',
						help="Maximum percentage of Ns allowed in a consensus [1.0]")
	parser.add_argument('--readlength', type=int, default=84, dest='read_length',
						help="Length of the input read that is being used. [84]")
	parser.add_argument('--barcode_length', type=int, default=12, dest='blength',
						help='Length of the duplex tag sequence. Should match the value in tag_to_header.[12]')
	parser.add_argument('--read_out', type=int, default=1000000, dest='rOut',
						help='How often you want to be told what the program is doing. [1000000]')
	parser.add_argument('--gzip-fqs', action="store_true", default=False, dest='gzip_fastqs',
						help='Output gzipped fastqs [False]')
	o = parser.parse_args()

	# Initialization of all global variables, main input/output files, and main iterator and dictionaries.
	in_bam = pysam.Samfile(o.infile, "rb")  # Open the input BAM file
	out_bam = pysam.Samfile(o.outfile, "wb", template=in_bam)  # Open the output BAM file
	fastq_file1 = fastq_open(o.outfile, o.gzip_fastqs, 'r1')
	fastq_file2 = fastq_open(o.outfile, o.gzip_fastqs, 'r2')

	read_num = 0
	duplexes_made = 0
	uP = 0
	nC = 0

	file_done = False  # Initialize end of file bool
	finished = False
	read_one = True

	bam_entry = in_bam.fetch(until_eof=True)  # Initialize the iterator
	first_read = bam_entry.next()  # Get the first read
	read_dict = {}  # Initialize the read dictionary
	first_tag = first_read.qname.split(":")[0]
	qual_score = first_read.qual  # Set a dummy quality score
	consensus_dict = {}
	cig_dum = first_read.cigar  # set a dummy cigar score

	# Start going through the input BAM file, one position at a time.
	for line in bam_entry:
		# Reinitialize first line
		read_num += 1
		if read_one is True and first_read.is_unmapped is False:

			read_dict[first_tag] = [first_read.flag, first_read.rname, first_read.pos, first_read.mrnm,
									first_read.mpos, first_read.isize, first_read.seq]
			read_one = False

		while line.pos == first_read.pos and file_done is False:
			tag = line.qname.split(":")[0]  # Extract the barcode
			# Add the sequence to the read dictionary

			if line.is_unmapped is False:
				read_dict[tag] = [line.flag, line.rname, line.pos, line.mrnm, line.mpos, line.isize, line.seq]
			try:  # Keep StopIteration error from happening
				line = bam_entry.next()  # Iterate the line
				read_num += 1
			except:
				file_done = True  # Tell the program that it has reached the end of the file
				read_num += 1

			if read_num % o.rOut == 0:
				sys.stderr.write("%s reads processed\n" % read_num)
		else:
			# Send reads to dcs_maker
			first_read = line  # Store the present line for the next group of lines
			first_tag = first_read.qname.split(":")[0]
			read_one = True
			dict_keys = read_dict.keys()

			for dict_tag in read_dict.keys():  # Extract sequences to send to the dcs_maker
				switch_tag = dict_tag[o.blength:] + dict_tag[:o.blength]

				try:
					consensus = dcs_maker([read_dict[dict_tag][6], read_dict[switch_tag][6]],  o.read_length)
					duplexes_made += 1
					# Filter out consensuses with too many Ns in them
					if consensus.count("N")/len(consensus) > o.Ncutoff:
						nC += 1
					else:
						# Write a line to the consensus_dictionary
						a = pysam.AlignedRead()
						a.qname = dict_tag
						a.flag = read_dict[dict_tag][0]

						if a.is_reverse is True:
							tmp_seq = Seq(consensus, IUPAC.unambiguous_dna)
							a.seq = str(tmp_seq.reverse_complement())
						else:
							a.seq = consensus

						a.rname = read_dict[dict_tag][1]
						a.pos = read_dict[dict_tag][2]
						a.mapq = 255
						a.cigar = cig_dum
						a.mrnm = read_dict[dict_tag][3]
						a.mpos = read_dict[dict_tag][4]
						a.isize = read_dict[dict_tag][5]
						a.qual = qual_score

			# Write DCSs to output BAM file in read pairs.
						if dict_tag in consensus_dict:
							if a.is_read1 is True:
								fastq_file1.write('@:%s\n%s\n+\n%s\n' % (a.qname, a.seq, a.qual))
								out_bam.write(a)
								fastq_file2.write('@:%s\n%s\n+\n%s\n' % (consensus_dict[dict_tag].qname,
																		consensus_dict[dict_tag].seq,
																		consensus_dict[dict_tag].qual))
								out_bam.write(consensus_dict.pop(dict_tag))
							else:
								fastq_file1.write('@:%s\n%s\n+\n%s\n' % (consensus_dict[dict_tag].qname,
																		consensus_dict[dict_tag].seq,
																		consensus_dict[dict_tag].qual))
								out_bam.write(consensus_dict.pop(dict_tag))
								fastq_file2.write('@:%s\n%s\n+\n%s\n' % (a.qname, a.seq, a.qual))
								out_bam.write(a)
						else:
							consensus_dict[dict_tag] = a

					del read_dict[dict_tag]
					del read_dict[switch_tag]

				except:
					pass

		read_dict = {}  # Reset the read dictionary

	# Close BAM files
	in_bam.close()

	# Write unpaired DCSs
	for consTag in consensus_dict.keys():
		a = pysam.AlignedRead()
		a.qname = consTag
		a.flag = 5
		a.seq = '.' * o.read_length
		a.rname = consensus_dict[consTag].rname
		a.pos = consensus_dict[consTag].pos
		a.mapq = 255
		a.cigar = cig_dum
		a.mrnm = consensus_dict[consTag].mrnm
		a.mpos = consensus_dict[consTag].pos
		a.isize = consensus_dict[consTag].isize
		a.qual = qual_score

		if consensus_dict[consTag].is_read1 is False:
			fastq_file1.write('@:%s\n%s\n+\n%s\n' % (a.qname, a.seq, a.qual))
			out_bam.write(a)
			fastq_file2.write('@:%s\n%s\n+\n%s\n' % (consensus_dict[consTag].qname, consensus_dict[consTag].seq,
													consensus_dict[consTag].qual))
			out_bam.write(consensus_dict.pop(consTag))
		else:
			fastq_file1.write('@:%s\n%s\n+\n%s\n' % (consensus_dict[consTag].qname, consensus_dict[consTag].seq,
													consensus_dict[consTag].qual))
			out_bam.write(consensus_dict.pop(consTag))
			fastq_file2.write('@:%s\n%s\n+\n%s\n' % (a.qname, a.seq, a.qual))
			out_bam.write(a)

		uP += 1

	fastq_file1.close()
	fastq_file2.close()
	out_bam.close()

	# Write summary statistics.  Duplexes made includes unpaired duplexes
	sys.stderr.write("Summary Statistics: \n")
	sys.stderr.write("Reads Processed: %s\n" % read_num)
	sys.stderr.write("Duplexes Made: %s\n" % duplexes_made)
	sys.stderr.write("Unpaired Duplexes: %s\n" % uP)
	sys.stderr.write("N-clipped Duplexes: %s\n" % nC)

if __name__ == "__main__":
	main()
