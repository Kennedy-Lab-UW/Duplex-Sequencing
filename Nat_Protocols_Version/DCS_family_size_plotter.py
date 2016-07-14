import matplotlib.pyplot as plt
import sys
from collections import defaultdict
from argparse import ArgumentParser
import numpy as np


parser=ArgumentParser()
parser.add_argument("--tagfile", action="store", dest="tagfile", required=True)
parser.add_argument("--name", action="store", dest="name", required=True)
o = parser.parse_args()

bam_line = sys.stdin
tag_list = set()
noMT_tag_list = set()
good_fam_size_dict = defaultdict(lambda:0)
failed_fam_size_dict = defaultdict(lambda:0)
total_fam_size_dict = defaultdict(lambda:0)

for line in bam_line:
	if line.strip('\n').split('\t')[2] != "chrM":
		tag_list.add(line.strip('\n').split('\t')[0].split(':')[0])

tag_file = open(o.tagfile, 'r')

total_reads = 0
for line in tag_file:
	splitline = line.strip('\n').split('\t')
	total_reads += int(splitline[1])

	if splitline[0].split(':')[0] in tag_list:

		good_fam_size_dict[int(splitline[1])] += 1
		total_fam_size_dict[int(splitline[1])] += 1
	elif int(splitline[1]) > 2:
		failed_fam_size_dict[int(splitline[1])] += 1
		total_fam_size_dict[int(splitline[1])] += 1

for family_size in total_fam_size_dict.keys():
	total_fam_size_dict[family_size] *= int(family_size)
for family_size in good_fam_size_dict.keys():
	good_fam_size_dict[family_size] *= int(family_size)
for family_size in failed_fam_size_dict.keys():
	failed_fam_size_dict[family_size] *= int(family_size)

tag_file.close()
x_value = range(1,300)
good_y_value = []
failed_y_value = []
total_y_value = []

for family_size in x_value:

	if family_size in good_fam_size_dict.keys():
		good_y_value.append(float(good_fam_size_dict[family_size])/float(total_reads))
	else:
		good_y_value.append(0)

	if family_size in failed_fam_size_dict.keys():
		failed_y_value.append(float(failed_fam_size_dict[family_size])/float(total_reads))
	else:
		failed_y_value.append(0)

	if family_size in total_fam_size_dict.keys():
		total_y_value.append(float(total_fam_size_dict[family_size])/float(total_reads))
	else:
		total_y_value.append(0)

failed_y = np.array(failed_y_value)
good_y = np.array(good_y_value)

#plt.bar(x_value, failed_y, alpha=0.3, bottom=good_y)
plt.bar(x_value, good_y, color='r')
plt.xlabel('Family Size')
plt.ylabel('Proportion of Total Reads')
plt.savefig(o.name + '_combined_fam_size.png', bbox_inches='tight')



#plt.bar(failed_x_value, failed_y_value)
#plt.xlabel('Family Size')
#plt.ylabel('Proportion of Total Reads')
#plt.savefig('failed_reads_fam_size.png', bbox_inches='tight')

#plt.bar(x_value, total_y_value)
#plt.xlabel('Family Size')
#plt.ylabel('Proportion of Total Reads')
#plt.savefig(o.name + '_total_reads_fam_size.png', bbox_inches='tight')

