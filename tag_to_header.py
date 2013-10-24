'''
Tag To Header
Version 2.0
By Joe Hiatt, Scott Kennedy, Brendan Kohrn and Mike Schmitt
October 23, 2013
'''

import sys
from argparse import ArgumentParser
from Bio import SeqIO


class fastQRead:
    def __init__(self, in1, in2, in3, in4):
        '''This class is meant to hold a single fastQ read.
        '''
        self.name=in1.strip().split("@")[1] if "@" in in1 else in1
        self.seq=in2.strip()
        self.spacer=""
        self.qual=in4.strip()
        if len(self.seq)!=len(self.qual):
            raise ValueError("Sequence and quality scores of different lengths!/n%s/n%s/n%s/n%s" % (in1, in2, "+", in4))

    def __getitem__(self, index):
        '''This should allow slicing of the read to proceed properly.
        '''
        if isinstance(index, int):
            return self.seq[index]
        elif isinstance(index, slice):
            answer = self.__class__(self.name, self.seq[index], self.spacer, self.qual[index])
            return answer
        raise ValueError("Invalid index")


class fastQItterator:
    def __init__(self, inFile):
        '''This class will go through a fastQ file one line at a time.
        '''
        self.source=inFile
        self.eof=False
    
    def next(self):
        new=[]
        for j in xrange(4):
            try:
                tmp=self.souce.next()
            except StopIteration:
                self.eof=True
                return("EOF") 
            new.append(self.source.next())
        newRead=fastQRead(new[0],new[1],new[2],new[3])
        return(newRead)

    def close(self):
        self.source.close()
        return(True)


class fastqWriter:
    def __init__(self, outFile):
        self.file=outFile
        self.firstLine=True
    
    def write(self, read):
        if self.firstLine==True:
            self.file.write("@" + read.name)
            self.firstLine=False
        else:
            self.file.write("\n@" + read.name)
        self.file.write("\n" + read.seq)
        self.file.write("\n" + read.spacer)
        self.file.write("\n" + read.qual)
        return(True)
    
    def close(self):
        self.file.close()
        return(True)


def tagExtractFxn(x, blen):
    '''this is the function that extracts the UID tags from both the 
    forward and reverse read.  Assigns read1 the sequence from some 
    position to the end, then read2 from some position to the end, 
    then assigns tag1 from the 5'-end to length of the UID tag for 
    read1 and then read 2.
    '''
    return(x[0][:blen], x[1][:blen])


def hdrRenameFxn(x, y, z):
    '''this function renames the header with the formatting of 
    *header coordinates,etc*, *tag from read1*, *tag from read2*, 
    *read designation from original header*
    '''
    return("%s#%s%s/%s" % (x.split("#")[0], y, z, x.split("/")[-1]))

def main():
    parser = ArgumentParser()
    parser.add_argument('--adapter1',  default = None,  dest = 'adapterSeq1')
    parser.add_argument('--adapter2',  default = None,  dest = 'adapterSeq2')
    parser.add_argument('--infile1', default = None, dest = 'infile1')
    parser.add_argument('--infile2', default = None, dest = 'infile2')
    parser.add_argument('--outfile1', default = None, dest = 'outfile1')
    parser.add_argument('--outfile2', default = None, dest = 'outfile2')
    parser.add_argument('--barcode_length', type = int, default = 12, dest = 'blength')
    parser.add_argument('--spacer_length', type = int, default = 5, dest = 'slength')
    parser.add_argument('--read_out', type = int, default = 1000000, dest = 'rOut')
    o=parser.parse_args()


    in1=fastQItterator(open(o.infile1, 'rU'))
    in2=fastQItterator(open(o.infile2, 'rU'))
    out1=fastqWriter(open(o.outfile1, 'w'))
    out2=fastqWriter(open(o.outfile2, 'w'))

    ctr=0
    nospacer = 0
    goodreads = 0
    badtag = 0
    oldBad = 0
    isEOF=False

    while isEOF==False:
        read1 = in1.next()
        read2 = in2.next()
        if read1 == "EOF" or read2 == "EOF":
            isEOF = True
        else:
            
            ctr += 1
            if o.adapterSeq1 != None and o.adapterSeq2 != None and (read1.seq[o.blength:o.blength + o.slength] != o.adapterSeq1 or read2[o.blength:o.blength + o.slength] != o.adapterSeq2):
        #        sys.stderr.write('Error: something is wrong with the spacers')
        #        print(read1)
        #        print(read2)
                nospacer += 1

            else :
                #extract tags
                tag1, tag2 = tagExtractFxn((read1.seq, read2.seq),o.blength)
                
                #header reconstruction
                read1.name = hdrRenameFxn(read1.name, tag1, tag2) 
                read2.name = hdrRenameFxn(read2.name, tag1, tag2)
                
                #fastq reconstruction
                if (tag1.isalpha() and tag1.count('N') == 0) and (tag2.isalpha() and tag2.count('N') == 0) :
                    rOut1 = read1[o.blength + o.slength:]
                    rOut2 = read2[o.blength + o.slength:]
                    
                    out1.write(rOut1)
                    out2.write(rOut2)
                    goodreads += 1
                else: 
                    #sys.stderr.write('Error: something is wrong with the tags')
                    #print(read1)
                    #print(read2)
                    badtag += 1
            
                if ctr%o.rOut==0:
                    sys.stderr.write("Total sequences processed: %s\n" % (ctr))
                    sys.stderr.write("Sequences passing filter: %s\n" % (goodreads))
                    sys.stderr.write("Missing spacers: %s\n" % (nospacer))
                    sys.stderr.write("Bad tags: %s\n\n" % (badtag))
                    if badtag == oldBad+o.rOut:
                        raise IOError("Error between lines %s and %s." % ((ctr-o.rOut)*4,(ctr-o.rOut)*4))
                    else:
                        oldBad = badtag

    in1.close()
    in2.close()
    out1.close()
    out2.close()

    sys.stderr.write("Summary statistics:\n")
    sys.stderr.write("Total sequences processed: %s\n" % (ctr))
    sys.stderr.write("Sequences passing filter: %s\n" % (goodreads))
    sys.stderr.write("Missing spacers: %s\n" % (nospacer))
    sys.stderr.write("Bad tags: %s\n\n" % (badtag))

if __name__ == "__main__":
    main()
