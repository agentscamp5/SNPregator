#!/usr/bin/env python2
'''
SNPregator Version 1.0.0
Author : Alexander Paul (agentscamp5)
'''
import time
import json
import os
import math
import sys
import decimal
import multiprocessing
import argparse
import scipy.stats
import matplotlib.pyplot as plt


def p_val(cells,numerator):
    '''calculate the hypergeometric probability mass function given a contingency table and the numerator for
    pmf calculation, cells is an array of int values and numerator should be a Decimal object
    '''
    if not isinstance(numerator,decimal.Decimal):
        numerator = decimal.Decimal(numerator)
    cells[0] = math.factorial(cells[0])
    denominator = reduce(lambda x,y : x * math.factorial(y) , cells)
    p_value = numerator / denominator
    #convert Decimal object into float, having numerator as Decimal allows for divison of numerator and denominator
    #to create a very small Decimal value instead of just rounding down to 0L
    return float(p_value)


def recursive_fisher(fixedcells,rowlength,Row1temp,Row2temp,Columnstemp,N,p_cutoff,p_significance,numerator):
    '''recursively compute pmf for all tables with row and column totals equal to table being tested'''
    if not isinstance(numerator,decimal.Decimal):
        numerator = decimal.Decimal(numerator)
    originalRow1temp = Row1temp
    originalRow2temp = Row2temp
    p_total = 0.0
    x_used = []
    #test all possible values for cells in the first variable column
    for x in xrange(originalRow1temp+1):
        x_used.append(x)
        Row1temp = originalRow1temp - x
        y = Columnstemp[0] - x
        Row2temp = originalRow2temp - y
        if (Row1temp > originalRow1temp or Row2temp > originalRow2temp or Row1temp < 0 or Row2temp < 0):
            continue
        #if degrees of freedom k != 1, then we lock one degree of freedom in place and recurse on the new table with k-1 degrees
        #of freedom
        if rowlength != 2:
                p_partial = recursive_fisher(fixedcells + [x,y],rowlength-1,Row1temp,Row2temp,Columnstemp[1:],N,p_cutoff,p_significance,numerator)
                p_total += p_partial
        else:
            #calculate hypergeometric pmf value for this table
            p_partial = p_val(fixedcells+[x,Row1temp,y,Row2temp,N],numerator)
            #if table is at least as unlikely as the original table, add its pmf value to our total
            if p_partial <= p_cutoff:
                p_total += p_partial
            # first tail of distribution explored, move to other tail of distribution
            else:
                break
        if p_total >= p_significance:
            return p_total
    #same as above code, but working from opposite end of hypergeometric distribution
    for x in xrange(originalRow1temp+1,0,-1):
        if x in x_used:
            break
        Row1temp = originalRow1temp - x
        y = Columnstemp[0] - x
        Row2temp = originalRow2temp - y
        if (Row1temp > originalRow1temp or Row2temp > originalRow2temp or Row1temp < 0 or Row2temp < 0):
            continue
        if rowlength != 2:
                p_total += recursive_fisher(fixedcells + [x,y],rowlength-1,Row1temp,Row2temp,Columnstemp[1:],N,p_cutoff,p_significance,numerator)
        else:
            p_partial = p_val(fixedcells+[x,Row1temp,y,Row2temp,N],numerator)
            if p_partial <= p_cutoff:
                p_total += p_partial
            else:
                break
        if p_total >= p_significance:
            return p_total
    return p_total

def twobymfishersexact(row1,row2,p_significance):
    '''calculate p-value for contingency table using Fisher's Exact Test, row1 is the first
    row of values in the contingency table, row2 is the second, and p_significance is a float 
    representing the level of significance required for a table to be accepted as significant'''
    if len(row1) <= 1:
        return 1.0
    R1 = sum(row1)
    R2 = sum(row2)
    N = R1 + R2
    #sort column totals so that smaller columns are chosen first, minimizing the number of recursive calls
    Columns = sorted([row1[x] + row2[x] for x in xrange(len(row1))])
    numerator = decimal.Decimal(reduce(lambda x, y: x * math.factorial(y), [math.factorial(Columns[0])] + Columns[1:] + [R1,R2]))
    p_cutoff = p_val(row1 + row2 + [N],numerator)
    #make recursive call to sum 
    p_total = recursive_fisher([],len(row1),R1,R2,Columns,N,p_cutoff,p_significance,numerator)
    return p_total

def make_parser():
    '''creates an ArgumentParser object that will handle user input from the command line and call the selected tool'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description="Tools to test association between groups in VCF files and provides filters for SNPs")
    subparsers = parser.add_subparsers(help="tool list",title="tools",dest="parser")
    
    #association tool subparser
    association_parser = subparsers.add_parser('association',help="Create contingency tables of aggregate allele counts for a VCF file and perform tests of assocation")
    NUM_PROCS = multiprocessing.cpu_count()
    association_parser.add_argument('-n','--numprocs', type=int,default=NUM_PROCS, help="Number of worker processes to launch, default is number of CPUs in computer")
    association_parser.add_argument("vcf",metavar="VCF",help="input vcf file")
    association_parser.add_argument("samples", metavar="SAMPLES",help="file with sample names")
    association_parser.add_argument("-a","--acceptlow",action='store_false',help="if set, chi-square test will not be performed if table has cell with value < 5")
    association_parser.add_argument("-o","--output",default="out.vcf",help="output file name, default is 'out.vcf'")
    association_parser.add_argument("-b","--bonferroni",action='store_true',help="if set, apply the bonferroni correction for multiple comparisons")
    association_parser.add_argument("-q","--condense",choices=[1,2,3],default=1,help="Condense SNPs to [1: raw input(unchanged); 2: ref-ref,alt-* ;3: ref-ref,ref-alt,alt-alt], default is 1")
    association_parser.add_argument("-f","--fisher",action='store_true',help="if set, use Fisher Exact Test metric not Pearson Chi-Square Test")
    association_parser.add_argument("-g","--groupsize",type=int,default=None,help="sample size for group 1, default is half of total sample size")
    association_parser.add_argument("-p","--pvalue",type=float,default=0.05,help="cutoff pvalue for significance, default is 0.05")
    association_parser.add_argument("-d","--addheader",action='store_true',help="if set, new header lines describing association test metrics and group names are added")
    association_parser.add_argument("-k","--keepheader",action='store_false',help="if set, VCF header lines not included in output file")

    #density tool subparser
    dense_parser = subparsers.add_parser("density",help="Filter VCF file so only chunks of SNPs of sufficient size and closeness remain")
    dense_parser.add_argument("vcf",metavar="VCF",help="input vcf file")
    dense_parser.add_argument("-s","--size",type=int,default=1,help="minimum size for chunk to be accepted,default is 1")
    dense_parser.add_argument("-d","--distance",type=int,default=100,help="maximum distance between SNP i and i+1 for both to be in same group,default is 100")
    dense_parser.add_argument("-o","--output",default="out.vcf",help="output file name, default is 'out.vcf'")

    #filter tool subparser
    filter_parser = subparsers.add_parser("filter",help="Filter VCF file by value of metrix in INFO file")
    filter_parser.add_argument("vcf",metavar="VCF",help="input vcf file")
    filter_parser.add_argument("-m","--metric",default="CHI2",help="metric in INFO column to be used, default is CHI2")
    filter_parser.add_argument("-v","--value",type=float,default=0.05,help="cutoff value for metric, default=0.05")
    filter_parser.add_argument("-o","--output",default="out.vcf",help="output file name, default is 'out.vcf'")
    filter_parser.add_argument("-r","--greater",action='store_true',help="If set, tool will only select SNPs with metric value greater than or equal to cutoff value, default is lesser than or equal to")

    #graph tool subparser
    graph_parser = subparsers.add_parser("graph",help="Graphical representation of SNPs based value in INFO field and location on chromosome")
    graph_parser.add_argument("vcf",metavar="VCF",help= "input vcf file")
    graph_parser.add_argument("-m","--metric",default="CHI2",help="metric in INFO column to be used,default=CHI2")
    graph_parser.add_argument("-c","--chrom",default=None,help="which chromosome from the VCF file to graph, default is to plot each on separate graph")
    return parser

class VCFWorker(object):
    '''class that performs association testing on a VCF file using multiple processes to analyze multiple
        SNPs simultaneously'''
    def __init__(self, opts):
        '''initialize object and start worker processes'''
        self.numprocs = opts.numprocs
        if self.numprocs < 1:
            self.numprocs = 1
        self.pvalue = opts.pvalue
        if self.pvalue < 0.0:
            self.pvalue = 0.0
        elif self.pvalue > 1.0:
            self.pvalue = 1.0
        self.infile = open(opts.vcf,'r')
        self.inq = multiprocessing.Queue()
        self.outq = multiprocessing.Queue()
        self.usefisher = opts.fisher
        self.condense = opts.condense
        self.outfile = opts.output
        self.acceptlow = opts.acceptlow
        sample_handle= open(opts.samples)
        samples = [row for i,row in enumerate(sample_handle)]
        sample_handle.close()
        if not opts.groupsize:
            opts.groupsize = len(samples) / 2
        if opts.groupsize < 1:
            opts.groupsize = 1
        elif opts.groupsize >= len(samples):
            opts.groupsize = len(samples) - 1
        group1samples = [x[:-2] for x in samples[:int(opts.groupsize)]]
        group2samples = [x[:-2] for x in samples[int(opts.groupsize):]]
        self.group1samples = []
        self.vcfprefix = []
        newheadersadded = False
        for line in self.infile:
            #option to keep or leave out vcf header
            if opts.keepheader:
                #option to include new info headers for new data (not yet tested if it affects use in other tools)
                if "##INFO" in line and not newheadersadded and opts.addheader:
                    newheadersadded = True
                    fet_info_row = '##INFO=<ID=FET,Number=1,Type=Integer,Description="Fisher\'s Exact Test p-value based on GROUP1 and GROUP2\n'
                    CHI2_info_row = '##INFO=<ID=CHI2,Number=1,Type=Integer,Description="Pearson\'s Chi-Square Test of Independence p-value based on GROUP1 and GROUP2\n'
                    group1_info_row = '##GROUP1 = ' + " ".join(group1samples) + '\n'
                    group2_info_row = '##GROUP2 = ' + " ".join(group2samples) + '\n'
                    self.vcfprefix.append(fet_info_row)
                    self.vcfprefix.append(CHI2_info_row)
                    self.vcfprefix.append(group1_info_row)
                    self.vcfprefix.append(group2_info_row)
                self.vcfprefix.append(line)
            #this line contains the names of each sample in the file
            if "#CHROM" in line:
                self.vcfprefix.append(line)
                line = line.split("\t")
                for x in xrange(9,len(line)):
                    if x == len(line) - 1:
                        line[x] = line[x][:-1]
                    if line[x] in group1samples:
                        self.group1samples.append(x)
                break
        if opts.bonferroni:
            SNPcount = 0.0
            for line in self.infile:
                SNPcount += 1.0
            if SNPcount == 0.0:
                SNPcount = 1.0
            self.pvalue /= SNPcount
            self.infile.seek(0)
            for line in self.infile:
                if "#CHROM" in line:
                    break

        #start the worker processes
        self.pin = multiprocessing.Process(target=self.parse_input_vcf, args=())
        self.pout = multiprocessing.Process(target=self.write_output_vcf, args=())
        self.ps = [ multiprocessing.Process(target=self.process_row, args=()) for i in range(self.numprocs)]

        self.pin.start()
        self.pout.start()
        for p in self.ps:
            p.start()

        #do not let worker processes finish before input process
        self.pin.join()
        c = 1
        for p in self.ps:
            p.join()
            print ("Done %i" % (c))
            c += 1

        self.pout.join()

        self.infile.close()
    
    def parse_input_vcf(self):
            '''read VCF file as input and parse into chunks to enqueue for worker processes'''
            chunk = []
            counter = 0
            for i, row in enumerate(self.infile):
                if(self.inq.qsize() > 100):
                    while not self.inq.empty():
                        time.sleep(.1)
                chunk.append(row)
                if len(chunk) == 1000:
                    self.inq.put( (counter, chunk) )
                    chunk = []
                    counter += 1
            if chunk != []:
                self.inq.put((counter,chunk))
            for i in range(self.numprocs):
                self.inq.put("STOP")
    
    def process_row(self):
        '''create contingency table and perform assocation test on a chunk of the VCF file'''
        outputchunk = []
        for i,chunk in iter(self.inq.get, "STOP"):
            for row in chunk:
                line = row.split('\t')
                #determine refence allele and alternate allele(s)
                alternates = line[4].split(",")
                allele_array = [line[3]] + [val for val in alternates]
                alleles1 = {}
                alleles2 = {}
                #calculate biallelic genotype for each sample for this SNP and aggregate the counts
                #in bins depending on which group it belongs in
                for index in range(9,len(line)):
                    if line[index][0] == ".":
                        if line[index] == ".":
                            allele = "0"
                        else:
                            allele = ["0","0"]
                    else:
                        if "/" in line[index]:
                            allele = line[index][:line[index].find(":")].split("/")
                        elif "|" in line[index]:
                            allele = line[index][:line[index].find(":")].split("|")
                        else:
                            allele = line[index]
                    if self.condense == 3:
                        if len(allele) == 1:
                            continue
                        if allele == ["0","0"]:
                            allele = "ref-ref"
                        elif allele[0] == allele[1]:
                            allele = "alt-alt"
                        else:
                            allele = "ref-alt"
                    elif self.condense == 2:
                        if allele == ["0","0"] or allele == "0":
                            allele = "ref"
                        else:
                            allele = "alt"
                    else:
                        if len(allele) == 2:
                            allele = allele_array[int(allele[0])] + "-" + allele_array[int(allele[1])]
                        else:
                            allele = allele_array[int(allele)]

                    if index in self.group1samples:
                        if allele not in alleles1:
                            alleles1[allele] = 1
                        else:
                            alleles1[allele] += 1
                    else:
                        if allele not in alleles2:
                            alleles2[allele] = 1
                        else:
                            alleles2[allele] += 1
                #turn genotype dictionaries into 2 by M table (allows for 0 values)
                matrix_group1 = []
                matrix_group2 = []
                allelesinorder = []
                for allele in alleles1:
                    allelesinorder.append(allele)
                    matrix_group1.append(alleles1[allele])
                    if allele not in alleles2:
                        matrix_group2.append(0)
                    else:
                        matrix_group2.append(alleles2[allele])
                for allele in alleles2:
                    if allele not in alleles1:
                        allelesinorder.append(allele)
                        matrix_group2.append(alleles2[allele])
                        matrix_group1.append(0)
                #calculate fisher's exact test for table if set
                if self.usefisher:
                    fet_score = 1.0
                    if len(matrix_group1) != 1:
                        fet_score = twobymfishersexact(matrix_group1,matrix_group2,self.pvalue)
                    if fet_score <=self.pvalue:
                        line[7] += ";FET=" + str(fet_score)
                    else:
                        continue
                else:
                    #allow user to ignore entries where the contingency table has small cell values
                    toosmall = False
                    for val in (matrix_group1 + matrix_group2):
                        if val < 5 and not self.acceptlow:
                            toosmall = True
                            break
                    if toosmall:
                        continue
                    chi_score = 1.0
                    if len(matrix_group1) != 1:
                        a=1#chi_score = scipy.stats.chi2_contingency([matrix_group1,matrix_group2])[1]
                    if chi_score <= self.pvalue:
                        line[7] += ";CHI2=" + str(chi_score)
                    else:
                        continue
                
                #data is sent to output process for writing to output VCF file and table file
                str_matrix_group1 = [str(x) for x in matrix_group1]
                str_matrix_group2 = [str(x) for x in matrix_group2]
                tableentry = line[0] + '\t' + line[1] + '\tref:' + line[3] + '\n' + '\t'.join(allelesinorder) + '\n'+'\t'.join(str_matrix_group1) + '\n' + '\t'.join(str_matrix_group2) + '\n'
                # @@@ used as delimiter between VCF and table file data
                outputchunk.append('\t'.join(line) + "@@@" + tableentry)
            self.outq.put( (i, outputchunk ) )
            outputchunk = []
        self.outq.put("STOP")

    def write_output_vcf(self):
        '''process that writes data to output files'''
        cur = 0
        stop = 0
        buffer = {}
        outfile = open(self.outfile, "w")
        tablefile = open(self.outfile + '.table',"w")
        outfile.write(''.join(self.vcfprefix))
        #Keep running until we see numprocs STOP messages
        for works in range(self.numprocs):
            #use buffer in case chunks of VCF file are finished out of order so they are assembled in order 
            for i, chunk in iter(self.outq.get, "STOP"):
                if cur != i:
                    buffer[i] = chunk
                else:
                    cur += 1
                    for val in chunk:
                        val = val.split("@@@")
                        outfile.write(val[0])
                        tablefile.write(val[1])
                    while cur in buffer:
                        for val in buffer[cur]:
                            val = val.split("@@@")
                            outfile.write(val[0])
                            tablefile.write(val[1])
                        del buffer[cur]
                        cur += 1
        outfile.close()
        tablefile.close()


def filterbyinfo(infile,metric,cuttoff_value,outfile,greater):
    '''filter VCF file so that only SNPs with a value above or below a certain threshold
    for an item in the INFO column, where infile is the input VCF file, metric is the entry in the
    INFO column we are examining, cuttoff_value is a float that SNP values are compared too, outfile is the name
    of the output file, and greater is a boolean that decides if we accept SNPs with values greater or less than the cuttoff'''
    infile = open(infile,"r")
    outfile = open(outfile,"w")
    for line in infile:
        outfile.write(line)
        if "#CHROM" in line:
            break
    file_chunk = []
    for line in infile:
        line = line.split("\t")
        info = line[7].split(";")
        #find desired INFO metric, if not available, SNP is skipped
        for x in xrange(len(info)):
            if metric in info[x]:
                val = float(info[x][info[x].index('=')+1:])
                if (greater and val >= cuttoff_value) or (not greater and val <= cuttoff_value):
                    file_chunk.append('\t'.join(line))
                    if len(file_chunk) >= 10000:
                        outfile.write(''.join(file_chunk))
                        file_chunk = []
                    break
    outfile.write(''.join(file_chunk))
    del file_chunk
    outfile.close()
    infile.close()

def graphbyinfo(infile,metric,chrom):
    '''graph SNPs by location on chromosome and by chosen value in the INFO column, where infile 
    is the name of the input VCF file, metric is a string representing a name for a type of value to be graphed,
    and chrom can either be None, or a string representing which chromosome's SNPs to graph'''
    infile = open(infile,"r")
    for line in infile:
        if "#CHROM" in line:
            break
    current_group =''
    xaxis = []
    yaxis = []
    for line in infile:
        if "#CHROM" in line:
            continue
        line = line.split("\t")
        if chrom and chrom != line[0]:
            continue
        if current_group == '':
            current_group = line[0]
        elif current_group != line[0]:
            print current_group
            plt.scatter(xaxis,yaxis)
            plt.xlabel("location in chromosome")
            plt.ylabel("-log(value)")
            plt.show()
            xaxis = []
            yaxis = []
            current_group = line[0]
        info = line[7].split(";")
        for x in xrange(len(info)):
            if metric in info[x]:
                val = -1.0*math.log10(float(info[x][info[x].index('=')+1:]))
                yaxis.append(val)
                xaxis.append(int(line[1]))
    plt.scatter(xaxis,yaxis)
    plt.xlabel("location in chromosome")
    plt.ylabel("-log(value)")
    plt.show()
    infile.close()

def densityvcf(infile,width,density,outfile):
    '''select only those SNPs from a VCF file that are clustered together with 
    desired closeness and cluster size, infile is the name of the input file, width is an integer indicating
    how far two adjacent SNPs can be on a chromosome to be considered in the same group, density is an integer
    that represents the minimum size a group must be to be considered, and outfile is the output VCF file name'''
    infile = open(infile,"r")
    densityfile = open(outfile+".density",'w')
    outfile = open(outfile,"w")
    for line in infile:
        outfile.write(line)
        if "#CHROM" in line:
            break

    if density < 0:
        densityvcf = 0
    if width < 0:
        width = 0
    group = []
    groupvcf = []
    current_contig = ''
    for line in infile:
        line = line.split("\t")
        contig = line[0]
        loc = line[1]
        if current_contig == '':
            current_contig = contig
        elif current_contig != contig:
            if len(group) >= density:
                densityfile.write(contig + "\t" + "\t".join(group) + '\n')
                outfile.write(''.join(groupvcf))
            current_contig = contig
            group = []
            groupvcf = []
        #if new group, add SNP immediately
        if len(group) == 0:
            group.append(loc)
            groupvcf.append('\t'.join(line))
        #if next SNP is close enough to previous SNP, add to group
        elif int(loc) - int(group[len(group)-1]) <= width:
            group.append(loc)
            groupvcf.append('\t'.join(line))
        else:
            #if group is large enough, accept
            if len(group) >= density:
                densityfile.write(contig + "\t" + "\t".join(group)+ '\n' )
                outfile.write(''.join(groupvcf))
            group = []
            groupvcf = []
    if len(group) >= density:
        densityfile.write(contig + "\t" + "\t".join(group)+'\n')
        outfile.write(''.join(groupvcf))

    densityfile.close()
    outfile.close()
    infile.close()

def main(argv):
    '''Main function, takes in command line arguments and calls appropriate tool'''
    parser = make_parser()
    opts = parser.parse_args(argv)
    #parse the first argument to determine which function is being executed
    if opts.parser == "association":
        c = VCFWorker(opts)
    elif opts.parser == "density":
        densityvcf(opts.vcf,opts.distance,opts.size,opts.output)
    elif opts.parser == "filter":
        filterbyinfo(opts.vcf,opts.metric,opts.value,opts.output,opts.greater)
    elif opts.parser == "graph":
        graphbyinfo(opts.vcf,opts.metric,opts.chrom)
    return
   
if __name__ == '__main__':
    main(sys.argv[1:])
