SNPregator  V1.0.0

Author: Alexander Paul (agentscamp5)

SNPregator is a set of command-line tools used to collect aggregate data on Single Nucleotide Polymorphisms (SNPs) , as well as small Insertions and Deletions in the genome (INDELs). This data can be used to perform assocation tests like Fisher's Exact Test and Pearson's Chi-Square Test to determine if a SNP is associated with a particular group of samples being studied. The SNPs can also be filtered by these results, or by how density they are packed together, and can also be visualized graphically.

Requirements:
Python 2.7+
Numpy, Scipy, and Matplotlib libraries
Unix-based Operating System (for multiprocessing)

Example VCF file:

example.vcf:

...
...
...
"#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  A   B   C   D   E   F   G   H   I   J   K"
chrGroup1.1 1057    .   TATTTGGGA   T   999 PASS    AC=68;AF=0.773;AN=88;DP=0;FQ=999;INDEL;set=Intersection GT:GQ:PL    1/1:12:129,9,0  1/1:61:255,60,0 0/1:71:65,0,249 1/1:99:255,81,0 1/1:99:255,57,0 1/1:11:143,9,0  1/1:6:60,3,0    1/1:11:143,9,0  1/1:6:60,3,0    0/1:99:206,0,243    1/1:49:255,48,0
chrGroup1.1 1116    .   T   A   161 PASS    AC=17;AF=0.293;AN=58;DP=0;G3=6.694e-07,1,3.567e-91;HWE=0.00271;set=Intersection GT:GQ:PL    0/0:69:0,66,201 0/0:51:0,49,255 0/0:36:0,36,242 ./. 0/1:16:13,0,255 0/1:12:9,0,154  0/0:17:0,12,108 ./. ./. ./. ./.
chrGroup1.1 1152    .   G   A   999 PASS    AC=30;AF=0.385;AN=78;DP=0;set=Intersection  GT:GQ:PL    0/0:7:0,5,125   0/1:42:42,0,88  0/1:9:5,0,136   1/1:5:43,6,0    0/0:9:0,9,109   0/0:7:0,5,125   0/1:47:46,0,168 0/1:32:53,0,29  0/0:9:0,9,109   0/1:16:19,0,131 0/1:92:89,0,143
chrGroup1.1 1198    .   TATTTAGATAA T   195 PASS    AC=3;AC1=3;AF=0.30;AF1=0.3193;AN=10;DP=0;DP4=13,15,1,7;FQ=196;INDEL;MQ=56;PV4=0.11,4.4e-06,0.048,1;VDB=0.0186;set=Scout GT:GQ:PL    ./. ./. ./. ./. ./. ./. ./. ./. ./. ./. ./.
chrGroup1.1 1199    .   ATTTAGATAAT A   999 PASS    AC=37;AF=0.420;AN=88;DP=0;INDEL;set=Intersection    GT:GQ:PL    0/1:99:129,0,105    0/1:3:29,3,0    0/0:7:0,6,97    0/1:82:80,0,255 0/1:99:113,0,255    0/1:54:52,0,255 0/1:75:98,0,72  0/1:99:140,0,255    0/0:3:0,3,60    0/1:42:39,0,255 0/1:54:52,0,255
chrGroup1.1 1348    .   C   T   999 PASS    AC=47;AF=0.534;AN=88;DP=0;set=Intersection  GT:GQ:PL    0/1:73:70,0,130 0/1:87:84,0,112 0/1:56:57,0,56  0/1:56:53,0,155 0/1:34:31,0,83  0/1:99:108,0,249    0/1:50:49,0,84  0/1:47:45,0,50  0/1:56:53,0,155 0/0:4:0,3,40    0/1:87:84,0,112
chrGroup1.1 1523    .   A   T   999 PASS    AC=69;AF=0.784;AN=88;DP=0;set=Intersection  GT:GQ:PL    1/1:15:64,12,0  1/1:11:57,9,0   1/1:48:161,48,0 0/1:40:31,0,138 0/1:78:71,0,102 1/1:36:71,18,0  1/1:11:56,9,0   0/1:29:23,0,156 1/1:72:163,54,0 0/1:32:48,0,37  1/1:30:90,27,0
chrGroup1.1 2300    .   C   CATATAT,CATAT   159 PASS    AC=16,5;AC1=10;AF=0.533,0.167;AF1=1;AN=30;DP=0;DP4=0,0,1,0;FQ=-33;INDEL;MQ=60;PV4=1,0.46,1,1;VDB=0.0188;set=Intersection    GT:GQ:PL    ./. ./. ./. ./. ./. 1/1:4   ./. ./. ./. 0/2:3   ./.
...
...
...

Example sample name file:

One sample name per line, list group 1 samples first, then group 2 samples, size of group 1 indicated when calling SNPregator.py association

Group 1 = A,B,C,F,K
Group 2 = D,E,G,H,I,J

example_samples.txt:
A
B
C
F
K
D
E
G
H
I
J

SNPregator:

Main function, if called, will inform you of available tools

python SNPregator.py -h
usage: SNPregator.py [-h] {association,density,filter,graph} ...

Tools to test association between groups in VCF files and provides filters for S
NPs

optional arguments:
  -h, --help            show this help message and exit

tools:
  {association,density,filter,graph}
                        tool list
    association         Create contingency tables of aggregate allele counts
                        for a VCF file and perform tests of assocation
    density             Filter VCF file so only chunks of SNPs of sufficient
                        size and closeness remain
    filter              Filter VCF file by value of metrix in INFO file
    graph               Graphical representation of SNPs based value in INFO
                        field and location on chromosome

association

perform a test of assocation on SNPs and INDELs in a VCF file based on provided grouping

python SNPregator.py association -h
usage: SNPregator.py association [-h] [-n NUMPROCS] [-a] [-o OUTPUT]
                                         [-b] [-q {1,2,3}] [-f] [-g GROUPSIZE]
                                         [-p PVALUE] [-d] [-k]
                                         VCF SAMPLES

positional arguments:
  VCF                   input vcf file
  SAMPLES               file with sample names

optional arguments:
  -h, --help            show this help message and exit
  -n NUMPROCS, --numprocs NUMPROCS
                        Number of worker processes to launch, default is
                        number of CPUs in computer
  -a, --acceptlow       if set, chi-square test will not be performed if table
                        has cell with value < 5
  -o OUTPUT, --output OUTPUT
                        output file name, default is 'out.vcf'
  -b, --bonferroni      if set, apply the bonferroni correction for multiple
                        comparisons
  -q {1,2,3}, --condense {1,2,3}
                        Condense SNPs to [1: raw input(unchanged); 2: ref-
                        ref,alt-* ;3: ref-ref,ref-alt,alt-alt], default is 1
  -f, --fisher          if set, use Fisher Exact Test metric not Pearson Chi-
                        Square Test
  -g GROUPSIZE, --groupsize GROUPSIZE
                        sample size for group 1, default is half of total
                        sample size
  -p PVALUE, --pvalue PVALUE
                        cutoff pvalue for significance, default is 0.05
  -d, --addheader       if set, new header lines describing association test
                        metrics and group names are added
  -k, --keepheader      if set, VCF header lines not included in output file

Example Usage:

    python SNPregator association example.vcf example_samples.txt -g 5 -bkf

    This will calculate the Fisher's Exact Test p-value for each SNP in example.vcf, grouping the first 5 sample names
    in example_samples.txt as group 1 and the rest as group 2, it will apply the bonferonni correction, and any headers in the 
    original VCF file will not be included in the output file, which will be out.vcf, and also produce an output file called
    out.vcf.table, that has the contigency tables for each SNP stored in tab-delimited format

    example *.table file:
    ...
    ...
    ...
    chrGroup1.1 2413    ref:T
    T-T T-TAN   TAN-TAN
    17  5   0
    20  0   2
    chrGroup1.1 5357    ref:C
    C-G C-C
    16  6
    7   15
    ...
    ...
    ...

    python SNPregator association example.vcf example_samles.txt -g 5 -o otherout.vcf -p 0.001 -dn 2

    This will calculate the Chi-square test p-value for each  SNP in example.vcf, grouping the first 5 sample names
    in example_samples.txt as group 1 and the rest as group 2, with a significance cutoff level of 0.001, it will add new header lines
    into the output VCF file describing the group composition and the new test meanings in the iNFO column, the file will be outputted to
    otherout.vcf, and only 2 processes will be used congruently to process the data in the input VCF file.Also produced is an output file called
    otherout.vcf.table, that has the contigency tables for each SNP stored in tab-delimited format

density

Filter a VCF file so that only SNPs that are grouped together with desired closeness and group size are left represented.

python SNPregator.py density -h
usage: SNPregator.py density [-h] [-s SIZE] [-d DISTANCE] [-o OUTPUT]
                                     VCF

positional arguments:
  VCF                   input vcf file

optional arguments:
  -h, --help            show this help message and exit
  -s SIZE, --size SIZE  minimum size for chunk to be accepted,default is 1
  -d DISTANCE, --distance DISTANCE
                        maximum distance between SNP i and i+1 for both to be
                        in same group,default is 100
  -o OUTPUT, --output OUTPUT
                        output file name, default is 'out.vcf'

Example Usage:

    python SNPregator.py density example.vcf -s 5 -d 40  -o denseoutput.vcf

    This will accept example.vcf and search for groups of SNPs with each SNP no more than 40 base pairs away from the
    next SNP on the chromosome and a minimum size (# of snps in group) of 5. Two output files will be produced, denseoutput.vcf,
    which will be a regular VCF file that only contains SNPs found in these groups, and another file named denseoutput.vcf.dense,
    which will contain the SNP locations for each group on a per chromosome and per group basis

    exaple *.dense file:
    chrGroup1.1
    1234,1300,1345,1350,1400
    4566,5500,6544,7000,7200,8000,8303
    chrGroup1.2
    3,430,500,654,1000
    2000,2034,2055,2200,2425

filter

Filter a VCF file so that only SNPs with a certain type of value in the VCF INFO field within a certain range of values
are selected.

python SNPregator.py filter -h
usage: SNPregator.py filter [-h] [-m METRIC] [-v VALUE] [-o OUTPUT]
                                    [-r]
                                    VCF

positional arguments:
  VCF                   input vcf file

optional arguments:
  -h, --help            show this help message and exit
  -m METRIC, --metric METRIC
                        metric in INFO column to be used, default is CHI2
  -v VALUE, --value VALUE
                        cutoff value for metric, default=0.05
  -o OUTPUT, --output OUTPUT
                        output file name, default is 'out.vcf'
  -r, --greater         If set, tool will only select SNPs with metric value
                        greater than or equal to cutoff value, default is
                        lesser than or equal to

Example Usage:
    
    python SNPregator.py filter example.vcf -v 0.001 -m FET

    This will accept example.vcf as an input file and produce an output file out.vcf that contains only
    SNPs from example.vcf who had a Fisher's Exact Test p-value in their INFO column, and whose Fisher's p-value
    was less than or equal to 0.001

SNPregator graph:

Graphically display the assocation test p-values for different SNPs in a VCF file based on their location on a particular
chromosome.

python SNPregator.py graph -h
usage: SNPregator.py graph [-h] [-m METRIC] [-c CHROM] VCF

positional arguments:
  VCF                   input vcf file

optional arguments:
  -h, --help            show this help message and exit
  -m METRIC, --metric METRIC
                        metric in INFO column to be used,default=CHI2
  -c CHROM, --chrom CHROM
                        which chromosome from the VCF file to graph, default
                        is to plot each on separate graph

Example Usage:

    SNPregator.py graph example.vcf -m FET -c chrGroup1

    This will produce a scatterplot of SNPs in example.vcf that lie in Chromosome "chrGroup1" with location
    on chromosome on the x-axis and it's Fisher's Exact Test p-value on the y-axis.
