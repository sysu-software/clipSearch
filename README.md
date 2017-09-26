# clipSearch
A tool for identifying miRNA-target interactions from CLIP-seq peaks

Overview:
---------
clipSearch is a tool for for identifying miRNA-target interactions from CLIP-seq peaks<BR>

Usage:
---------
Usage:  clipSearch [options] <genome file> <genome fai> <mir file, fasta> <peak, bed><BR>
peak format is bed6+<BR>
[options]<BR>
-v/--verbose                   : verbose information<BR>
-V/--version                   : clipSearch version<BR>
-h/--help                      : help informations<BR>
-o/--output <string>           : output file<BR>
-m/--max-mfe <double>          : maximum mfe for miR-target duplex [default < 0]<BR>
-s/--min-score <int>           : minimum score for miR-target duplex [default > 0]<BR>

Installation:
---------
Download clipSearch.tar.gz from http://starbase.sysu.edu.cn/clipSearch/; unpack it, and make:<BR>
tar -xvf clipSearch.tar.gz <BR>
cd clipSearch<BR>
make<BR>
The newly compiled binary (clipSearch) is in the clipSearch /bin directory.<BR>

At this point you should have everything to run a built-in test data set.<BR>
cd test_data<BR>
./run_test.sh<BR>

System requirements:
---------
Operating system: clipSearch is designed to run on POSIX-compatible platforms, including UNIX, Linux and Mac OS/X. We have tested  most extensively on Linux and MacOS/X because these are the machines we develop on.<BR>
Compiler: The source code is compiled with  the C++ compiler g++. We test the code using the g++ compilers. Libraries and other installation requirements: clipSearch includes one software library: the RNAfold library package. All will automatically compile during clipSearch installation process. By default, clipSearch does not require any additional libraries to be installed by you.<BR>

Prerequisites:
---------
If everything goes well you can get started with your real data! :)<BR>
You need to have the reference genome<BR>
As an example, let's assume you use human genome (version hg19).<BR>
(1)	Genome: wget -c 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'<BR>
rm -f chrUn* *hap*.fa.gz *random.fa.gz<BR>
mkdir genome<BR>
gzip ¨Cd *.gz<BR>
cat chr*.fa > genome/hg19.fa; cd..<BR>
(2) build the fai index:<BR>
samtools faidx genome/hg19.fa<BR>
(3) miRNA file<BR>
You can get the miR file from mirbase.org<BR>
(4)AGO CLIP-seq peaks<BR>
You can download the CLIP-seq peaks from starBase download center<BR>
http://starbase.sysu.edu.cn/download.php<BR>

run clipSearch:
---------
bin/clipSearch ./test_data/testGenome.fa ./test_data/testGenome.fa.fai ./test_data/testMir.fa ./test_data/testPeak.bed >./test_data/test_clipSearch_mtis.txt<BR>

Output:
---------
#chrom, chromStart, chromEnd, miR:peak, score, strand, seedType, mfe, alignScore, miR(3'-5'), pairs, target(5'-3'), <BR>
\>PDCD4	112	134	hsa-miR-21-5p:AGOPeak1	11.00	+	8mer	-9.00	43.00<BR>
miRNA  3'--AGTTGTAGTCAGACTATTCGAT-5'<BR>
          -.:::.||-|.|||.||||||||<BR>
target 5'-AGTGGAAT-ATTCTAATAAGCTA-3'<BR>

Acknowledgements:
---------
Thanks a lot to everyone who contributed to the public code used by clipSearch.<BR>

Contact :
---------
/*******************************************************************************
 *	clipSearch - A tool for identifying miRNA-target interactions from CLIP-seq peaks
 *
 *	Author : Jian-Hua Yang yangjh7@mail.sysu.edu.cn
 * 
 *	School of Life Sciences, Sun Yat-Sen University
 *
 *  Create Time: 18/09/2010
 *******************************************************************************/
