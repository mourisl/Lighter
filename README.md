Lighter
=======

Described in: 
 
Song, L., Florea, L. and Langmead, B., [Lighter: Fast and Memory-efficient Sequencing Error Correction without Counting](http://genomebiology.com/2014/15/11/509/). Genome Biol. 2014 Nov 15;15(11):509.

Copyright (C) 2012-2013, and GNU GPL, by Li Song, Liliana Florea and Ben Langmead

Lighter includes source code from the [`bloom` C++ Bloom filter library](http://code.google.com/p/bloom/).  `bloom` is distributed under the [Mozilla Public License (MPL)](http://www.mozilla.org/MPL/) v1.1. 

### What is Lighter?

Lighter is a kmer-based error correction method for whole genome sequencing data.  Lighter uses sampling (rather than counting) to obtain a set of kmers that are likely from the genome.  Using this information, Lighter can correct the reads
containing sequence errors.

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/lighter), e.g. with `git clone https://github.com/mourisl/Lighter.git`
2. Run `make` in the repo directory

Lighter is small and portable, with [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and [zlib](http://en.wikipedia.org/wiki/Zlib) being the only library dependency. 


### Usage

    Usage: ./lighter [OPTIONS]
    OPTIONS:
    Required parameters:
	    -r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files
                         The file can be fasta and fastq, and can be gzip'ed with extension *.gz. 
                         When the input file is *.gz, the corresponding output file will also be gzip'ed.
	    -k kmer_length genome_size alpha
	    		or
	    -K kmer_length genome_size
    Other parameters:
	    -od output_file_directory: (default: ./)
	    -t num_of_threads: number of threads to use (default: 1)
	    -maxcor INT: the maximum number of corrections within a 20bp window (default: 4)
	    -trim: allow trimming (default: false)
	    -discard: discard unfixable reads. Will LOSE paired-end matching when discarding (default: false)
	    -noQual: ignore the quality socre (default: false)
	    -newQual ascii_quality_score: set the quality for the bases corrected to the specified score (default: not used)
            -saveTrustedKmers file: save the trusted kmers to specified file then stop (default: not used)
            -loadTrustedKmers file: directly get solid kmers from specified file (default: not used)
	    -zlib compress_level: set the compression level(1-9) of gzip (default: 1)
	    -h: print the help message and quit
	    -v: print the version information and quit

    NOTICE: genome_size does not need to be accurate, but it should be at least as large as the size of the sequenced genome.
            alpha is the sampling rate and decided by the user. A rule of thumb: alpha=(7/C), where C is the average coverage of the data set.
	    When using "-K" instead of "-k", Lighter will go through the reads an extra pass to decide C. And for "-K", genome_size should be relative accurate.
	    Lighter may adjust -maxcor and bad quality threshold if it decides the error rate of the data set is high. You can specify "-maxcor" explicitly to turn off this feature.
	    If you want to use longer kmer, you can adjust the constraints "MAX_KMER_LENGTH" in "utils.h" file before compiling lighter.

### Output

For each input file (specified by -r), Lighter will create a corresponding file containing the corrected reads in the directory specified by "-od". If the name of a read file is like A.fq or A.fastq, the corresponding output file name is A.cor.fq. For A.fa or A.fasta, the output file is A.cor.fa. For A.fastq.gz or A.fq.gz, the output file is A.cor.fq.gz. For A.fasta.gz or A.fa.gz, the output file is A.cor.fa.gz.

For each header line, Lighter will append some information regarding the correction result. 
	
	"cor": some bases of the sequence are corrected
	"bad_suffix", "bad_prefix": the errors in the suffix or prefix of indicated length could not be corrected
	"unfixable_error": the errors could not be corrected
	There are some two-characters string to indicate the reason why Lighter could not correct the errors:
		"lc": low coverage
		"ak": ambiguous candidate kmer
		"oc": over corrected

### Example

Suppose the data set is from E.Coli whose genome size is about 4.7M. If the data set has 3.3M reads in total with read length 100bp, then the average coverage is about 70x and we can set the alpha to be 7/70=0.1.

Single-end data set:

    ./lighter -r read.fq -k 17 5000000 0.1 -t 10

Paired-end data set:

    ./lighter -r left.fq -r right.fq -k 17 5000000 0.1 -t 10

Instead of "-k 17 5000000 0.1", you can use "-K 17 4700000" so that Lighter will go through the data set an extra pass to decide alpha. 

### Miscellaneous

Lighter is available as conda package and can be installed with:

```bash
conda install lighter --channel bioconda
```

A Docker container is available with:

```bash
docker run quay.io/mulled/lighter:1.1.1--1 lighter -h
```

If you have a [Galaxy](https://galaxyproject.org/) instance you can install Lighter from the [Galaxy Tool Shed](https://toolshed.g2.bx.psu.edu/view/bgruening/lighter).

If you are using Brew or Linuxbrew, you can install Lighter with:

```bash
brew install lighter
```

### Terms of use

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received (LICENSE.txt) a copy of the GNU General
Public License along with this program; if not, you can obtain one from
http://www.gnu.org/licenses/gpl.txt or by writing to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
### Support

Create a [GitHub issue](https://github.com/mourisl/lighter/issues).
