Lighter
=======

Described in: 
 
Song, L., Florea, L. and Langmead, B (2014). Lighter: Fast and Memory-efficient Error Correction without Counting

Copyright (C) 2012-2013, and GNU GPL, by Li Song, Liliana Florea and Ben Langmead

Lighter includes source code from the [`bloom` C++ Bloom filter library](http://code.google.com/p/bloom/).  `bloom` is distributed under the [Mozilla Public License (MPL)](http://www.mozilla.org/MPL/) v1.1. 

### What is Lighter?

Lighter is a kmer-based error correction method for whole genome sequencing data.  Lighter uses sampling (rather than counting) to obtain a set of kmers that are likely from the genome.  Using this information, Lighter can correct the reads
containing sequence errors.

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/lighter), e.g. with `git clone https://github.com/mourisl/Lighter.git`
2. Run `make` in the repo directory

Lighter is small and portable, with [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) being the only library dependency.  We have successfully built it on Linux, Mac OS X, and Windows.  To build on Windows, you will need to download the `pthreadsGC2.dll` library from the [Pthreads Win32 library](http://www.sourceware.org/pthreads-win32/) and copy it to the repo directory with the Lighter sources.

### Usage

    Usage: ./lighter [OPTIONS]
    OPTIONS:
    Required parameters:
	    -r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files
	    -k kmer_length genome_size alpha
    Other parameters:
	    -od: output_file_directory (default: ./)
	    -t: number of threads to use (default: 1)
	    -trim: allow trimming (default: false)
	    -discard: discard unfixable reads (default: false)
	    -maxcor: the maximum number of correction for within a kmer_length window (default: 4)

    NOTICE: genome_size does not need to be accurate, but it should be at least as large as the size of the sequenced genome.
            alpha is decided by the user. A rule of thumb: alpha=(3.5/C), where C is the coverage of the data set.

### Example

Suppose the data sets' coverage is about 70x:

Single-end data set:

    ./lighter -r read.fq -k 17 5000000 0.05 -t 10

Paired-end data set:

    ./lighter -r left.fq -r right.fq -k 17 5000000 0.05 -t 10

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
