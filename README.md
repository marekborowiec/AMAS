# AMAS
Alignment manipulation and summary statistics

## Installation

You can download this repository zipped (button on the right-hand side of the screen) and use `AMAS.py` in the `amas` directory as a stand-alone program or clone it if you have [git installed](http://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your system.

To use `AMAS` as a Python package in your programs you can get it with [pip](https://pip.pypa.io/en/latest/installing.html) from the [Python Package Index](https://pypi.python.org/pypi/amas/):
```shell
pip install amas
```

## Usage
`AMAS` can be run from the command line. Here is the general usage (you can view this in your command line with `python3 AMAS.py -h`):

```
usage: AMAS.py [-h] -i [IN_FILE [IN_FILE ...]] -f
               {fasta,phylip,nexus,phylip-int,nexus-int} -d {aa,dna} [-c] [-s]
               [-v] [-l SPLIT] [-r REPLICATE REPLICATE] [-p CONCAT_PART]
               [-t CONCAT_OUT] [-o SUMMARY_OUT]
               [-u {fasta,phylip,nexus,phylip-int,nexus-int}]

Alignment manipulation and summary statistics

optional arguments:
  -h, --help            show this help message and exit
  -c, --concat          Concatenate input alignments
  -s, --summary         Print alignment summary
  -v, --convert         Convert to other file format
  -l SPLIT, --split SPLIT
                        File name for partitions to be used for alignment
                        splitting.
  -r REPLICATE REPLICATE, --replicate REPLICATE REPLICATE
                        Create replicate data sets for phylogenetic jackknife
                        [replicates, no alignments for each replicate]
  -p CONCAT_PART, --concat-part CONCAT_PART
                        File name for the concatenated alignment partitions.
                        Default: 'partitions.txt'
  -t CONCAT_OUT, --concat-out CONCAT_OUT
                        File name for the concatenated alignment. Default:
                        'concatenated.out'
  -o SUMMARY_OUT, --summary-out SUMMARY_OUT
                        File name for the alignment summary. Default:
                        'summary.txt'
  -u {fasta,phylip,nexus,phylip-int,nexus-int}, --out-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        File format for the output alignment. Default: fasta

required named arguments:
  -i [IN_FILE [IN_FILE ...]], --in-file [IN_FILE [IN_FILE ...]]
                        Alignment files to be taken as input. You can specify
                        multiple files using wildcards (e.g. --in-file *fasta)
  -f {fasta,phylip,nexus,phylip-int,nexus-int}, --in-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        The format of input alignment
  -d {aa,dna}, --data-type {aa,dna}
                        Type of data
```

## Examples
For every `AMAS.py` run on the command line you need to provide:

1) input file name(s) with `-i` (or in long version: `--in-file`),

2) format with `-f` (`--in-format`),

3) and data type with `-d` (`--data-type`). 

The options available for the format are `fasta`, `phylip`, `nexus` (sequential), `phylip-int`, and `nexus-int` (interleaved). Data types are `aa` for protein alignments and `dna` for nucleotide alignments. 

You also need to choose at least one action with `-c` (same as `--concat`), `-s` (`--summary`), `-v` (`--convert`), or `-r` (`--replicate`) for the input to be processed. The order in which arguments are given does not matter.

### Concatenating alignments
For example, if you want to concatenate all DNA phylip files in a directory and all of them have the `.phy` extension, you can run:
```
python3 AMAS.py -f phylip -d dna -i *phy -c
```
By default the output will be written to two files: `partitions.txt`, containing partitions from which your new alignment was constructed, and `concatenated.out` with the alignment itself in the fasta format. You can change the default names for these files with `-p` (`--concat-part`) and `-t` (`--concat-out`), respectively, followed by the desired name. The output format is specified by `-u` (`--out-format`) and can also be any of the following: `fasta`, `phylip`, `nexus` (sequential), `phylip-int`, or `nexus-int` (interleaved).

Below is a command specifying the concatenated file output format as nexus with `-u nexus`:
```
python3 AMAS.py -f fasta -d aa -i *phy -c -u nexus
```
Alignments to be concatenated need not have identical sets of taxa before processing: the concatenated alignment will be populated with missing data where a given locus is missing a taxon.

Note that it takes `AMAS` about 10x longer to write an interleaved file than a sequential one, which may be a factor if you are concatenating to a large (>50MB) alignment on a laptop or an older desktop computer.

### Getting alignment statistics
This is an example of how you can summarize two protein fasta alignments by running:
```
python3 AMAS.py -f fasta -d aa -i my_aln.fasta my_aln2.fasta -s
```
By default `AMAS` will write a file with the summary of the alignment in `summary.txt`. You can change the name of this file with `-o` or `--summary-out`. You can also summarize a single or multiple sequence alignments at once. 

The statistics calculated include the number of taxa, alignment length, total number of matrix cells, overall number of undetermined characters, percent of missing data, AT and GC contents (for DNA alignments), number and proportion of variable sites, number and proportion of parsimony informative sites, and proportions of all characters relative to the matrix size.

### Converting among formats
To convert all nucleotide fasta files with a `.fas` extension in a directory to nexus alignments, you could use:
```
python3 AMAS.py -d dna -f fasta -i *fas -v -u nexus
```
In the above, the required options are combined with `-v` (`--convert`) action to convert and `-u nexus` indicating the output format.

### Splitting alignment by partitions
If you have a raxml-style partition file, you can split a concatenated alignment and write a file for each partition:
```
python3 AMAS.py -f nexus -d dna -i concat.nex -l  partitions.txt -u nexus
```
In the above one input file `concat.nex` was provided for splitting with `-i` (can also use `--in-file`) and partitions file `partitions.txt` with `-l` (same as `--split`). For splitting you can only use one input and one partition file at a time. This is an example partition file:
```
  AApos1&2  =  1-604\3, 2-605\3
  AApos3  =  3-606\3
  28SAutapoInDels=7583, 7584, 7587, 7593
```
If this was the `partitions.txt` file from the example above, `AMAS` would write three output files called `concat_AApos1&2.nex`, `concat_AApos3.nex`, and `concat_28SautapoInDels.nex`. The partitions file will be parsed correctly as long as there is no text prior to the partition name (`CHARSET AApos1&2` or `DNA, AApos1&2` will not work) and commas separate ranges or individual sites in each partition.

### Creating replicate data sets
With `AMAS` you can create concatenated alignments from a proportion of randomly chosen alignments that can be used for, for example, a phylogenetic jackknife analysis. Say you have 1000 phylip files, each containing a single aligned locus, and you want to create 200 replicate phylip alignments, each built from 100 loci randomly chosen from all the input files. You can do this by supplying the `-r` or `--replicate` followed by the number of replicates (in this case `200`) and number of alignments (`100`). Remember to supply the output format with `-u` if you want it to be other than fasta:
```
python3 AMAS.py -d dna -f phylip -i *phy -r 200 100 -u phylip
```
### Combining options
You can get statistics for all input alignments, convert them to phylip, and concatenate (also to a phylip file) in one go by simply combining actions:
```
python3 AMAS.py -d aa -f fasta -i *fas -c -s -v -u phylip
```
### AMAS as a Python package
AMAS can be also imported to other Python modules:

```python
>>> from amas import AMAS

>>> in_file = 'fasta1.fas'
>>> in_format = 'fasta'
>>> data_type = 'dna'
>>> aln = AMAS.DNAAlignment(in_file, in_format, data_type)
>>> aln.summarize_alignment()
['fasta1.fas', '10', '100', '1000', '1', '0.1', '2', '0.02', '1', '0.01']
```
