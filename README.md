# AMAS
Alignment manipulation and summary statistics

## Installation

If you have [git installed](http://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your system, you can simply clone this repository and use the AMAS.py file as a stand-alone program.

You can also download AMAS from the [Python Package Index](https://pypi.python.org/pypi/amas/) or install it with [pip](https://pip.pypa.io/en/latest/installing.html):
```shell
pip install amas
```

## Usage
AMAS can be run from the command line. Here is the general usage:

```shell
usage: AMAS.py [-h] -i [IN_FILE [IN_FILE ...]] -f
               {fasta,phylip,nexus,phylip-int,nexus-int} -d {aa,dna} [-c] [-s]
               [-p CONCAT_PART] [-t CONCAT_OUT] [-o SUMMARY_OUT]

optional arguments:
  -h, --help            show this help message and exit
  -c, --concat          Concatenate input alignments
  -s, --summary         Print alignment summary
  -p CONCAT_PART, --concat-part CONCAT_PART
                        File name for the concatenated alignment partitions
  -t CONCAT_OUT, --concat-out CONCAT_OUT
                        File name for the concatenated alignment
  -o SUMMARY_OUT, --summary-out SUMMARY_OUT
                        File name for the concatenated alignment

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

You need to choose at least one action with `-c` (same as `--concat`) or `-s` (`--summary`) for the input to be processed. If you want to concatenate all DNA phylip files in a directory and all of them have the `.phy` extension, you can run:
```shell
python3 AMAS.py --in-format phylip --data-type dna --in-file *phy --concat
```
By default the output will be written to two file `partitions.txt`, containing partitions for your new alignment, and `concatenated-fasta.out` with the alignment itself in fasta format. You can change the default names for these files with `-p` (`--concat-part`) and `-t` (`--concat-out`), respectively, followed by the desired name.

Below is the same command using short versions of options:
```shell
python3 AMAS.py -f phylip -d dna -i *phy -c
```
You can summarize a protein fasta alignment by running:
```shell
python3 AMAS.py -f fasta -d aa -i my_alignment.nex -s
```
By default AMAS will write a file with the summary of the alignment in `summary.txt`. You can change the name of this file with `-o` or `--summary-out`.

You can perform concatenation and at the same time write summaries of the input alignments. Order in which options are supplied does not matter:
```shell
python3 AMAS.py -s -c -t all_gene_matrix.fas -f nexus-int -d dna -i *.nex
``` 

Also, AMAS can be imported to other Python modules:

```python
>>> from amas import AMAS

>>> in_file = 'fasta1.fas'
>>> in_format = 'fasta'
>>> data_type = 'dna'
>>> aln = AMAS.DNAAlignment(in_file, in_format, data_type)
>>> aln.summarize_alignment()
['fasta1.fas', '10', '100', '1000', '1', '0.1', '2', '0.02', '1', '0.01']
```
