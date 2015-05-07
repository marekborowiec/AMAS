# AMAS
Calculate various summary statistics on a multiple sequence alignment

## Usage
AMAS can be run from the command line:

```shell
usage: AMAS.py [-h] --in-file [IN_FILE [IN_FILE ...]] --in-format
               {fasta,phylip,nexus,phylip-int,nexus-int} --data-type {aa,dna}

Calculate various statistics on a multiple sequence alignment

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  --in-file [IN_FILE [IN_FILE ...]]
                        Alignment files to be taken as input. You can specify
                        multiple input files using wildcards (e.g. --in-file
                        *fasta)
  --in-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        The input alignment format
  --data-type {aa,dna}  Type of data
```

Also, AMAS can be imported from other Python modules:

```python
>>> from amas import AMAS

>>> in_file = 'fasta1.fas'
>>> in_format = 'fasta'
>>> data_type = 'dna'
>>> aln = AMAS.DNAAlignment(in_file, in_format, data_type)
>>> aln.summarize_alignment()
['fasta1.fas', '10', '100', '1000', '1', '0.1', '2', '0.02', '1', '0.01']
```
