# AMAS
Calculate various summary statistics on a multiple sequence alignment

## Usage
AMAS can be run from the command line:

```shell
python AMAS.py <input_file> <format> <alphabet>
```
The supported formats are `fasta`, `phylip`, `nexus`, `phylip-int`, and `nexus-int`. The alphabets are `aa` or `dna`.

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
