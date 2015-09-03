# AMAS
Alignment manipulation and summary statistics

## Installation

You can download this repository zipped (button on the right-hand side of the screen) and use `AMAS.py` in the `amas` directory as a stand-alone program or clone it if you have [git installed](http://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your system.

If your system doesn't have a Python version 3.0 or newer, you will need to [download and install it](http://www.python.org/downloads/). On Linux-like systems (including Ubuntu) you can install it from the command line using

```
sudo apt-get install python3
```
To use `AMAS` as a Python package you can get it with  [pip](https://pip.pypa.io/en/latest/installing.html) from the [Python Package Index](https://pypi.python.org/pypi/amas/):
```
pip install amas
```
See below for the instructions on how to use this program as a Python module.

# Command line interface
`AMAS` can be run from the command line. Here is the general usage (you can view this in your command line with `python3 AMAS.py -h`):

```
usage: AMAS.py [-h] -i [IN_FILES [IN_FILES ...]] -f
               {fasta,phylip,nexus,phylip-int,nexus-int} -d {aa,dna} [-c] [-s]
               [-v] [-l SPLIT] [-r REPLICATE REPLICATE] [-p CONCAT_PART]
               [-t CONCAT_OUT] [-o SUMMARY_OUT]
               [-u {fasta,phylip,nexus,phylip-int,nexus-int}] [-e]

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
  -e, --check-align     Check if input sequences are aligned. Default: no
                        check

required named arguments:
  -i [IN_FILES [IN_FILES ...]], --in-files [IN_FILES [IN_FILES ...]]
                        Alignment files to be taken as input. You can specify
                        multiple files using wildcards (e.g. --in-files *fasta)
  -f {fasta,phylip,nexus,phylip-int,nexus-int}, --in-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        The format of input alignment
  -d {aa,dna}, --data-type {aa,dna}
                        Type of data
```

## Examples
For every `AMAS.py` run on the command line you need to provide:

1) input file name(s) with `-i` (or in long version: `--in-files`),

2) format with `-f` (`--in-format`),

3) and data type with `-d` (`--data-type`). 

The options available for the format are `fasta`, `phylip`, `nexus` (sequential), `phylip-int`, and `nexus-int` (interleaved). Data types are `aa` for protein alignments and `dna` for nucleotide alignments. 

You also need to choose at least one action with `-c` (same as `--concat`), `-s` (`--summary`), `-v` (`--convert`), or `-r` (`--replicate`) for the input to be processed. The order in which arguments are given does not matter.

IMPORTANT! `AMAS` is fast and powerful, but be careful: it assumes you know what you are doing and will not prevent you overwriting a file. It will, however, print out a warning if this has happened. You also need to be mindful of the input format specified, as incorrect format may result in unexpected program behavior. `AMAS` performs some basic checks to see if the parsing was successful but there is a trade-off between automated file format detection and computation times. In short, you can expect your output to be garbled if you don't get the input format right. `AMAS` was designed specifically to work with aligned sequences and it should not be used on unaligned files. You can turn on checking whether your files contain only aligned sequences with `-e` (`--check-align`) but this will increase computation times. 

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
Alignments to be concatenated need not have identical sets of taxa before processing: the concatenated alignment will be populated with missing data where a given locus is missing a taxon. However, if every file to be concatenated includes only unique names (for example species name plus sequence name: `D_melanogaster_NW_001845408.1` in one alignment, `D_melanogaster_NW_001848855.1` in other alignment etc.), you will first need to trim those names so that sequences from one taxon have equivalents in all files.   

Note that it takes `AMAS` longer to write an interleaved file than a sequential one, which may be an issue if you are concatenating to a large alignment on a laptop or an older desktop computer.

### Getting alignment statistics
This is an example of how you can summarize two protein fasta alignments by running:
```
python3 AMAS.py -f fasta -d aa -i my_aln.fasta my_aln2.fasta -s
```
By default `AMAS` will write a file with the summary of the alignment in `summary.txt`. You can change the name of this file with `-o` or `--summary-out`. You can also summarize a single or multiple sequence alignments at once. 

The statistics calculated include the number of taxa, alignment length, total number of matrix cells, overall number of undetermined characters, percent of missing data, AT and GC contents (for DNA alignments), number and proportion of variable sites, number and proportion of parsimony informative sites, and counts of all characters of the relevant amino acid or nucleotide alphabet.

### Converting among formats
To convert all nucleotide fasta files with a `.fas` extension in a directory to nexus alignments, you could use:
```
python3 AMAS.py -d dna -f fasta -i *fas -v -u nexus
```
In the above, the required options are combined with `-v` (`--convert`) action to convert and `-u nexus` indicating the output format.

`AMAS` will not overwrite over input here but will create new files instead, automatically appending appropriate extensions to the input file's name: `-out.fas`, `-out.phy`, `-out.int-phy`, `-out.nex`, or `-out.int-nex`.

### Splitting alignment by partitions
If you have a partition file, you can split a concatenated alignment and write a file for each partition:
```
python3 AMAS.py -f nexus -d dna -i concat.nex -l  partitions.txt -u nexus
```
In the above one input file `concat.nex` was provided for splitting with `-i` (can also use `--in-files`) and partitions file `partitions.txt` with `-l` (same as `--split`). For splitting you can only use one input and one partition file at a time. This is an example partition file:
```
  AApos1&2  =  1-604\3, 2-605\3
  AApos3  =  3-606\3
  28SAutapoInDels=7583, 7584, 7587, 7593
```
If this was the `partitions.txt` file from the example command above, `AMAS` would write three output files called `concat_AApos1&2.nex`, `concat_AApos3.nex`, and `concat_28SautapoInDels.nex`. The partitions file will be parsed correctly as long as there is no text prior to the partition name (`CHARSET AApos1&2` or `DNA, AApos1&2` will not work) and commas separate ranges or individual sites in each partition.

### Creating replicate data sets
With `AMAS` you can create concatenated alignments from a number of randomly chosen alignments that can be used for, for example, a phylogenetic jackknife analysis. Say you have 1000 phylip files, each containing a single aligned locus, and you want to create 200 replicate phylip alignments, each built from 100 loci randomly chosen from all the input files. You can do this by supplying the `-r` or `--replicate` followed by the number of replicates (in this case `200`) and number of alignments (`100`). Remember to supply the output format with `-u` if you want it to be other than fasta:
```
python3 AMAS.py -d dna -f phylip -i *phy -r 200 100 -u phylip
```
### Combining options
You can get statistics for all input alignments, convert them to phylip, and concatenate (also to a phylip file) in one go by simply combining actions:
```
python3 AMAS.py -d aa -f fasta -i *fas -c -s -v -u phylip
```
### Checking if input is aligned
By specifying optional argument `-e` (`--check-align`), you can make `AMAS` check if your input files contain only aligned sequences. This option is disabled by default because it can substantially increase computation times in files with many taxa. Enabling this option also provides an additional check against misspecified input file format.

# AMAS as a Python module
Using `AMAS` inside your Python pipeline gives you much more flexibility in how the input and output are being processed. All the major functions of the command line interface can recreated using `AMAS` as a module. Following installation from [pip](https://pip.pypa.io/en/latest/installing.html) you can import it with:

```python
from amas import AMAS
```
The class used to manipulate alignments in `AMAS` is `MetaAlignment`. This class has to be instantiated with the same, named arguments as on the command line: `in_files`, `data_type`, `in_format`. MetaAlignment holds one or multiple alignments and its `in_files` option must be a list, even if only one file is being read.
```python

meta_aln = AMAS.MetaAlignment(

in_files=["gene1.phy"], data_type="dna",in_format="phylip"

)
```
Creating MetaAlignment with multiple files is easy:
```python
multi_meta_aln = AMAS.MetaAlignment(

in_files=["gene1.phy", "gene1.phy"], data_type="dna", in_format="phylip"

)
```
Now you can call the various methods on your alignments. `.get_summaries()` method will compute summaries for your alignments and produce headers for them as a tuple with first element being the header and the second element a list of lists with the statistics:
```python
summaries = meta_aln.get_summaries()
```
The header is different for nucleotide and amino acid data. You may choose to skip it and print only the second element of the tuple, that is a list of summary statistics:
```python
statistics = summaries[1]
```
`.get_parsed_alignments()` returns a list of dictionaries where each dictionary is an alignment and where taxa are the keys and sequences are the values. This allows you to, for example, print only taxa names in each alignment:
```python3 
# get parsed dictionaties
aln_dicts = multi_meta_aln.get_parsed_alignments()

# print only taxa names in the alignments:
for alignment in aln_dicts:
    for taxon_name in alignment.keys():
        print(taxon_name)
```
To split alignment use `.get_partitioned("your_partitions_file")` on a `MetaAlignment` with a single input file. `.get_partitioned()` returns a list of dictionaries of dictionaries, with `{ partition_name : { taxon : sequence } }` structure for each partition:
```python
partitions = meta_aln.get_partitioned("partitions.txt")
```
`AMAS` uses `.get_partitions("your_partitions_file")` to parse the partition file:
```python
parsed_parts = meta_aln.get_partitions("partitions.txt")
print(parsed_parts)
```

`.get_replicate(no_replicates, no_loci)` gives a list of parsed alignments (dictionaries), each a replicate constructed from the specified number of loci:
```python
replicate_sets = multi_meta_aln.get_replicate(2, 2)
```
To concatenate multiple alignments first parse them with `.get_parsed_alignments()`, then pass to `.get_concatenated(your_parsed_alignments)`. This will return a tuple where the first element is the `{ taxon : sequence }` dict
of concatenated alignment and the second element is the partitions dict with `{ name : range }`.
```python
parsed_alns = multi_meta_aln.get_parsed_alignments()
concat_tuple = multi_meta_aln.get_concatenated(parsed_alns)
concatenated_alignments = concat_tuple[0]
concatenated_partitions = concat_tuple[1]
```
To print to file or convert among file formats use one of the `.print_format(parsed_alignment)` methods called with a parsed dictionary as an argument. These methods include `.print_fasta()`, `.print_nexus()`, `.print_nexus_int`, `print_phylip()`, and `.print_phylip_int()`. They return an apporpriately formatted string.
```python
for alignment in aln_dicts:
    nex_int_string = meta_aln.print_nexus_int(alignment)
    print(nex_int_string)
```
