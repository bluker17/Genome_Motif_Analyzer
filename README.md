# Genom Motif Analyzer
> Please use the latest stable release. Download here: [Latest Release URL](hhttps://github.com/bluker17/Genome_Motif_Analyzer/releases/latest)

## Author
**Bobby Luker**    
rluker@charlotte.edu    
UNCC ID: 801484356

## Program Description

For a draft genome assembly, when given multiple BLAST result TSV files and a contig TXT file containing contig IDs and lengths, the contigs are filtered by bin priority, bitscore, coverage threshold, and contig size threshold. Summary statistics are generated for each bin, along with two bar plots showing the number of contigs per bin and the total base pairs per bin. The default coverage threshold is set to 0.9 to ensure that only the highest-quality contigs are used. However, the user can adjust the coverage threshold to relax this restriction if needed.


[Github Project URL](https://github.com/bluker17/Genome_Motif_Analyzer)

## License
**GNU General Public License v3.0**

Review `LICENSE` for details.

## Project File Structure:
```
└── 📁Genome_Motif_Analyzer
    └── 📁data
    └── 📁output
    └── 📁src
        └── 📁csv_output
            ├── __init__.py
            ├── output.py
        └── 📁file_reader
            ├── __init__.py
            ├── reader.py
        └── 📁motif_locator
            ├── __init__.py
            ├── locator.py
        └── 📁statistic_analysis
            ├── __init__.py
            ├── statistics.py
    └── 📁testing_materials
        └── 📁example_data
        └── 📁example_outputs
        └── 📁expected_example_outputs
        ├── run_test_1.sh
        ├── run_test_2.sh
        ├── run_test_3.sh
        ├── run_test_4.sh
    ├── dependencies.txt
    ├── environment.yml
    ├── LICENSE
    ├── main.py
    └── README.md
```

## Overview
`example_data`: Contains the example genome FASTA files and example format of motif information CSV files.

`src`: Contains multiple subdirectories leading to modules handling the input and output data. 
1. `file_reader` contains module `reader.py` which collects all the provided FASTA files and motifs to search from the provided motif CSV file.
2. `motif_locator` contains module `locator.py` which searches for all motif instances in the FASTA files using the Aho-Corasick algorithm. Results are passed onto the `statistics.py` for further analysis in the form of a dictionary.
3. `statistic_analysis` contain module `statistics.py` which performs a one-sided proportion tests on motif occurrences per strand and per FASTA file as provided from the `locator.py` findings. Statistical and locator results are saved into a dictionary and passed onto `output.py`.
4. `csv_output`contains module `output.py` which generates a CSV file for the overall results with the following information:
    - FASTA Organism
    - FASTA File
    - Strand Searched ('forward' or 'reverse')
    - Motif Enzyme Name
    - Enzyme Source Organism
    - Motif Sequence
    - Significance ('+' for significant, '-' for non-significant, 0 for failed test)
    - p-value (per FASTA file)
    - z-stat (per FASTA file)
    - Observed Motif Count
    - Possible Motif Positions Count
    - Expected Motif Matches Count
    - Expected Motif Probability
    - Genome Length
    - Genome GC Content
    - Observed Genome Base Probabilities (A, C, G, T)

`main.py` executes all modules to produce results. 

`run_test.sh`: Bash script that executes multiple test runs at different coverage threshold values for the user.

`example_runs`: Contains output file directories for different coverage threshold values when `run_test.sh` is executed.

`dependencies.txt`: contains all the neccesary packages for program to execute. 

## Program Instructions

### Installation

1. Download the latest release or clone the repository:
```bash
git clone https://github.com/bluker17/Genome_Motif_Analyzer.git
```

2. Create the conda environment. Conda will automatically create an environment named `motif_analyzer` with all the specified packages and versions required to run the program.
```bash
conda env create -f environment.yml
```

### Usage

1. Activate the environment:
```bash
conda activate motif_analyzer
```
2. Run following command to test:
```bash
testing_materials/run_test.sh
```
3. Example single-line terminal command to execute the program:
```
./main.py -f data/cattle/ -m data/sa_motif.csv -c output/20260520_sa_results.csv -s forward
```

### Command-Line Arguments
| Argument  | Description |
| --- | --- |
| -c, --csv_output | Output CSV file path and filename |
| -f, --fasta_files | Directory containing FASTA files. |
| -m, --motif_file | CSV file containing motif information. |
| -s, --strand_to_search | Strand to search for motifs. Options: forward, reverse, or both. |

### Expected Output:

- Runs the program with multiple coverage thresholds. 
- Prints output locations of each example run.
- Prints completion statement upon successful execution of all runs. 

## References

### Python Standard Library

**`argparse`**    
Python Software Foundation. (2024). *argparse — Parser for command-line options, arguments and sub-commands*. Python 3 Documentation.
https://docs.python.org/3/library/argparse.html

Used for parsing command-line arguments and handling CLI input configuration.

---

**`collections.Counter`**    
Python Software Foundation. (n.d.). *collections — Container datatypes. Python 3 Documentation*.
https://docs.python.org/3/library/collections.html#collections.Counter

Used for counting hashable objects, producing frequency distributions, and efficiently aggregating occurrences in iterable data.

---

**`csv`**    
Python Software Foundation. (n.d.). *csv — CSV File Reading and Writing*. Python 3 Documentation.
https://docs.python.org/3/library/csv.html

Used for reading and writing tabular data in Comma-Separated Values (CSV) format, supporting parsing, dialect handling, and data serialization.

---

**`math`**    
Python Software Foundation. (n.d.). math — Mathematical functions. Python 3 Documentation.
https://docs.python.org/3/library/math.html

Used for mathematical operations such as trigonometry, logarithms, exponentials, rounding, combinatorics, and numeric constants like π and e.

---

**`pathlib`**    
Python Software Foundation. (2024). *pathlib — Object-oriented filesystem paths*. Python 3 Documentation.  
https://docs.python.org/3/library/pathlib.html  

Used for platform-independent file and directory path handling.

---

**`sys`**    
Python Software Foundation. (2024). *sys — System-specific parameters and functions*. Python 3 Documentation.
https://docs.python.org/3/library/sys.html

Used for interacting with interpreter-level functionality such as command-line arguments and program exit handling.

---

**`typing`**    
Python Software Foundation. (n.d.). *typing — Support for type hints. Python 3 Documentation*.
https://docs.python.org/3/library/typing.html

Used for type annotations, including generics, unions, optional types, and static type checking support in Python code.

---

### Third-Party Libraries

**`Biopython`**    
Cock, P. J. A., Antao, T., Chang, J. T., et al. (2009). *Biopython: freely available Python tools for computational molecular biology and bioinformatics*. Bioinformatics, 25(11), 1422–1423.
https://biopython.org/

Used for file parsing (FASTA/GenBank),

---

**`pyahocorasick`**    
Appier Ltd. (n.d.). *pyahocorasick documentation*.
https://pyahocorasick.readthedocs.io/en/latest/

Used for efficient multi-pattern string matching using the Aho–Corasick algorithm, applied to motif searching.

---

**`statsmodels`**    
Seabold, S., & Perktold, J. (2010). *Statsmodels: Econometric and statistical modeling with Python*. Proceedings of the 9th Python in Science Conference.
https://www.statsmodels.org/stable/

Used for 

---

### Bioinformatics Concepts & Standards

**IUPAC Nucleotide Ambiguity Codes**    
Nomenclature Committee of the International Union of Biochemistry (NC-IUB). (1985). Nomenclature for incompletely specified bases in nucleic acid sequences. *European Journal of Biochemistry*, 150(1), 1–5.
https://doi.org/10.1111/j.1432-1033.1985.tb08977.x
The `IUPAC` dictionary maps ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) to their corresponding regex character classes for motif matching.

---

### AI assistance

This project was developed with the help of [ChatGPT-5.5](https://chatgpt.com) (`chatgpt-5.5`) by [OpenAI](https://openai.com).

ChatGPT assisted with:
- Code architecture and implementation
- Docstring and documentation writing
- Debugging and code review

All generated code was reviewed and tested by the author.

