# Genome Motif Analyzer
> Please use the latest stable release. Download here: [Latest Release URL](hhttps://github.com/bluker17/Genome_Motif_Analyzer/releases/latest)

## Author
**Bobby Luker**    
rluker@charlotte.edu    
UNCC ID: 801484356

## Program Description

From user provided DNA or RNA FASTA file(s) and a motif information CSV file, the program counts all motif instances within the FASTA file(s) using the Aho-Corasick method. One-sided proportion tests are then performed per motif per genome to determine motif enrichment or depletion. These results can reveal evolutionary genetic changes in an organism over time. The results are saved to a CSV containing motif information, FASTA file information, base probabilities of a FASTA file, and statistic values used to perform the one-sided proportion test.

The program handles eukaryotic and prokaryotic FASTA files. When running the program with a eukaryotic organism, the motif search and statisic calculations are performed per chromosome.

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
`data/`: Contains the example genome FASTA files and example format of motif information CSV files.

`output/`: Contains the generated CSV output files from running the program.

`src/`: Contains multiple subdirectories leading to modules handling the input and output data. 
1. `file_reader/` contains module `reader.py` which collects all the provided FASTA files and motifs to search from the provided motif CSV file.
2. `motif_locator/` contains module `locator.py` which searches for all motif instances in the FASTA files using the Aho-Corasick algorithm. Results are passed onto the `statistics.py` for further analysis in the form of a dictionary.
3. `statistic_analysis/` contain module `statistics.py` which performs a one-sided proportion tests on motif occurrences per strand and per FASTA file as provided from the `locator.py` findings. Statistical and locator results are saved into a dictionary and passed onto `output.py`.
4. `csv_output/`contains module `output.py` which generates a CSV file for the overall results with the following information:
    - FASTA Organism (genetic entry in FASTA file)
    - FASTA File
    - Strand Searched ('forward' or 'reverse')
    - Motif Enzyme Name
    - Enzyme Source Organism
    - Motif Sequence
    - Significance ('+' for observing signficantly more motifs than expected, '-' for observing significanty less motifs than expected, 0 no signifcant difference of motifs from expected)
    - p-value (per genetic entry in FASTA file)
    - z-stat (per genetic entry in FASTA file)
    - Observed Motif Count
    - Possible Motif Positions Count
    - Expected Motif Matches Count
    - Expected Motif Probability
    - Genome Length
    - Genome GC Content
    - Observed Genome Base Probabilities (A, C, G, T)

`main.py` executes all modules to produce results. 

`run_test.sh`: Bash script that executes multiple test runs at different coverage threshold values for the user.

`testing_materials/example_runs/`: Contains output CSV files for different parameters and data when a `run_test_X.sh` is executed.

`testing_materials/expected_example_runs/`: Contains the expected CSV output files for each `run_test_X.sh`.

`environment.yml`: Contains the full Conda environment specification, including Python version, required packages, dependency versions, and channels.

`dependencies.txt`: Contains all the neccesary packages for program to execute. 

`LICENSE`: Contains the full text of the GNU General Public License v3.0, which defines the terms under which the software may be used, modified, and redistributed.

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
2. Run any or all of the following commands to test the program:
```bash
testing_materials/run_test_1.sh
testing_materials/run_test_2.sh
testing_materials/run_test_3.sh
testing_materials/run_test_4.sh
```
3. Results from `testing_materials/run_test_X.sh` can be found in `testing_materials/example_outputs/`. The generated test output CSV files can be compared to the expected results in `testing_materials/expected_example_outputs/`.

4. To run the program with specific data, place the FASTA files in a sub-directory of `data/` and place the motif information CSV file into `data/`. **Please refer to `testing_materials/example_data/` for the required format of motif information CSV files.** The following is an example of all provided files in the `data/`:
```
└── 📁Genome_Motif_Analyzer
    └── 📁data
        └── 📁cattle
        ├── sa_motifs.csv
```

5. Following the example in 4, the program can then be run as the following single-line command:
```bash
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

- Searches the provided FASTA file(s) with the provided motif sequence(s).
- Performs a one-sided proportion test per motif per strand of the provided FASTA file(s).
- Results saved to a CSV file.

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

**`collections.defaultdict`**    
Python Software Foundation. (n.d.). *collections — Container datatypes. Python 3 Documentation*.
https://docs.python.org/3/library/collections.html#collections.defaultdict

Used for dictionary-like objects that automatically initialize missing keys with a default factory function (e.g., lists, integers, sets), commonly used to simplify grouping and aggregation logic without needing explicit key existence checks.

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