# Genome Motif Analyzer
> Please use the latest stable release. Download here: [Latest Release URL](https://github.com/bluker17/Genome_Motif_Analyzer/releases/latest)

## Author
**Bobby Luker**    
rluker@charlotte.edu    

## Program Description

Given user-provided DNA or RNA FASTA file(s) and a motif information CSV file, the program counts all motif instances within all entries of the FASTA file(s). Utilizing NumPy and Numba vectorized search methods, the motif search can be performed on the forward, reverse, or both strands of the genome, as specified by the user. Following the search, one-sided proportion tests are performed per motif per FASTA entry to determine motif enrichment or depletion. The results are then saved to a CSV file containing motif parameters, FASTA metadata, base probabilities of each FASTA entry, and statistical values used to perform the one-sided proportion test.

Please keep in mind that as the number and size of the FASTA entries and motifs increase, the  computational complexity scales accordingly. Consequently, larger input datasets will prolong the program's execution time.

[GitHub Project URL](https://github.com/bluker17/Genome_Motif_Analyzer)

## License
**GNU General Public License v3.0**

Review `LICENSE` for details.

## Project File Structure:
```
└── 📁Genome_Motif_Analyzer
    └── 📁bin
        ├── main.py
        ├── run_summary.py
    └── 📁data
    └── 📁modules
        ├── run_analyzer.nf
        ├── summarize_stats.nf
    └── 📁output
    └── 📁src
        └── 📁alphabet
            ├── __init__.py
            ├── macromolecule_alphabet.py
            ├── result_alphabet.py
        └── 📁csv_output
            ├── __init__.py
            ├── output.py
        └── 📁file_reader
            ├── __init__.py
            ├── reader.py
        └── 📁motif_locator
            ├── __init__.py
            ├── numba_locator.py
            ├── numpy_locator.py
        └── 📁statistic_analysis
            ├── __init__.py
            ├── statistics.py
        └── 📁summary_statistics
            ├── __init__.py
            ├── summary.py
    └── 📁testing_materials
        └── 📁example_data
            └── 📁DNA_test_phages
            └── 📁RNA_test_phage
            ├── cattle_test_motifs.csv
            ├── DNA_test_motifs.csv
            ├── RNA_test_motifs.csv
        └── 📁example_outputs
            └── 📁cattle_test
            └── 📁DNA_test
            └── 📁RNA_test
        └── 📁expected_example_outputs
            └── 📁cattle_test
                ├── cattle_test_summary_statistics.csv
                ├── cattle_test.csv
            └── 📁DNA_test
                ├── DNA_phage_test_summary_statistics.csv
                ├── DNA_phage_test.csv
            └── 📁RNA_test
                ├── RNA_phage_test_summary_statistics.csv
                ├── RNA_phage_test.csv
    ├── .gitignore
    ├── environment.yml
    ├── LICENSE
    ├── main.nf
    ├── nextflow.config
    └── README.md
```

## Overview
`data/`: An empty directory for the user to place their FASTA and motif data files.

`output/`: An empty directory for the user to place their generated CSV output files from the program.

`bin/`: Contains the two python scripts for the program.
1. `main.py`: Executes most modules to search the FASTA files for motif instances and generate a CSV containing statistical findings upon analysis.
2. `run_summary.py`: Executes the `summary_statistics/` module to generate summary statistics CSV based on the statistical analysis CSV generated from `main.py`.

`src/`: Contains multiple subdirectories leading to modules that handle the input and output data. 

1. `alphabet/` contains modules that define dataclasses containing information passed throughout the program modules.
    - `macromolecule_alphabet.py` defines the DNA and RNA alphabet dataclasses, which contain degenerate and complement base information, and bit maps for converting genetic and motif sequences for faster motif searches. 
    - `result_alphabet.py` has three dataclasses that build upon each other:
        1. `MotifObservations` which stores the search results for a specific motif.
        2. `StrandResults` which stores the base probabilities, and proportion test results per strand of a FASTA entry using the information from `MotifObservations`.
        3. `EntryResults` which stores all gathered results per FASTA entry using `StrandResults`.

2. `file_reader/` contains the module `reader.py` which collects all the provided FASTA files and motifs for searching from the provided motif CSV file.
3. `motif_locator/` contains modules that search the FASTA file(s) for motif instances and pass the results to `statistics.py` via dataclasses found in `result_alphabet.py`.
    - `numpy_locator.py`  searches for all instances of motifs with a length of 4 bp or less in the FASTA files using NumPy vectorization. 
    - `numba_locator.py` searches for all instances of motifs with a length greater than 4 bp in the FASTA files using Numba vectorization. 
4. `statistic_analysis/` contains the module `statistics.py` which performs one-sided proportion tests on motif occurrences per strand and per FASTA file as provided from the `result_alphabet.py` dataclasses. Statistical and locator results are saved into a dictionary and passed onto `output.py`.
5. `csv_output/`contains the module `output.py` which generates a CSV file for the overall results with the following information:
    - FASTA Organism (genetic entry in FASTA file)
    - FASTA File
    - Strand Searched ('forward' or 'reverse')
    - Motif Enzyme Name
    - Enzyme Source Organism
    - Motif Sequence
    - Significance ('+' for observing significantly more motifs than expected, '-' for observing significantly fewer motifs than expected, '0' for no significant difference between observed and expected motifs)
    - p-value (per genetic entry in FASTA file)
    - z-stat (per genetic entry in FASTA file)
    - Observed Motif Count
    - Possible Motif Positions Count
    - Expected Motif Matches Count
    - Expected Motif Probability
    - Genome Length
    - Genome GC Content
    - Observed Base Probabilities (A, C, G, T/U)
6. `summary_statistics/` contains the module `summary.py` which generates a CSV file containing summary statistics based on the CSV generated from `output.py`. The CSV contains the following information:
    - FASTA File
    - Strand
    - Enzyme
    - FASTA Entry Count
    - p-value mean
    - p-value median
    - z-stat mean
    - z-stat median

`modules`: Contains the Nextflow process modules of the program.
1. `run_analyzer.nf` executes `main.py`.
2. `summarize_stats.nf` executes `run_summary.py`.

`main.nf`: The main Nextflow pipeline that executes the Nextflow modules to run the program.

`nextflow.config`: Contains the default testing parameters of the program and establishes  engine profiles that automatically build the necessary program environment.

`testing_materials/example_outputs/`: Contains output CSV files for different parameters and data when a test run is executed.

`testing_materials/expected_example_outputs/`: Contains the expected CSV output files for each test run.

`environment.yml`: Contains the full Conda environment specification, including Python version, required packages, dependency versions, and channels.

`LICENSE`: Contains the full text of the GNU General Public License v3.0, which defines the terms under which the software may be used, modified, and redistributed.

## Program Instructions

### Installation
1. Download the latest release or clone the repository:
```bash
git clone https://github.com/bluker17/Genome_Motif_Analyzer.git
```

2. If Nextflow is not already installed on your machine, then please follow the installation instructions here: [Nextflow Installation](https://docs.seqera.io/nextflow/install)

### Usage
1. Genome Motif Analyzer supports Conda, Micromamba, and Mamba package managers. Execute the default DNA phage test run with the profile that matches your system:
```bash
# If you use standard Conda
nextflow run main.nf -profile conda

# If you use Micromamba
nextflow run main.nf -profile micromamba

# If you use Mamba
nextflow run main.nf -profile mamba
```
**Please keep in mind that the first run will take longer than expected as Nextflow must initialize the environment. After the initial run, Nextflow caches the environment for future runs. Thus, subsequent runs will execute faster.**

If all dependencies are installed locally, then there is no reason to include the -profile parameter when executing runs.

2. Run any of the following commands to further test the program and the time it takes to complete runs with the test FASTA entries and motifs.
```bash
# RNA TEST
# Please update '-profile' with the appropriate profile used on your machine.
nextflow run main.nf -profile ProfileUsed \
    --csv_output testing_materials/example_outputs/RNA_test/RNA_phage_test.csv \
    --fasta_files testing_materials/example_data/RNA_test_phage/ \
    --motif_file testing_materials/example_data/RNA_test_motifs.csv \
    --strand both \
    --macromolecule RNA
```

```bash
# EUKARYOTIC DNA TEST
# Please follow instructions below to download test data before executing.
# Please update '-profile' with the appropriate profile used on your machine.
nextflow run main.nf -profile micromamba \
    --csv_output testing_materials/example_outputs/cattle_test/cattle_test.csv \
    --fasta_files testing_materials/example_data/cattle_test_genome/ \
    --motif_file testing_materials/example_data/cattle_test_motifs.csv \
    --strand both \
    --macromolecule DNA
```
- To execute the 'EUKARYOTIC DNA TEST', please download the "Genome Sequences (FASTA)" option of the [*Bos taurus*](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/) genome from NCBI.

- Save the file as `cattle_test_genome.zip` and unzip the contents in `testing_materials/example_data`. Ensure the contents are located in a directory labeled `cattle_test_genome`. Once complete, the test file can be executed.

3. Results from the test runs can be found in their respective directories in `testing_materials/example_outputs/`. The generated test output CSV files can be compared to the expected results in `testing_materials/expected_example_outputs/`.


4. To run the program with specific data, place the FASTA files in a sub-directory of `data/` and place the motif information CSV file into `data/`. **Please refer to any of the motif CSV files in `testing_materials/example_data/` for the required format of motif information CSV files.** The following is an example of all provided files in the `data/`:
```
└── 📁Genome_Motif_Analyzer
    └── 📁data
        └── 📁cattle
            ├── cattle_1.fasta
            ├── cattle_2.fasta
        ├── sa_motifs.csv
```

5. Following the example in step 4, the program can then be run as the following single-line command:
```bash
nextflow run main.nf -profile ProfileUsed \
    --csv_output output/20260528_cattle_sa_phages.csv \
    --fasta_files data/cattle/ \
    --motif_file data/sa_motifs.csv \
    --strand forward \
    --macromolecule DNA
```

### Command-Line Arguments
| Scope | Argument  | Description |
| --- | --- | --- |
| Python | --csv_output | Output CSV file path and filename |
| Python | --fasta_files | Directory containing FASTA file(s). |
| Python | --motif_file | CSV file containing motif information. |
| Python | --strand | Strand to search for motifs. Options: forward, reverse, or both. |
| Python | --macromolecule | Macromolecule provided in FASTA directory. |
| Nextflow | -profile | Package manager used to establish the program environment |

### Expected Output:

- Searches the provided FASTA file(s) with the provided motif sequence(s).
- Performs a one-sided proportion test per motif per strand of the provided FASTA file(s).
- Results are saved to a CSV file.
- Summary statistics are generated into a CSV file based on the initial CSV file.

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

Used for counting nucleotide bases within a FASTA entry.

---

**`csv`**    
Python Software Foundation. (n.d.). *csv — CSV File Reading and Writing*. Python 3 Documentation.
https://docs.python.org/3/library/csv.html

Used for reading and writing data in Comma-Separated Values (CSV) format.

---

**`math`**    
Python Software Foundation. (n.d.). *math — Mathematical functions*. Python 3 Documentation.
https://docs.python.org/3/library/math.html

Used for establishing Not a Number constants for edge cases.

---

**`pathlib`**    
Python Software Foundation. (2024). *pathlib — Object-oriented filesystem paths*. Python 3 Documentation.  
https://docs.python.org/3/library/pathlib.html  

Used for platform-independent file and directory path handling.

---

**`sys`**    
Python Software Foundation. (2024). *sys — System-specific parameters and functions*. Python 3 Documentation.
https://docs.python.org/3/library/sys.html

Used for interacting with command-line arguments and program exit handling.

---

**`traceback`**    
Python Software Foundation. (n.d.). *traceback — Print or retrieve a stack traceback. Python 3 Documentation*.
https://docs.python.org/3/library/traceback.html

Used for extracting, formatting, and printing the call stack of an exception, enabling detailed debugging and logging of errors in Python applications.

---

**`typing`**    
Python Software Foundation. (n.d.). *typing — Support for type hints. Python 3 Documentation*.
https://docs.python.org/3/library/typing.html

Used for type annotations in Python code.

---

### Third-Party Python Libraries

**`Biopython`**    
Cock, P. J. A., Antao, T., Chang, J. T., et al. (2009). *Biopython: freely available Python tools for computational molecular biology and bioinformatics*. Bioinformatics, 25(11), 1422–1423.
https://biopython.org/

Used for parsing FASTA files to extract sequences from FASTA entries. 

---

**`statsmodels`**    
Seabold, S., & Perktold, J. (2010). *Statsmodels: Econometric and statistical modeling with Python*. Proceedings of the 9th Python in Science Conference.
https://www.statsmodels.org/stable/

Used to calculate a one-sided proportion test per motif per FASTA entry.

---

### Bioinformatics Concepts & Standards

**IUPAC Nucleotide Ambiguity Codes**    
Nomenclature Committee of the International Union of Biochemistry (NC-IUB). (1985). Nomenclature for incompletely specified bases in nucleic acid sequences. *European Journal of Biochemistry*, 150(1), 1–5.
https://doi.org/10.1111/j.1432-1033.1985.tb08977.x
The `IUPAC` dictionary maps ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) to their corresponding regex character classes for motif matching.

---

**Nextflow Workflow Management**    
P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:10.1038/nbt.3820

---

### AI assistance

This project was developed with the help of [ChatGPT-5.5](https://chatgpt.com) by [OpenAI](https://openai.com) and [Gemini-3.1-Pro](https://gemini.google.com/app) by [Google](https://gemini.google/about/).

ChatGPT assisted with:
- Code architecture and implementation of Numba search methods
- Code architecture and implementation of bitmask conversion for search methods
- Docstring and documentation writing

Gemini assisted with:
- Code architecture and implementation of Alphabet dataclasses and Nextflow pipeline
- Pipeline orchestration debugging and configuration code review
- Debugging and code review

All generated code was reviewed and tested by the author.