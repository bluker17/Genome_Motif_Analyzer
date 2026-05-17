# Program Title
> Please use the latest stable release. Download here: [Latest Release URL](https://github.com/bluker17/Grid_Overlay/releases/latest)

## Author:
**Bobby Luker**
rluker@charlotte.edu
UNCC ID: 801484356

## Program Description:

For a draft genome assembly, when given multiple BLAST result TSV files and a contig TXT file containing contig IDs and lengths, the contigs are filtered by bin priority, bitscore, coverage threshold, and contig size threshold. Summary statistics are generated for each bin, along with two bar plots showing the number of contigs per bin and the total base pairs per bin. The default coverage threshold is set to 0.9 to ensure that only the highest-quality contigs are used. However, the user can adjust the coverage threshold to relax this restriction if needed.


[Github Project URL](https://github.com/bluker17/Template_repository)

## License: 
**MIT License**

Review `License` for details.

## Project File Structure:
***envrionment.yml files need to be created***
```
└── 📁testing_materials
    └── 📁example_data
        ├── data.txt
    └── 📁example_outputs
        ├── test_output.txt
    ├── run_test.sh
└── 📁output
    ├── output.txt
└── 📁src
├── environment-alternative.yml
├── environment.yml
├── LICENSE
├── main.py
├── README.md
├── dependencies.txt
```

## Overview:
`example_data`: Contains the example BLAST result TSV files and the contig size TXT file.

`src`: Contains multiple subdirectories leading to modules handling the input and output data. 
1. `reader` contains module `read_files.py` which reads in the BLAST result TSV files and contig size TXT file to merge all files into one main data frame for further analysis.
2. `prioritization` contains module `priority.py` which adds a bin column to the data frame. The data frame is then filtered by the highest priority bin and bitscore for each contig id. 

`main.py` executes all modules to produce results. 

`run_test.sh`: Bash script that executes multiple test runs at different coverage threshold values for the user.

`example_runs`: Contains output file directories for different coverage threshold values when `run_test.sh` is executed.

`dependencies.txt`: contains all the neccesary packages for proram to execute. 

## Program Instructions:

### Installation

1. Download the latest release or clone the repository:
```bash
git clone https://github.com/bluker17/Grid_Overlay.git
```

2. Create the conda environment. Conda will automatically create an environment named `grid-app` with all the specified packages and versions.
```bash
conda env create -f environment.yml
```

### Usage

1. Activate the environment:
```bash
conda activate env_name_to_be_created
```
2. Run following command to test:
```bash
bash testing_materials/run_test.sh
```
3. Example single-line terminal command to execute the program:
```
python3 main.py -p data/data.txt --coverage_threshold 0.75 --contig_size_threshold 1000
```

#### Command-Line Arguments:
| Argument                          | Description                                  | Default               |
| --------------------------------- | -------------------------------------------- | --------------------- |
| `-p`, `--priority_file` | TSV file containing bins and priority. 'Bin' corresponds to files assigned to contigs. 'priority' corresponds to numerical priority classification for each bin. The lowest numerical value roeesponds to the highest priority bin.                           | example_data/prioritization.tsv|
| `--coverage_threshold` | Coverage threshold float. Used to filter out contigs less than the threshold.                           | 0.9|
| `--contig_size_threshold` | Coverage threshold int. Used to filter out contigs less than a specific bp size.                           | 3000|




Expected Output:

- Runs the program with multiple coverage thresholds. 
- Prints output locations of each example run.
- Prints completion statement upon successful execution of all runs. 


## Contributions (if there are multiple authors)


## References:

### Python Standard Library

**`argparse`**  
Python Software Foundation. (2024). *argparse — Parser for command-line options, arguments and sub-commands*. Python 3 Documentation.  
https://docs.python.org/3/library/argparse.html  

Used for parsing command-line arguments and handling CLI input configuration.

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

### Third-Party Libraries
#### Pillow (PIL)

**`Pillow`**  
Pillow Contributors. (n.d.). *Pillow (PIL Fork) documentation*. Pillow Documentation.  
https://pillow.readthedocs.io/en/stable/  

Used for image loading, manipulation, drawing operations, colour handling, and alpha compositing.

---

**`ImageColor`**  
Pillow Contributors. (n.d.). *ImageColor module*. Pillow Documentation.  
https://pillow.readthedocs.io/en/stable/reference/ImageColor.html  

Used for parsing and converting colour specifications into RGB values.

---

**`ImageColor.getrgb`**  
Pillow Contributors. (n.d.). *PIL.ImageColor.getrgb*. Pillow Documentation.  
https://pillow.readthedocs.io/en/stable/reference/ImageColor.html#PIL.ImageColor.getrgb  

Used to convert colour strings and hexadecimal colour codes into RGB tuples.

---

**`Image.open`**  
Pillow Contributors. (n.d.). *PIL.Image.open*. Pillow Documentation.  
https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.open  

Used to load and open image files for processing.

---

**`Image.convert`**  
Pillow Contributors. (n.d.). *PIL.Image.Image.convert*. Pillow Documentation.  
https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.Image.convert  

Used to convert images between colour modes (e.g., RGB, RGBA, grayscale).

---

**`Image.new`**  
Pillow Contributors. (n.d.). *PIL.Image.new*. Pillow Documentation.  
https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.new  

Used to create new blank image canvases with specified dimensions and colour modes.

---

**`ImageDraw.Draw`**  
Pillow Contributors. (n.d.). *PIL.ImageDraw.Draw*. Pillow Documentation.  
https://pillow.readthedocs.io/en/stable/reference/ImageDraw.html#PIL.ImageDraw.Draw  

Used to create drawable image contexts for rendering shapes, lines, and annotations.

---

**`Image.alpha_composite`**  
Pillow Contributors. (n.d.). *PIL.Image.alpha_composite*. Pillow Documentation.  
https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.alpha_composite  

Used to combine RGBA images using alpha-channel compositing.

### Bioinformatics Concepts & Standards

**IUPAC Nucleotide Ambiguity Codes**
Nomenclature Committee of the International Union of Biochemistry (NC-IUB). (1985). Nomenclature for incompletely specified bases in nucleic acid sequences. *European Journal of Biochemistry*, 150(1), 1–5.
https://doi.org/10.1111/j.1432-1033.1985.tb08977.x
The `IUPAC` dictionary maps ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) to their corresponding regex character classes for motif matching.
 
**Restriction Enzyme Cut Notation**
Rebase - The Restriction Enzyme Database. Roberts, R. J., Vincze, T., Posfai, J., & Macelis, D. (2023). REBASE: a database for DNA restriction and modification: enzymes, genes and genomes. *Nucleic Acids Research*, 51(D1), D629–D630.
https://doi.org/10.1093/nar/gkac975
https://rebase.neb.com/
The `^` (top-strand cut) and `_` (bottom-strand cut) notation parsed by `_top_cut_offset()` and `_bot_cut_offset()` follows the REBASE convention for describing staggered and blunt restriction enzyme cleavage sites.
 
**DNA Complementarity**
Watson, J. D., & Crick, F. H. C. (1953). Molecular structure of nucleic acids: A structure for deoxyribose nucleic acid. *Nature*, 171, 737–738.
https://doi.org/10.1038/171737a0
The `COMPLEMENT` dictionary (A↔T, G↔C) used to generate the bottom strand in `DoubleStrandedMap` is based on Watson–Crick base-pairing rules.



### AI assistance:

This project was developed with the help of [ChatGPT-5.5](https://chatgpt.com) (`chatgpt-5.5`) by [OpenAI](https://openai.com).

ChatGPT assisted with:
- Code architecture and implementation
- Docstring and documentation writing
- Debugging and code review

All generated code was reviewed and tested by the author.

