# README.md
### Step 1
Activating a virtual environment (venv)
A virtual environment is a separate Python environment, where you can install packages without affecting the global Python installation on your system.

To activate a venv on Windows:

```
cd myproject
path\to\python -m venv venv
venv\Scripts\activate
```

To activate a venv on Unix or macOS:
```
cd myproject
python3 -m venv venv
source venv/bin/activate
```

### Step 2
Installing requirements from requirements.txt
To install the required packages for your project, you can use the following command

```
pip install -r requirements.txt
```
This command reads the requirements.txt file and installs the listed packages and their dependencies in the active virtual environment.


### Step 3

Important input files to have inorder to run the full main.py pipline:
1. 'input_file': "data/raw_data/sequencing_output_fastq/all_sequencing_output.fastq".
   
   all_sequencing_output.fastq file contains all the sequences after combining all the sequencing output fastq files.  
2. 'design_before_conversion_file': "data/raw_data/design/all_design_before_parse.fasta".
   
   all_design_before_parse.fasta file contains all the sequences after combining all the sequencing output fasta files. 
3. 'const_design_file': "config/design.csv".
   design.csv file contains the design and location of the universals in the sequences.
4. 'barcodes_design_file': "config/barcodes_design.csv".

   barcodes_design.csv file contains the design and location of the barcodes in the sequences.
5. 'payload_design_file': "config/payload_design.csv".

   payload_design.csv file contains the design and location of the payloads in the sequences.

Note:
Make sure that the config files are with the correct locations in the sequence.

### Step 4
Run pipline:
```
cd dna_storage_helix_experiment\analyze_sequencing_data
import main
python main.main()
```