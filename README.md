This Python script is a comprehensive tool designed to extract, process, and analyze genomic features from GenBank files. It provides advanced functionality for parsing GenBank files, extracting specific features (e.g., coding sequences (CDS), tRNA, rRNA, genes, or mRNA), and generating detailed outputs in both FASTA and HTML formats. Below is a description of its key components and functionality:

---

### **Script Description: Advanced GenBank Feature Extractor**

#### **Purpose**
The script automates the extraction of genomic features from GenBank files, enabling researchers to focus on specific types of features (e.g., CDS, tRNA) or gene names. It also generates detailed reports and exports data in standardized formats for downstream bioinformatics analyses.

---

#### **Key Features**
1. **Feature Extraction**:
   - Extracts genomic features such as CDS, tRNA, rRNA, genes, or mRNA from GenBank files.
   - Supports filtering by specific gene or product names, allowing users to target particular genes or products of interest.

2. **Duplicate Handling**:
   - Automatically skips duplicate sequences based on unique IDs, ensuring clean and non-redundant datasets.

3. **Sanitization of IDs**:
   - Cleans sequence IDs to ensure compatibility with phylogenetic software by removing illegal characters (e.g., spaces, colons, parentheses).

4. **Metadata Enrichment**:
   - Captures rich metadata, including species information, gene names, product descriptions, and sequence lengths, which are embedded in the output files.

5. **Output Formats**:
   - **FASTA Files**: Exports extracted sequences into separate FASTA files for each gene, organized in a user-specified output directory.
   - **HTML Report**: Generates a detailed, interactive HTML report summarizing:
     - Gene statistics (e.g., most common genes, total sequences per gene).
     - Species distribution (e.g., number of sequences per species).
     - A complete table of all extracted sequences with metadata.

6. **Error Handling**:
   - Includes robust error handling to gracefully manage issues during file processing and feature extraction.

7. **Command-Line Interface**:
   - Provides a user-friendly command-line interface with options to specify:
     - Input path (single file or directory containing multiple GenBank files).
     - Type of feature to extract (e.g., CDS, tRNA).
     - Specific gene or product names to filter.
     - Output directory for results.

---

#### **Workflow Overview**
1. **Input Parsing**:
   - Accepts a single GenBank file or a directory containing multiple GenBank files.
   - Reads and processes each file using the `Bio.SeqIO` module from Biopython.

2. **Feature Extraction**:
   - Iterates through genomic features in each file, identifying those matching the specified type (e.g., CDS) and optional gene/product names.
   - Extracts relevant sequences and associated metadata (e.g., gene name, product description, species).

3. **Data Processing**:
   - Organizes extracted features into a structured format, grouping sequences by gene name.
   - Tracks statistics such as the total number of unique genes, species, and sequences.

4. **Output Generation**:
   - Writes sequences to FASTA files, one per gene, with sanitized IDs and descriptive metadata.
   - Creates an interactive HTML report summarizing key findings, including gene and species statistics, and a table of all extracted sequences.

---

#### **Usage Example**
To extract coding sequences (CDS) for specific genes (e.g., "COX1", "ATP6") from a directory of GenBank files and save the results in a folder named `output_dir`, run the following command:

```bash
python genbankparser.py /path/to/genbank_files --feature_type CDS --feature_names COX1 ATP6 --output_dir output_dir
```

---

#### **Dependencies**
- **Python Libraries**:
  - `argparse`: For command-line argument parsing.
  - `os`: For file and directory operations.
  - `collections.defaultdict`: For efficient data organization.
  - `Bio.SeqIO` and `Bio.SeqRecord` (from Biopython): For parsing and handling GenBank files.

---

#### **Applications**
This script is particularly useful for:
- Phylogenetic studies requiring curated sequence datasets.
- Comparative genomics research focusing on specific genes or features.
- Generating summaries of genomic data for publication or further analysis.

---

By combining powerful feature extraction capabilities with intuitive reporting, this script serves as a versatile tool for genomic data analysis and exploration.
