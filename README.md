# RP\_blast

This repository contains a Perl script for running `tblastn` searches to identify and process homologous genes across multiple genomes. The tool extracts genes of interest, performs BLAST searches, and organizes results into meaningful alignments and ortholog groups.

---

## Features

- **Gene BLAST Search:** Searches for homologous genes using `tblastn`.
- **Customizable Parameters:** Adjustable e-value thresholds, query coverage, and identity thresholds.
- **Ortholog Identification:** Groups homologous hits and extracts orthologous sequences.
- **Sequence Alignment:** Prepares aligned CDS and protein sequences for downstream analyses.

---

## Requirements

### Dependencies

- Perl 5 or later
- `tblastn` (part of the BLAST+ suite)

Install Perl modules via CPAN if necessary.

### Inputs

1. **Gene list file (********`genes_list.txt`********)**: Contains genes of interest (names and accessions).
2. **Drosophila CDS FASTA file (********`dmel_cds.fas`********)**: Reference coding sequences.
3. **Genome directory (********`path_to_genomes`********)**: Directory containing genome FASTA files for target species.

---

## Usage

### Example Command

```bash
perl RP_blast.pl genes_list.txt dmel_cds.fas path_to_genomes
```

### Options

| Option          | Description                               | Default |
| --------------- | ----------------------------------------- | ------- |
| `max_eval`      | Maximum e-value threshold for BLAST hits. | 3       |
| `min_qcov`      | Minimum query coverage for BLAST hits (%) | 50      |
| `min_hit_ident` | Minimum identity threshold for hits (%)   | 0       |
| `min_grp_ident` | Minimum identity threshold for groups (%) | 0       |

---

## Output

### Generated Files

1. **Raw BLAST Results (********`CDS_seqs/Raw_blast`********)**
   - Contains raw `tblastn` outputs for each gene-genome combination.
2. **Ortholog Sequences (********`CDS_seqs/Orthologs`********)**
   - Filtered and grouped orthologous CDS and protein sequences.
3. **Alignments (********`CDS_seqs/Alignments`********)**
   - Aligned sequences for homologous genes across species.
4. **Summary Table (********`homologous_genes.txt`********)**
   - Comprehensive summary of homologous genes and their alignment statistics.

---

## Notes

- **BLAST Configuration:** Adjust `tblastn` options directly in the script via `$tbnopt`.
  ```perl
  my $tbnopt = '-seg no -word_size 2';
  ```
- Ensure all required genome files are in the specified directory and follow consistent naming conventions.
- Outputs will overwrite existing files with the same names.

---

## License

This project is licensed under the MIT License.

---

## Author

Daniel Gebert\
dg572\@cam.ac.uk

Feel free to open an issue or contact me with questions or suggestions!

