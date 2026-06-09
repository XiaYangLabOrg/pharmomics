# Running app3_analysis_local.R

This script performs Jaccard similarity analysis between disease genes and drug signatures across multiple species (Human, Mouse, Rat).

## Prerequisites

- R (version 3.6 or higher recommended)
- Required R packages:
  - `GeneOverlap`

Install the required package if needed:
```r
install.packages("GeneOverlap")
```

## Required Input Files

The script expects the following files to be present in the working directory:

**Note:** The required gene input files depend on the `consider_direction` setting:

- **When `consider_direction = TRUE`**: Requires `up_genes.txt` (required) and `down_genes.txt` (optional, can be empty but must have GENE column if present)
- **When `consider_direction = FALSE`**: Requires `genes.txt` (single file with all genes)

1. **`up_genes.txt`** - Tab-delimited file containing upregulated genes (required when `consider_direction = TRUE`)
   - Must have a column named `GENE`
   
2. **`down_genes.txt`** - Tab-delimited file containing downregulated genes (optional when `consider_direction = TRUE`)
   - Must have a column named `GENE`
   - When `consider_direction = TRUE`, this file is optional and can be empty (just the header with GENE column)
   - When `consider_direction = FALSE`, this file is not used (the script loads `genes.txt` instead)
   
3. **`genes.txt`** - Tab-delimited file containing all genes (required when `consider_direction = FALSE`)
   - Must have a column named `GENE`
   - Used for both up and down gene lists when direction is not considered

4. **`Jaccard_app_databasev3_No_GeoDE.rda`** - R data file containing:
   - `Humanframe` - Human drug signature data
   - `Ratframe` - Rat drug signature data
   - `Mouseframe` - Mouse drug signature data
   - `Humangenesup` - Human upregulated gene lists
   - `Humangenesdown` - Human downregulated gene lists
   - `Ratgenesup` - Rat upregulated gene lists
   - `Ratgenesdown` - Rat downregulated gene lists
   - `Mousegenesup` - Mouse upregulated gene lists
   - `Mousegenesdown` - Mouse downregulated gene lists
   - `HUGO_symbols2` - Human gene symbol mapping
   - `RAT_symbols2` - Rat gene symbol mapping
   - `Mouse_symbols2` - Mouse gene symbol mapping

5. **`ADR_scores.rds`** - RDS file containing ADR (Adverse Drug Reaction) scores
   - Must contain `hepatotox_complete_signature_Jaccard` field

6. **`hepatotox_genes.rda`** - R data file containing hepatotoxicity gene signatures

## Command Line Arguments

The script requires command line arguments:

### Required Arguments

1. **`consider_direction`** (first argument, mandatory)
   - Whether to consider gene direction (up vs down) in the analysis
   - Acceptable values: `TRUE`, `FALSE`, `T`, `F`, `1`, `0`
   - If `TRUE`: Loads `up_genes.txt` (required) and `down_genes.txt` (optional, can be empty) and uses directional Jaccard score calculation (up-up + down-down - up-down - down-up)
   - If `FALSE`: Loads `genes.txt` and uses standard Jaccard score (ignores direction)

### Optional Arguments

2. **`GENECOUNT_THRESHOLD`** (second argument, optional)
   - Minimum number of genes required in a signature to be included in analysis
   - Default: `40`
   - Must be a numeric value

## Usage Examples

### Basic Usage

**With direction consideration (default threshold):**
```bash
Rscript app3_analysis_local.R TRUE
```

**Without direction consideration (default threshold):**
```bash
Rscript app3_analysis_local.R FALSE
```

### Custom Threshold

**With direction consideration and custom threshold:**
```bash
Rscript app3_analysis_local.R TRUE 50
```

**Without direction consideration and custom threshold:**
```bash
Rscript app3_analysis_local.R FALSE 30
```

### Alternative Boolean Formats

The script accepts various formats for the `consider_direction` argument:
```bash
# All of these are equivalent:
Rscript app3_analysis_local.R TRUE
Rscript app3_analysis_local.R T
Rscript app3_analysis_local.R 1

# All of these are equivalent:
Rscript app3_analysis_local.R FALSE
Rscript app3_analysis_local.R F
Rscript app3_analysis_local.R 0
```

## Output Files

The script generates two output files:

1. **`_app3result.txt`** - Main results file (tab-delimited) containing:
   - Database information
   - Method, Drug, Species, Tissue/Cell Line
   - Study details (Dose, Time)
   - Jaccard Score, Odds Ratio, P value
   - Within Species Rank
   - Gene overlap information
   - Signature details

2. **`_app3result_hepatotox.txt`** - Hepatotoxicity analysis results (tab-delimited) containing:
   - Adverse drug reaction name
   - Jaccard score
   - Jaccard score rank
   - Jaccard p value
   - ADR genes and input gene overlap

## What the Script Does

1. **Loads input data**: Reads gene lists and required database files
2. **Filters signatures**: Removes signatures with fewer genes than `GENECOUNT_THRESHOLD`
3. **Species detection**: Automatically detects if input genes are from Human or non-Human species
4. **Gene mapping**: Maps genes across species (Human ↔ Rat ↔ Mouse) as needed
5. **Jaccard analysis**: Calculates Jaccard similarity scores for each drug signature:
   - If `consider_direction = TRUE`: Uses directional Jaccard score (up-up + down-down - up-down - down-up)
   - If `consider_direction = FALSE`: Uses standard Jaccard score
6. **Statistical analysis**: Calculates odds ratios and p-values using GeneOverlap
7. **Ranking**: Ranks results within each species
8. **Hepatotoxicity assessment**: Compares input genes against hepatotoxicity signatures
9. **Output generation**: Writes results to tab-delimited files

## Progress Indicators

The script prints progress percentages during execution:
- Starts at 5%
- Updates as it processes Rat, Mouse, and Human frames
- Completes at 100%

## Error Handling

The script will stop with an error message if:
- `consider_direction` argument is missing
- `consider_direction` argument is not a valid boolean value
- Required input files are missing
- Required R packages are not installed

## Notes

- The script processes data for all three species (Human, Mouse, Rat) sequentially
- Results are sorted by Jaccard score in descending order
- Empty gene lists are handled gracefully
- The script includes debug output for Rat species when `consider_direction = TRUE`

