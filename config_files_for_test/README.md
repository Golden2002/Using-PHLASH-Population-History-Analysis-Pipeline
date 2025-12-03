
# PHLASH Test Configuration Module

This directory provides example configuration files (`config.sh.example`) for testing the PHLASH analysis pipeline, covering both autosomes and sex/mitochondrial chromosomes. Users can copy and modify it to suit their environment.

## Quick Start

1. **Copy the example config**

   ```bash
   cp config.sh.example config.sh
   ```

2. **Edit `config.sh`**

   ### Required settings

   * **VCF files**

     * For autosomes: `VCF` (phased VCF, e.g., combined chr1–chr22)
     * For sex/mitochondrial chromosomes: `VCF_X`, `VCF_Y`, `VCF_MT`
   * **Sample information file**: `INFO_FILE` (two columns: `SampleID` `PopulationID`)
   * **Male sample IDs file** (if using sex chromosomes): `MALE_SAMPLES_FILE`
   * **PHLASH environment**: `PHLASH_ENV` (path to Conda or virtual environment activate script)
   > Notice that the coordinates in my scripts are based on GRCh38 reference genome, where the Chromosomes are named as follows:
   > Autosomes: chr1, chr2, …, chr22 (or chr1_23 if you combine all autosomes)
   > Sex chromosomes: chrX, chrY
   > Mitochondrial genome: chrM

   ### Optional settings

   * Output directory: `RESULTS_DIR` (default: `results`)
   * Analysis parameters:

     * `MUTATION_RATE` (default: 1.25e-8)
     * `KNOTS` (default: 20)
     * `REG_PENALTY` (default: 6.0)
     * `SPLINE_TYPE` (default: piecewise)
     * `SEED` (random seed file)
     * `SAMPLE_PER_POP` (number of samples per population, default: 10)
     * `MIN_CONTIGS` (minimum number of contigs, default: 50)
     * `MAX_CHROMOSOMES` (max autosomes to analyze, default: 22)
     * `MASK_FILE` (optional BED mask file)

3. **Run the test pipeline**
   Submit the SLURM job using your modified `config.sh`:

   ```bash
   sbatch run_phlash.sh
   ```

## Notes

* This configuration is for testing and demonstration purposes.
* Ensure all paths exist and files are correctly formatted.
* The pipeline starts directly from the specified VCF files and is fully compatible with SLURM HPC environments.
* You can customize parameters to match your dataset size and analysis goals.

---
