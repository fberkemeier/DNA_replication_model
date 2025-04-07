# A Model of DNA Replication Timing

## Description

We introduce a comprehensive toolkit for analysing DNA replication timing, origin firing rates, and genomic stability across cell lines and chromosomal regions. Complementing the work in [Berkemeier et al. (2024)](https://www.biorxiv.org/content/10.1101/2024.11.25.625090v2), it includes functions for loading and processing data across whole-genome regions, telomeres, centromeres, and specific loci of interest. By fitting origin firing rates to replication timing data, the toolkit efficiently predicts and compares experimental and modelled timing profiles. The resulting error distributions between predicted and experimental data help pinpoint regions of interest. With datasets spanning diverse chromosomes and genomic features, this toolkit enables detailed visualisation and analysis of replication dynamics and genomic stability.

<!---
⚠️ **Note:** This repository is a work in progress and represents an ongoing project. The code and documentation are subject to updates and refinements, and while we strive for accuracy, they may not yet reflect the final, polished version. Datasets for minimal working examples are currently stored in this repository via Git Large File Storage (LFS) but will soon be relocated to an external hosting platform for better accessibility. Your understanding and feedback are appreciated.
-->

## Mathematical model

Consider a DNA molecule with $n$ discrete genomic loci, where each locus $j$ can potentially act as an origin that fires at rate $f_j$ (indepedently of other origins) to initiate a fork that progresses bidirectionally with speed $v$, typically measured in kilobases per minute (kb/min). Let $T_j$ be the time a site $j$ takes to fire or be passively replicated by a fork. Then, for a sufficiently large genome (such as human), the expected time of replication at $j$ is given by

$$
\mathbb{E}[T_j]=\sum_{k=0}^{\infty}  \frac{e^{-\sum_{|i|\leq k}(k-|i|)f_{j+i}/v}-e^{-\sum_{|i|\leq k}(k+1-|i|)f_{j+i}/v}}{\sum_{|i|\leq k} f_{j+i}}.
$$

While this expression holds true for an infinitely large genome, in practical terms this series can be limited to $0\leq k\leq  R<n/2$, for some large enough $R$. This parameter represents the radius of replication influence: the distance within which neighbouring origins $\{j-R,...,j-1,j+1,...,j+R\}$ are assumed to affect the timing of a focal origin $j$. In other words, while every firing origin does theoretically affect replication timing at any other location, this effect decays rapidly with distance from the origin of interest $j$. Numerically, the finite version of the equation should mimic the average replication timing obtained from computational simulations.

## Usage

Clone the repository to a local directory. To explore the functionality and understand the workflow, we recommend using the provided dataset (`data.zip`), which includes all necessary dependencies. Extract it to the main directory. Note that datasets for minimal working examples are currently stored in this repository via Git Large File Storage (LFS), and `data.zip` must be downloaded and extracted.

Once extracted, place the files into a subfolder named `data` within the main project directory. The dataset includes `.bedgraph` files representing timing errors and origin firing rates. These files can be uploaded to the Genome Browser for visualisation and comparison with other genomic data. For example, one should be able to access the bedgraph files for the error misfits at `data/whole-genome_error/bedgraph_files/error_HUVEC.bedgraph`, from the directory where `plots.ipynb` is located.

In all subsequent analyses, first execute the code below as a prerequisite.
```
from ipynb.fs.full.utilities import *
```

### 1. Data generation

If you want to upload your own data, our tool allows you to process timing data and convert it to origin firing rates using our model. This can be run locally or on an HPC platform. To handle bigWig files, we recommend installing [`pybigtools`](https://pypi.org/project/pybigtools/). All required utilities and dependencies are documented in `utilities.ipynb`. DNA replication simulations are performed using the [Beacon Calculus](https://github.com/MBoemo/bcs).

#### Examples

As an example, consider importing a bigWig file for HUVEC cells from the [ENCODE database](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeUwRepliSeq) (wavelet smooth signal). A local copy is also provided at `data/bigwig_files/HUVEC.bw`. Before fitting the data, it must be pre-processed. This can be done using `model.ipynb` (Data generation), which converts Repli-seq bigWig files into text files at the desired resolution (e.g., 1 kb). The processed files are stored in `data/whole-genome_timing_data`.

The script `code_fit.py` contains the main fitting and data generation functions:
- **`fitfunction`**: Performs model fitting.
- **`datagenfs`**: Processes the input data for fitting.

You can customise the fitting process by modifying the following parameters:
- `cell_line`: Specify the cell line name.
- `chr_number`: Select the chromosome number.
- `chrpos_min` and `chrpos_max`: Define the fitting region (or enable whole-genome fitting with `all_dataQ`—note this is computationally intensive).
- Additional options:
  - **Fork speed**: `fork_speed`
  - **Sampling intervals**: `resolution` (default: 1000 for 1 kb)
  - **Timing data scaling**: `scale-factor`
  - **Fitting iterations**
  - **Radius of influence**: `int_width` (in bp)
 
The simplest way to understand how our fitting algorithm works is to follow the example in `model.ipynb`. To generate the timing data files (`*.txt`) for a given cell line, place its bigWig file in the directory `data/bigwig_files/{cell_line}.bw`. If you already have `data/whole-genome_timing_data` and `data/whole-genome_missing_data` for that cell line, you can skip this step. Ensure the `data` folder and any required subdirectories exist before proceeding. Note that `spec_fileQ` refers to fitting a specific region. For a whole-genome fit (`all_dataQ = True`), HPC resources are recommended, and a `.txt` file is generated for each chromosome.

### 2. Visualization

To reproduce the plots from the paper, open `plots.ipynb`. This notebook imports functions from `utilities.ipynb`, so ensure that the utilities are executed first and all dependencies are installed. You can run the cells in any order to explore dynamics across different cell lines and genomic regions. To save plots, ensure that the `Figures` folder exists in the main directory.

#### Examples

Within `plots.ipynb`, one can run the following code

```
cell_line = "HUVEC"
chr_number = 12
chrpos_min = 80000
chrpos_max = 90000
scale_factor = 6
file_name = f'{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}'
spec_fileQ = False
saveQ = False
rt_plotf(cell_line,chr_number,chrpos_min,chrpos_max,scale_factor,file_name,spec_fileQ,saveQ)
```
![image](https://github.com/user-attachments/assets/6cc62ce9-497f-4a83-b190-3122c1bc2f0e)

If `saveQ = True`, the plot is saved in `figures/file_name.pdf`.

## System requirements and performance

This codebase was developed and tested on Python 3.12.3 under both Windows 10 and Windows 11. No installation procedure is required beyond installing standard Python 3 and the key dependencies. The main libraries include:

- NumPy
- pandas
- matplotlib
- seaborn
- SciPy
- pybigtools
- pyideogram

All scripts should remain compatible with standard Python 3 distributions on other operating systems; however, minor modifications (e.g., file paths) may be necessary. Any compatibility or dependency issues can be reported through the GitHub repository’s issue tracker.

Example scripts and small-scale analyses generally complete in seconds to a few minutes on a standard desktop computer. Visualization routines, in particular, are typically very responsive. Full whole-genome fitting for large datasets can be computationally intensive. We recommend using a high-performance computing (HPC) environment to achieve reasonable run times, especially when processing multiple human cell lines at 1 kb resolution.

## License

This project is openly distributed under the MIT License. This license allows unrestricted use, redistribution, and modification, provided that proper attribution to the original creators is maintained.

## Contact information

For further information, contributions, or queries, please contact:

- **Email**: [fp409@cam.ac.uk](mailto:fp409@cam.ac.uk)
- **GitHub**: [fberkemeier](https://github.com/fberkemeier)

We welcome issues and discussions via GitHub to improve the model or address potential enhancements.

## References

[Berkemeier, F., Cook, P. R., & Boemo, M. A. (2024). DNA replication timing reveals genome-wide features of transcription and fragility. bioRxiv](https://www.biorxiv.org/content/10.1101/2024.11.25.625090v2)
