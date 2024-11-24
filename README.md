# A Model of DNA Replication

## Description

This project offers a comprehensive toolkit for analysing DNA replication timing, origin firing rates, and genomic stability across cell lines and chromosomal regions. Complementing the work in Berkemeier et al. (2024), it includes functions for loading and processing data across whole-genome regions, telomeres, centromeres, and specific loci of interest. By fitting origin firing rates to replication timing data, the toolkit efficiently predicts and compares experimental and modelled timing profiles. The resulting error distributions between predicted and experimental data help pinpoint regions of interest. With datasets spanning diverse chromosomes and genomic features, this toolkit enables detailed visualisation and analysis of replication dynamics and genomic stability.

⚠️ **Note:** This repository is a work in progress and represents an ongoing project. The code and documentation are subject to updates and refinements, and while we strive for accuracy, they may not yet reflect the final, polished version. Your understanding and feedback are appreciated.



## Mathematical model

Consider a DNA molecule with $n$ discrete genomic loci, where each locus $j$ can potentially act as an origin that fires at rate $f_j$ (indepedently of other origins) to initiate a fork that progresses bidirectionally with speed $v$, typically measured in kilobases per minute (kb/min). Let $T_j$ be the time a site $j$ takes to fire or be passively replicated by a fork. Then, for a sufficiently large genome (such as human), the expected time of replication at $j$ is given by

$$
\mathbb{E}[T_j]=\sum_{k=0}^{\infty}  \frac{e^{-\sum_{|i|\leq k}(k-|i|)f_{j+i}/v}-e^{-\sum_{|i|\leq k}(k+1-|i|)f_{j+i}/v}}{\sum_{|i|\leq k} f_{j+i}}.
$$

While this expression holds true for an infinitely large genome, in practical terms this series can be limited to $0\leq k\leq  R<n/2$, for some large enough $R$. This parameter represents the radius of replication influence: the distance within which neighbouring origins $\{j-R,...,j-1,j+1,...,j+R\}$ are assumed to affect the timing of a focal origin $j$. In other words, while every firing origin does theoretically affect replication timing at any other location, this effect decays rapidly with distance from the origin of interest $j$. Numerically, the finite version of the equation should mimic the average replication timing obtained from computational simulations.

## Usage

Clone the repository to a local directory. All data generated in Berkemeier et al. (2024) should be download from here. Add this to a subfolder 'data' within the main directory folder. Bedgraph files of timing error and origin firing rates can be uploaded to the Genome Browser for comparison with other data.

### 1. Data generation

If you want to upload your own data, you can use this tool to upload timing data and convert it to origin firing rates via our model. You can run this locally or on a HPC platform. We recommend installing [`pybigtools`](https://pypi.org/project/pybigtools/) to deal with bigWig files. All utilities and necessary packages are listed in `utilities.ipynb`. DNA replication is simulated using the [Beacon Calculus](https://github.com/MBoemo/bcs).

#### Examples

As an example, we look at importing a bigWig file for HUVEC, from the [ENCODE database](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeUwRepliSeq), also available at 'data/bigwig_files/HUVEC.bw' (in the data folder).

### 2. Visualization

To reproduce the plots in the paper, open `plots.ipynb`. This notebook imports the functions defined in `utilities.ipynb`, so make sure these are called first, and all dependencies are installed. Run the notebook in any order, and explore different dynamics for different cell lines and genomic regions. In order to save plots, make sure the folder `figures` is present in the main directory.

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


## License

This project is openly distributed under the MIT License. This license allows unrestricted use, redistribution, and modification, provided that proper attribution to the original creators is maintained.

## Contact information

For further information, contributions, or queries, please contact:

- **Email**: [fp409@cam.ac.uk](mailto:fp409@cam.ac.uk)
- **GitHub**: [fberkemeier](https://github.com/fberkemeier)

We welcome issues and discussions via GitHub to improve the model or address potential enhancements.
