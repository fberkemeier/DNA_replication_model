# A Model of DNA Replication

## Description

This project provides a comprehensive toolkit for analyzing DNA replication timing, origin firing rates, and genomic stability across different cell lines and chromosome regions. The toolkit includes functions to load and process data for various genomic regions, including whole-genome, telomeres, centromeres, and specific genomic sites. By fitting origin firing rates to replication timing data, the toolkit enables efficient prediction and comparison of experimental and predicted timing data. The resulting error distribution between predicted and experimental timings is instrumental in identifying regions of interest. Leveraging datasets that encompass different chromosomes and genomic features, the project facilitates detailed visualization and analysis of replication dynamics and genomic stability.

## Mathematical model

$$\mathbb[T_j]=\sum_{k=0}^{\infty}  \frac{e^{-\sum_{|i|\leq k}(k-|i|)f_{j+i}/v}-e^{-\sum_{|i|\leq k}(k+1-|i|)f_{j+i}/v}}{\sum_{|i|\leq k} f_{j+i}}.$$

## Usage

## Examples

## License

This project is openly distributed under the MIT License. This license allows unrestricted use, redistribution, and modification, provided that proper attribution to the original creators is maintained.

## Contact information

For further information, contributions, or queries, please contact:

- **Email**: [fp409@cam.ac.uk](mailto:fp409@cam.ac.uk)
- **GitHub**: [fberkemeier](https://github.com/fberkemeier)

We welcome issues and discussions via GitHub to improve the model or address potential enhancements.
