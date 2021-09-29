# ReGSP

## A visualized application for homology-based gene searching and plotting using multiple reference sequences

A multiple-Reference–based Gene Search and Plot (ReGSP) is a simple and convenient web tool that accepts multiple reference sequences for homology-based gene search. The tool incorporates cPlot, a novel dot plot tool, for illustrating nucleotide sequence similarity between the query and the reference sequences.

#### files in this repository
- Website (production level file built using Angular 10)
- api (Flask api for web communication with back-end)
- regsp (phyton file for gene searching and MUMmer plot)
- cPlot (python file nucleotide sequence similarity analysis using *k-mer*)

The `ReGSP` web tool for analysis can be accessed at this [site](https://ds.mju.ac.kr/regsp)