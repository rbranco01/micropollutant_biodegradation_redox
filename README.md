# micropollutant_biodegradation_redox

This is repository contains the analysis of micropollutant removal, microbial activity and microbial communities of an experiment studying micropollutant removal by two soil and two activated sludge communities. The dataset describes the temporal changes of micropollutant and electron acceptor concentrations and soil and activated sludge microbial communities under 5 different redox conditions (i.e. aerobic, nitrate reducing, iron reducing, sulfate reducing and methanogenic).

Chemical input data were obtained with liquid chromatograph-mass spectrometer, gas chromatograph and ion chromatograph. The micropollutant removal efficiencies table and micropollutant and electron acceptor concentrations table were imported into R (in input_data/). The unprocessed/untreated chemical data used to create the R input data has also been included in the repository (in raw_data/).

Molecular biology input data were created by sequencing 16S rRNA gene amplicons using primers 515F and 926R on Illumina Miseq (300 bp PE). Fastq sequence files were analysed using QIIME2 (scripts can be found in the scripts/QIIME2/ directory) and the feature table, taxonomic assignments, phylogeny and metadata were imported into R (in input_data/).

This analysis was published as [**Branco et al. (2023)**](https://doi.org/10.1016/j.scitotenv.2023.165233) in Science of the Total Environment, titled *Influence of redox condition and inoculum on micropollutant biodegradation by soil and activated sludge communities*.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.