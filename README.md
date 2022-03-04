# Vigilin/HDLBP analysis scripts

This is a collection of analysis scripts that accompany the publication "HDLBP binds ER-targeted mRNAs by multivalent interactions to promote protein synthesis of transmembrane and secreted proteins." (under revision).

The scripts allow reproduction of key figures obtained from RNA-seq, PAR-CLIP, ribo-seq, BioID and pSILAC datasets.

[Raw data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148262) are available from GEO. 

Please contact milekm@gmail.com for further information.

Authors:
Ulrike Zinnall<sup>1</sup><sup>,</sup><sup>3</sup><sup>,</sup><sup>+</sup>, Miha Milek<sup>1</sup><sup>,</sup><sup>2</sup><sup>,</sup><sup>+</sup><sup>,</sup><sup>#</sup>, Carlos H. Vieira-Vieira<sup>1</sup><sup>,</sup><sup>3</sup>, Simon Müller<sup>4</sup>, Guido Mastrobuoni<sup>1</sup>, Orsalia-Georgia Hazapis<sup>1</sup>, Igor Minia<sup>1</sup>, Simone del Guidice<sup>1</sup>, Nadine Bley<sup>4</sup>, Stefan Kempa<sup>1</sup>, Stefan Hüttelmaier<sup>4</sup>, Matthias Selbach<sup>1</sup><sup>,</sup><sup>5</sup>, Markus Landthaler<sup>1</sup><sup>,</sup><sup>6</sup><sup>*</sup>

<sup>1</sup>Max Delbrück Center for Molecular Medicine in the Helmholtz Association, Berlin Institute for Medical Systems Biology, Germany.
<sup>2</sup>National Institute of Chemistry, Slovenia.
<sup>3</sup>Faculty of Life Sciences, Humboldt-Universität zu Berlin, Germany.
<sup>4</sup>Institute of Molecular Medicine, Medical Faculty, Martin-Luther-University, Germany.
<sup>5</sup>Charite-Universitätsmedizin Berlin, Germany.
<sup>6</sup>IRI Life Sciences, Institute of Biology, Humboldt-Universität zu Berlin, Germany


Please start in the figs directory, which enables easy reproduction of key figures from processed data. Working directory should be set to project directory. This should also be enough to perform further exploratory analysis.


Documentation and reproduction of key tables in the data directory is available in the make_master_tables directory.


For the quantification of tRNAs from PAR-CLIP data, please see the scripts in the trna directory.


In the geo_processed_data directory we also provide the processed data that are otherwise also available from GEO. Required for some R scripts in the figs directory.


In the source_data directory the tables used for source data Excel file are for each panel.


In the suppl_tables directory the production of Supplemental Tables is documented.