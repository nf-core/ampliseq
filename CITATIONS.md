# nf-core/ampliseq: Citations

## [nf-core/ampliseq](https://pubmed.ncbi.nlm.nih.gov/33193131/)

> Straub D, Blackwell N, Langarica-Fuentes A, Peltzer A, Nahnsen S, Kleindienst S. Interpretations of Environmental Microbial Community Studies Are Biased by the Selected 16S rRNA (Gene) Amplicon Sequencing Pipeline. Front Microbiol. 2020 Oct 23;11:550420. doi: 10.3389/fmicb.2020.550420. PMID: 33193131; PMCID: PMC7645116.

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

### Core tools

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

  > Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].

- [Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200/479)

  > Marcel, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal 17.1 (2011): pp-10. doi: 10.14806/ej.17.1.200.

- [Barrnap](https://github.com/tseemann/barrnap)

  > Seemann T. barrnap 0.9 : rapid ribosomal RNA prediction.

- [DADA2](https://pubmed.ncbi.nlm.nih.gov/27214047/)
  > Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods. 2016 Jul;13(7):581-3. doi: 10.1038/nmeth.3869. Epub 2016 May 23. PMID: 27214047; PMCID: PMC4927377.

### Taxonomic classification and databases

- Classification by [QIIME2 classifier](https://pubmed.ncbi.nlm.nih.gov/29773078/)

  > Bokulich NA, Kaehler BD, Rideout JR, Dillon M, Bolyen E, Knight R, Huttley GA, Gregory Caporaso J. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin. Microbiome. 2018 May 17;6(1):90. doi: 10.1186/s40168-018-0470-z. PMID: 29773078; PMCID: PMC5956843.

- default: [SILVA](https://pubmed.ncbi.nlm.nih.gov/23193283/)

  > Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013 Jan;41(Database issue):D590-6. doi: 10.1093/nar/gks1219. Epub 2012 Nov 28. PMID: 23193283; PMCID: PMC3531112.

- [Greengenes2](https://doi.org/10.1038/s41587-023-01845-1)

  > McDonald, D., Jiang, Y., Balaban, M. et al. Greengenes2 unifies microbial data in a single reference tree. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01845-1

- [PR2 - Protist Reference Ribosomal Database](https://pubmed.ncbi.nlm.nih.gov/23193267/)

  > Guillou L, Bachar D, Audic S, Bass D, Berney C, Bittner L, Boutte C, Burgaud G, de Vargas C, Decelle J, Del Campo J, Dolan JR, Dunthorn M, Edvardsen B, Holzmann M, Kooistra WH, Lara E, Le Bescot N, Logares R, Mahé F, Massana R, Montresor M, Morard R, Not F, Pawlowski J, Probert I, Sauvadet AL, Siano R, Stoeck T, Vaulot D, Zimmermann P, Christen R. The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. Nucleic Acids Res. 2013 Jan;41(Database issue):D597-604. doi: 10.1093/nar/gks1160. Epub 2012 Nov 27. PMID: 23193267; PMCID: PMC3531120.

- [GTDB - Genome Taxonomy Database](https://pubmed.ncbi.nlm.nih.gov/30148503/)

  > Parks DH, Chuvochina M, Waite DW, Rinke C, Skarshewski A, Chaumeil PA, Hugenholtz P. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. Nat Biotechnol. 2018 Nov;36(10):996-1004. doi: 10.1038/nbt.4229. Epub 2018 Aug 27. PMID: 30148503.

- [SBDI-GTDB](https://scilifelab.figshare.com/articles/dataset/SBDI_Sativa_curated_16S_GTDB_database/14869077)

  > Lundin D, Andersson A. SBDI Sativa curated 16S GTDB database. FigShare. doi: 10.17044/scilifelab.14869077.v1

- [RDP - Ribosomal Database Project](https://pubmed.ncbi.nlm.nih.gov/24288368/)

  > Cole JR, Wang Q, Fish JA, Chai B, McGarrell DM, Sun Y, Brown CT, Porras-Alfaro A, Kuske CR, Tiedje JM. Ribosomal Database Project: data and tools for high throughput rRNA analysis. Nucleic Acids Res. 2014 Jan;42(Database issue):D633-42. doi: 10.1093/nar/gkt1244. Epub 2013 Nov 27. PMID: 24288368; PMCID: PMC3965039.

- [UNITE - eukaryotic nuclear ribosomal ITS region](https://pubmed.ncbi.nlm.nih.gov/15869663/)

  > Kõljalg U, Larsson KH, Abarenkov K, Nilsson RH, Alexander IJ, Eberhardt U, Erland S, Høiland K, Kjøller R, Larsson E, Pennanen T, Sen R, Taylor AF, Tedersoo L, Vrålstad T, Ursing BM. UNITE: a database providing web-based methods for the molecular identification of ectomycorrhizal fungi. New Phytol. 2005 Jun;166(3):1063-8. doi: 10.1111/j.1469-8137.2005.01376.x. PMID: 15869663.

- [MIDORI2 - a collection of reference databases](https://doi.org/10.1002/edn3.303/)

  > Leray, M., Knowlton, N., & Machida, R. J. (2022). MIDORI2: A collection of quality controlled, preformatted, and regularly updated reference databases for taxonomic assignment of eukaryotic mitochondrial sequences. Environmental DNA, 4, 894– 907. doi: https://doi.org/10.1002/edn3.303.

- [COIDB - CO1 Taxonomy Database](https://doi.org/10.17044/scilifelab.20514192.v2)

  > Sundh J, Manoharan L, Iwaszkiewicz-Eggebrecht E, Miraldo A, Andersson A, Ronquist F. COI reference sequences from BOLD DB. doi: https://doi.org/10.17044/scilifelab.20514192.v2.

- [PhytoRef plastid 16S rRNA database for photosynthetic eukaryotes](https://pubmed.ncbi.nlm.nih.gov/25740460/)

  > Decelle J, Romac S, Stern RF, Bendif el M, Zingone A, Audic S, Guiry MD, Guillou L, Tessier D, Le Gall F, Gourvil P, Dos Santos AL, Probert I, Vaulot D, de Vargas C, Christen R. PhytoREF: a reference database of the plastidial 16S rRNA gene of photosynthetic eukaryotes with curated taxonomy. Mol Ecol Resour. 2015 Nov;15(6):1435-45. doi: 10.1111/1755-0998.12401. Epub 2015 Apr 6. PMID: 25740460.

- [Zehr lab nifH database](http://doi.org/10.5281/zenodo.7996213)

  > M. A. Moynihan & C. Furbo Reeder 2023. nifHdada2 GitHub repository, v2.0.5. Zenodo. doi: http://doi.org/10.5281/zenodo.7996213

### Phylogenetic placement

- [nf-core/phyloplace](https://nf-co.re/phyloplace)

  > Daniel Lundin. (2023). nf-core/phyloplace: First release (1.0.0). Zenodo. https://doi.org/10.5281/zenodo.7643948

- [HMMER](https://pubmed.ncbi.nlm.nih.gov/22039361/)

  > Eddy, Sean R. “Accelerated Profile HMM Searches.” PLoS Comput Biol 7, no. 10 (October 20, 2011): e1002195. https://doi.org/10.1371/journal.pcbi.1002195.

- [MAFFT](https://pubmed.ncbi.nlm.nih.gov/12136088/)

  > Katoh, Kazutaka, Kazuharu Misawa, Kei‐ichi Kuma, and Takashi Miyata. “MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based on Fast Fourier Transform.” Nucleic Acids Research 30, no. 14 (July 15, 2002): 3059–66. https://doi.org/10.1093/nar/gkf436.

- [EPA-NG](https://pubmed.ncbi.nlm.nih.gov/30165689/)

  > Barbera, Pierre, Alexey M Kozlov, Lucas Czech, Benoit Morel, Diego Darriba, Tomáš Flouri, and Alexandros Stamatakis. “EPA-Ng: Massively Parallel Evolutionary Placement of Genetic Sequences.” Systematic Biology 68, no. 2 (March 1, 2019): 365–69. https://doi.org/10.1093/sysbio/syy054.

- [Gappa](https://pubmed.ncbi.nlm.nih.gov/32016344/)

  > Czech, Lucas, Pierre Barbera, and Alexandros Stamatakis. “Genesis and Gappa: Processing, Analyzing and Visualizing Phylogenetic (Placement) Data.” Bioinformatics 36, no. 10 (May 1, 2020): 3263–65. https://doi.org/10.1093/bioinformatics/btaa070.

### Multi region analysis (also include Greengenes 13_8 or SILVA 128)

- [q2-sidle](https://doi.org/10.1101/2021.03.23.436606)

  > Debelius, J.W.; Robeson, M.; Lhugerth, L.W.; Boulund, F.; Ye, W.; Engstrand, L. "A comparison of approaches to scaffolding multiple regions along the 16S rRNA gene for improved resolution." Preprint in BioRxiv. doi: 10.1101/2021.03.23.436606

- [SMURF](https://doi.org/10.1186/s40168-017-0396-x)

  > Fuks, G.; Elgart, M.; Amir, A.; Zeisel, A.; Turnbaugh, P.J., Soen, Y.; and Shental, N. (2018). "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." Microbiome. 6: 17. doi: 10.1186/s40168-017-0396-x

- [RESCRIPt](https://doi.org/10.1371/journal.pcbi.1009581)

  > Robeson MS 2nd, O'Rourke DR, Kaehler BD, Ziemski M, Dillon MR, Foster JT, Bokulich NA. RESCRIPt: Reproducible sequence taxonomy reference database management. PLoS Comput Biol. 2021 Nov 8;17(11):e1009581. doi: 10.1371/journal.pcbi.1009581. PMID: 34748542; PMCID: PMC8601625.

- [SEPP](https://doi.org/10.1128/msystems.00021-18)

  > Janssen S, McDonald D, Gonzalez A, Navas-Molina JA, Jiang L, Xu ZZ, Winker K, Kado DM, Orwoll E, Manary M, Mirarab S, Knight R. Phylogenetic Placement of Exact Amplicon Sequences Improves Associations with Clinical Information. mSystems. 2018 Apr 17;3(3):e00021-18. doi: 10.1128/mSystems.00021-18. PMID: 29719869; PMCID: PMC5904434.

### Downstream analysis

- [QIIME2](https://pubmed.ncbi.nlm.nih.gov/31341288/)

  > Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS 2nd, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, Caporaso JG. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol. 2019 Aug;37(8):852-857. doi: 10.1038/s41587-019-0209-9. Erratum in: Nat Biotechnol. 2019 Sep;37(9):1091. PMID: 31341288; PMCID: PMC7015180.

- [MAFFT](https://pubmed.ncbi.nlm.nih.gov/23329690/)

  > Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013 Apr;30(4):772-80. doi: 10.1093/molbev/mst010. Epub 2013 Jan 16. PMID: 23329690; PMCID: PMC3603318.

- [ANCOM](https://pubmed.ncbi.nlm.nih.gov/26028277/)

  > Mandal S, Van Treuren W, White RA, Eggesbø M, Knight R, Peddada SD. Analysis of composition of microbiomes: a novel method for studying microbial composition. Microb Ecol Health Dis. 2015 May 29;26:27663. doi: 10.3402/mehd.v26.27663. PMID: 26028277; PMCID: PMC4450248.

- [Adonis](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x) and [VEGAN](https://CRAN.R-project.org/package=vegan)

  > Marti J Anderson. A new method for non-parametric multivariate analysis of variance. Austral ecology, 26(1):32–46, 2001.

  > Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O’Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs, and Helene Wagner. vegan: Community Ecology Package. 2018. R package version 2.5-3.

- [Phyloseq](https://doi.org/10.1371/journal.pone.0061217)

  > McMurdie PJ, Holmes S (2013). “phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data.” PLoS ONE, 8(4), e61217.

### Non-default tools

- [ITSx](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12073)

  > Bengtsson-Palme, J., Ryberg, M., Hartmann, M., Branco, S., Wang, Z., Godhe, A., De Wit, P., Sánchez-García, M., Ebersberger, I., de Sousa, F., Amend, A., Jumpponen, A., Unterseher, M., Kristiansson, E., Abarenkov, K., Bertrand, Y.J.K., Sanli, K., Eriksson, K.M., Vik, U., Veldre, V. and Nilsson, R.H.. Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods Ecol Evol 2013, 4: 914-919. doi: 10.1111/2041-210X.12073.

- [PICRUSt2](https://pubmed.ncbi.nlm.nih.gov/32483366/)

  > Douglas GM, Maffei VJ, Zaneveld JR, Yurgel SN, Brown JR, Taylor CM, Huttenhower C, Langille MGI. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol. 2020 Jun;38(6):685-688. doi: 10.1038/s41587-020-0548-6. PMID: 32483366; PMCID: PMC7365738.

- PICRUSt2 is by default using [EPA-ng](https://pubmed.ncbi.nlm.nih.gov/30165689/)

  > Barbera P, Kozlov AM, Czech L, Morel B, Darriba D, Flouri T, Stamatakis A. EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences. Syst Biol. 2019 Mar 1;68(2):365-369. doi: 10.1093/sysbio/syy054. PMID: 30165689; PMCID: PMC6368480.

- PICRUSt2 is by default using [MinPath](https://pubmed.ncbi.nlm.nih.gov/19680427/)

  > Ye Y, Doak TG. A parsimony approach to biological pathway reconstruction/inference for genomes and metagenomes. PLoS Comput Biol. 2009 Aug;5(8):e1000465. doi: 10.1371/journal.pcbi.1000465. Epub 2009 Aug 14. PMID: 19680427; PMCID: PMC2714467.

- [VSEARCH](https://peerj.com/articles/2584/)

  > Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: a versatile open source tool for metagenomics. PeerJ. 2016 4:e2584. doi: 10.7717/peerj.2584

- VSEARCH option usearch_global implements the [USEARCH](https://doi.org/10.1093/bioinformatics/btq461) algorithm

  > Edgar RC. Search and clustering orders of magnitude faster than BLAST. Bioinformatics. 2010 26(19) 2460-2461

- VSEARCH option sintax implements the [SINTAX](https://doi.org/10.1101/074161) algorithm

  > Edgar RC. (2016) SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences, BioRxiv, 074161. Preprint.

- [Kraken2](https://pubmed.ncbi.nlm.nih.gov/31779668/)

  > Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome biology, 20(1), 257. https://doi.org/10.1186/s13059-019-1891-0

### Summarizing software

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

## Data

- [Full-size test data](https://doi.org/10.3389/fmicb.2020.550420)
  > Straub D, Blackwell N, Langarica-Fuentes A, Peltzer A, Nahnsen S, Kleindienst S. Interpretations of Environmental Microbial Community Studies Are Biased by the Selected 16S rRNA (Gene) Amplicon Sequencing Pipeline. Front Microbiol. 2020 Oct 23;11:550420. doi: 10.3389/fmicb.2020.550420. PMID: 33193131; PMCID: PMC7645116.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
