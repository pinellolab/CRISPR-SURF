# CRISPR-SURF

CRISPR-SURF (**S**creening of **U**ncharacterized **R**egion **F**unction) is an exploratory and interactive computational framework for the analysis of CRISPR-Cas, CRISPRi, and CRISPRa tiling screens.

CRISPR-SURF is available as a user-friendly, open-source software and can be used interactively as a web application at [crisprsurf.pinellolab.org](http://crisprsurf.pinellolab.org/) or as a stand-alone command-line tool with Docker [https://github.com/pinellolab/CRISPR-SURF](https://github.com/pinellolab/CRISPR-SURF).

## Installation with Docker

With Docker, no installation is required - the only dependence is Docker itself. Users will not need to deal with installation and configuration issues. Docker will do all the dirty work for you!

Docker can be downloaded freely here: [https://store.docker.com/search?offering=community&type=edition](https://store.docker.com/search?offering=community&type=edition)

To get a local copy of CRISPR-SURF, simply execute the following command:
* ```docker pull pinellolab/crisprsurf```

## CRISPR-SURF Count

The CRISPR-SURF Count script generates a required input file, the sgRNA summary file, for both the CRISPR-SURF interactive website and command-line interface. You will need one of the following:

* **Option (1)** sgRNA Library File with FASTQs
* **Option (2)** sgRNA Library File with counts

**Option (1)**:
----------------
sgRNA Library File Format Example (.CSV):

| Chr           | Start         | Stop          | sgRNA_Sequence       | Strand | sgRNA_Type       |
| ------------- | ------------- | ------------- | -------------------- | ------ | ---------------- |
| chr2          | 60717499      | 60717519      | AGCTCTGGAATGATGGCTTA | -      | observation      |
| chr2          | 60717506      | 60717526      | ATTGTGGAGCTCTGGAATGA | +      | observation      |
| chr2          | 60717514      | 60717534      | GGAGTTGGATTGTGGAGCTC | +      | observation      |
| chr2          | 60717522      | 60717542      | AGAAAATTGGAGTTGGATTG | -      | negative_control |
| chr2          | 60717529      | 60717549      | CTGGAATAGAAAATTGGAGT | +      | positive_control |

Required Column Names:
* **Chr** - Chromosome
* **Start** - sgRNA Start Genomic Coordinate
* **Stop** - sgRNA Start Genomic Coordinate
* **sgRNA_Sequence** - sgRNA sequence not including PAM sequence
* **Strand** - Targeting strand of the sgRNA
* **sgRNA_Type** - Label for sgRNA type (observation, negative_control, positive_control)

Place the sgRNA Library File and FASTQs in the same directory. The control FASTQs represent the sgRNA distribution prior to selection, while the sample FASTQs represent the sgRNA distribution following selection. Assuming the sgRNA Library File is named sgRNA_library_file.csv and the FASTQs (2 replicates) are named rep1_control.fastq, rep2_control.fastq, rep1_sample.fastq, rep2_sample.fastq, the command-line call would look like:

``` 
docker run -v $PWD:/CRISPR-SURF/SURF -w /CRISPR-SURF/SURF pinellolab/crisprsurf SURF_count -f sgRNA_library_file.csv -control_fastqs rep1_control.fastq rep2_control.fastq -sample_fastqs rep1_sample.fastq rep2_sample.fastq
```
**IMPORTANT:** The number of control FASTQs must equal the number of sample FASTQs. If a single control FASTQ (i.e. plasmid count) is used for multiple sample FASTQs, just enumerate the ```-control_fastqs``` option with the same single control FASTQ.

**Option (2)**:
----------------
sgRNA Library File Format Example (.CSV):

| Chr           | Start         | Stop          | sgRNA_Sequence       | Strand | sgRNA_Type       | Replicate1_Control_Count | Replicate2_Control_Count | Replicate1_Sample_Count | Replicate2_Sample_Count |
| ------------- | ------------- | ------------- | -------------------- | ------ | ---------------- | ------------------------ | ------------------------ | ----------------------- | ----------------------- |
| chr2          | 60717499      | 60717519      | AGCTCTGGAATGATGGCTTA | -      | observation      | 322                      | 615                      | 131                     | 403                     |
| chr2          | 60717506      | 60717526      | ATTGTGGAGCTCTGGAATGA | +      | observation      | 365                      | 812                      | 448                     | 227                     |
| chr2          | 60717514      | 60717534      | GGAGTTGGATTGTGGAGCTC | +      | observation      | 86                       | 169                      | 13                      | 129                     |
| chr2          | 60717522      | 60717542      | AGAAAATTGGAGTTGGATTG | -      | negative_control | 1823                     | 381                      | 1923                    | 321                     |
| chr2          | 60717529      | 60717549      | CTGGAATAGAAAATTGGAGT | +      | positive_control | 54                       | 124                      | 355                     | 521                     |

Required Column Names:
* **Chr** - Chromosome
* **Start** - sgRNA Start Genomic Coordinate
* **Stop** - sgRNA Start Genomic Coordinate
* **sgRNA_Sequence** - sgRNA sequence not including PAM sequence
* **Strand** - Targeting strand of the sgRNA
* **sgRNA_Type** - Label for sgRNA type (observation, negative_control, positive_control)
* **Replicate1_Control_Count** - sgRNA Count in Replicate 1 Control FASTQ (pre-selection)
* **Replicate2_Control_Count** - sgRNA Count in Replicate 2 Control FASTQ (pre-selection)
* **Replicate1_Sample_Count** - sgRNA Count in Replicate 1 Sample FASTQ (post-selection)
* **Replicate2_Sample_Count** - sgRNA Count in Replicate 2 Sample FASTQ (post-selection)

**IMPORTANT:** Additional ReplicateN_Control_Count and ReplicateN_Sample_Count columns can be added depending on the number of replicates used in the experiment. The number of ReplicateN_Control_Count columns must equal ReplicateN_Sample_Count columns. If a single control column (i.e. plasmid count) is used for multiple sample counts, just duplicate the single control column with the appropriate column names.

## CRISPR-SURF Interactive Website

In order to make CRISPR-SURF more user-friendly and accessible, we have created an interactive website: [http://crisprsurf.pinellolab.org](http://crisprsurf.pinellolab.org). The website implements all the features of the command line version and, in addition, provides interactive and exploratory plots to visualize your CRISPR tiling screen data.

The website offers two functions: 1) Running CRISPR-SURF on CRISPR tiling screen data provided by the user and 2) Visualizing CRISPR-SURF analysis on several published data sets, serving as the first database dedicated to CRISPR tiling screen data.

The website can also run on a local machine using the provided Docker image we have created. To run the website on a local machine after the Docker installation, execute the following command from the command line:
* ```docker run -p 9993:9993 pinellolab/crisprsurf SURF_webapp```

After execution of the command, the user will have a local instance of the website accessible at the URL: 
[http://localhost:9993](http://localhost:9993)
