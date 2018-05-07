# CRISPR-SURF

CRISPR-SURF (**S**creening of **U**ncharacterized **R**egion **F**unction) is an exploratory and interactive computational framework for the analysis of CRISPR-Cas, CRISPRi, and CRISPRa tiling screens.

CRISPR-SURF is available as a user-friendly, open-source software and can be used interactively as a web application at [crisprsurf.pinellolab.org](http://crisprsurf.pinellolab.org/) or as a stand-alone command-line tool with Docker [https://github.com/pinellolab/CRISPR-SURF](https://github.com/pinellolab/CRISPR-SURF).

Installation with Docker
------------------------

With Docker, no installation is required - the only dependence is Docker itself. Users will not need to deal with installation and configuration issues. Docker will do all the dirty work for you!

Docker can be downloaded freely here: [https://store.docker.com/search?offering=community&type=edition](https://store.docker.com/search?offering=community&type=edition)

To get a local copy of CRISPR-SURF, simply execute the following command:
* ```docker pull pinellolab/crisprsurf```

CRISPR-SURF Interactive Website
--------------------------

In order to make CRISPR-SURF more user-friendly and accessible, we have created an interactive website: [http://crisprsurf.pinellolab.org](http://crisprsurf.pinellolab.org). The website implements all the features of the command line version and, in addition, provides interactive and exploratory plots to visualize your CRISPR tiling screen data.

The website offers two functions: 1) Running CRISPR-SURF on CRISPR tiling screen data provided by the user and 2) Visualizing CRISPR-SURF analysis on several published data sets, serving as the first database dedicated to CRISPR tiling screen data.

The website can also run on a local machine using the provided Docker image we have created. To run the website on a local machine after the Docker installation, execute the following command from the command line:
* ```docker run -p 9993:9993 pinellolab/crisprsurf SURF_webapp```

After execution of the command, the user will have a local instance of the website accessible at the URL: 
[http://localhost:9993](http://localhost:9993)
