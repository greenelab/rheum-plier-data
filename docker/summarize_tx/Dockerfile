FROM argrosso/rbioconductor:3.4.3

RUN apt-get update && apt-get install -y software-properties-common
RUN apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 0xF1656F24C74CD1D8
RUN add-apt-repository 'deb [arch=amd64,i386,ppc64el] http://ftp.utexas.edu/mariadb/repo/10.3/ubuntu xenial main'

RUN apt-get update && apt-get install -y libmariadbclient-dev

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/RMySQL_0.10.15.tar.gz', type = 'source', repos = NULL)"
RUN R -e "BiocInstaller::biocLite('ensembldb')"
