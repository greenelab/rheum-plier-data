FROM rocker/tidyverse:3.4.3

RUN apt-get update && apt-get install -y \
	git \
	wget \
	curl

RUN R install2.r --error \
	--deps TRUE \
	foreach \
	parallel \
	RCurl

# bioconductor packages
RUN R -e "BiocInstaller::biocLite(c('pd.hg.u133.plus.2', 'pd.hg.u133a', 'pd.hg.u133b', 'pd.hg.u95av2', 'pd.hugene.1.0.st.v1', 'pd.hugene.1.1.st.v1', 'SCAN.UPC', 'affy', 'affyio', 'preprocessCore', 'affxparser', 'illuminaHumanv4.db', 'hgug4112a.db', 'biomaRt', 'org.Hs.eg.db', 'qvalue'))"

# brainarray for affymetrix arrays (v22.0.0)
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133plus2hsentrezgprobe_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133ahsentrezgprobe_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133bhsentrezgprobe_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hgu95av2hsentrezgprobe_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hugene11sthsentrezgprobe_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hugene10sthsentrezgprobe_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133plus2hsentrezgcdf_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133ahsentrezgcdf_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133bhsentrezgcdf_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hgu95av2hsentrezgcdf_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hugene11sthsentrezgcdf_22.0.0.tar.gz')"
RUN R -e "devtools::install_url('http://mbni.org/customcdf/22.0.0/entrezg.download/hugene10sthsentrezgcdf_22.0.0.tar.gz')"


# packages from github
RUN R -e "devtools::install_github('wgmao/PLIER', ref = 'a2d4a2aa343f9ed4b9b945c04326bebd31533d4d', dependencies = TRUE)"
RUN R -e "devtools::install_github('greenelab/TDM', ref = '0bb5b7e4f2478295c9889aa1bf3cf91f1db11006', dependencies = TRUE)"

RUN R -e "install.packages('doParallel')"
