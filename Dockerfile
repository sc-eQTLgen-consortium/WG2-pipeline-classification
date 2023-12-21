
################## BASE IMAGE ######################

FROM ubuntu:22.04

################## METADATA ######################

LABEL base_image="ubuntu:22.04"
LABEL version="1.0.0"
LABEL software="WG2 Pipeline"
LABEL about.summary="WG2 sceQTLGen Consortium Classification Pipeline"
LABEL about.documentation="https://github.com/sc-eQTLgen-consortium/WG2-pipeline-classification"
LABEL about.tags="Genomics"

# Build syntax: docker build ./ -t wg2-pipeline-classification:2023.11.24.0 --progress=plain > build.log 2>&1
# Total build takes ? minutes and has a size of 1.49 GB.
# Use dive wg2-pipeline-classification:2023.11.24.0 to investigate memory usage.

################## MAINTAINER ######################

MAINTAINER Martijn Vochteloo <m.vochteloo@umcg.nl>

################## INSTALLATION ######################

ADD . /tmp/repo
WORKDIR /tmp/repo

ENV PATH=/opt:/usr/games:$PATH
ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

# Needed to prevent asking for geographic location when installing things.
RUN export TZ=Europe/Amsterdam \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

RUN apt-get update -y \
    # libc-bin libc6 libsystemd0 libudev1
    && apt-get upgrade -y \
    # binutils binutils-common binutils-x86-64-linux-gnu build-essential bzip2 cpp
    # cpp-11 dpkg-dev g++ g++-11 gcc gcc-11 gcc-11-base libasan6 libatomic1
    # libbinutils libc-dev-bin libc6-dev libcc1-0 libcrypt-dev libctf-nobfd0
    # libctf0 libdpkg-perl libgcc-11-dev libgdbm-compat4 libgdbm6 libgomp1
    # libisl23 libitm1 liblsan0 libmpc3 libmpfr6 libnsl-dev libperl5.34
    # libquadmath0 libstdc++-11-dev libtirpc-dev libtsan0 libubsan1 linux-libc-dev
    # lto-disabled-list make patch perl perl-modules-5.34 rpcsvc-proto xz-utils
    && apt-get install -y --no-install-recommends build-essential \
    # ca-certificates openssl
    && apt-get install -y --no-install-recommends ca-certificates \
    # libpsl5 wget
    && apt-get install -y --no-install-recommends wget

#############################
############# R #############
#############################

RUN apt-get install -y --no-install-recommends cmake \
    # install two helper packages we need: dirmngr and software-properties-common
    && apt-get install -y --no-install-recommends dirmngr \
    # dbus distro-info-data gir1.2-glib-2.0 gir1.2-packagekitglib-1.0 gpg gpgconf
    # iso-codes libapparmor1 libappstream4 libargon2-1 libassuan0 libcap2-bin
    # libcryptsetup12 libcurl3-gnutls libdbus-1-3 libdevmapper1.02.1 libdw1
    # libelf1 libgirepository-1.0-1 libglib2.0-0 libglib2.0-bin libglib2.0-data
    # libgstreamer1.0-0 libip4tc2 libjson-c5 libkmod2 libmpdec3
    # libpackagekit-glib2-18 libpam-systemd libpolkit-agent-1-0
    # libpolkit-gobject-1-0 libpython3-stdlib libpython3.10-minimal
    # libpython3.10-stdlib libreadline8 libsqlite3-0 libstemmer0d libunwind8
    # libxmlb2 libyaml-0-2 lsb-release media-types packagekit pkexec policykit-1
    # polkitd python-apt-common python3 python3-apt python3-blinker
    # python3-cffi-backend python3-cryptography python3-dbus python3-distro
    # python3-gi python3-httplib2 python3-importlib-metadata python3-jeepney
    # python3-jwt python3-keyring python3-launchpadlib python3-lazr.restfulclient
    # python3-lazr.uri python3-minimal python3-more-itertools python3-oauthlib
    # python3-pkg-resources python3-pyparsing python3-secretstorage python3-six
    # python3-software-properties python3-wadllib python3-zipp python3.10
    # python3.10-minimal readline-common software-properties-common systemd
    # systemd-sysv
    && apt-get install -y --no-install-recommends software-properties-common \
    # add the signing key (by Michael Rutter) for these repos
    # To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    # Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
    && wget -qO- "https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc" \
    && tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    # add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    # fontconfig fontconfig-config fonts-dejavu-core libblas3 libbsd0 libcairo2
    # libcurl4 libdatrie1 libdeflate0 libfontconfig1 libfreetype6 libfribidi0
    # libgfortran5 libgomp1 libgraphite2-3 libharfbuzz0b libice6 libjbig0
    # libjpeg-turbo8 libjpeg8 liblapack3 libmd0 libpango-1.0-0 libpangocairo-1.0-0
    # libpangoft2-1.0-0 libpaper-utils libpaper1 libpixman-1-0 libpng16-16
    # libquadmath0 libsm6 libtcl8.6 libthai-data libthai0 libtiff5 libtk8.6
    # libwebp7 libx11-6 libx11-data libxau6 libxcb-render0 libxcb-shm0 libxcb1
    # libxdmcp6 libxext6 libxft2 libxrender1 libxss1 libxt6 r-base r-base-core
    # r-cran-boot r-cran-class r-cran-cluster r-cran-codetools r-cran-foreign
    # r-cran-kernsmooth r-cran-lattice r-cran-mass r-cran-matrix r-cran-mgcv
    # r-cran-nlme r-cran-nnet r-cran-rpart r-cran-spatial r-cran-survival
    # r-recommended tzdata ucf unzip x11-common xdg-utils zip
    && apt-get install -y --no-install-recommends r-base \
    # binutils binutils-common binutils-x86-64-linux-gnu build-essential bzip2 cpp
    # cpp-11 dpkg-dev g++ g++-11 gcc gcc-11 gcc-11-base gfortran gfortran-11
    # icu-devtools libasan6 libatomic1 libbinutils libblas-dev libbz2-dev
	# libc-dev-bin libc6-dev libcc1-0 libcrypt-dev libctf-nobfd0 libctf0
	# libdpkg-perl libgcc-11-dev libgdbm-compat4 libgdbm6 libgfortran-11-dev
	# libicu-dev libisl23 libitm1 libjpeg-dev libjpeg-turbo8-dev libjpeg8-dev
	# liblapack-dev liblsan0 liblzma-dev libmpc3 libmpfr6 libncurses-dev
	# libncurses5-dev libnsl-dev libpcre16-3 libpcre2-16-0 libpcre2-32-0
	# libpcre2-dev libpcre2-posix3 libpcre3-dev libpcre32-3 libpcrecpp0v5
	# libperl5.34 libpng-dev libreadline-dev libstdc++-11-dev libtirpc-dev
	# libtsan0 libubsan1 libxmuu1 linux-libc-dev lto-disabled-list make patch perl
	# perl-modules-5.34 pkg-config r-base-dev rpcsvc-proto xauth xz-utils
	# zlib1g-dev
    && apt-get install -y --no-install-recommends r-base-dev \
    # libbrotli-dev libexpat1-dev libfontconfig-dev libfontconfig1-dev
    # libfreetype-dev libfreetype6-dev uuid-dev
    && apt-get install -y --no-install-recommends libfontconfig1-dev \
    # libcurl4-openssl-dev
    && apt-get install -y --no-install-recommends libcurl4-openssl-dev \
    # libxml2-dev
    && apt-get install -y --no-install-recommends libxml2-dev \
    # libssl-dev
    && apt-get install -y --no-install-recommends libssl-dev \
    # gir1.2-harfbuzz-0.0 libblkid-dev libffi-dev libglib2.0-dev
    # libglib2.0-dev-bin libgraphite2-dev libharfbuzz-dev libharfbuzz-gobject0
    # libharfbuzz-icu0 libmount-dev libselinux1-dev libsepol-dev python3-distutils
    # python3-lib2to3
    && apt-get install -y --no-install-recommends libharfbuzz-dev \
    # libfribidi-dev
    && apt-get install -y --no-install-recommends libfribidi-dev \
    # libdeflate-dev libjbig-dev libtiff-dev libtiff5-dev libtiffxx5
    && apt-get install -y --no-install-recommends libtiff5-dev \
    # Required for hdf5r
    && apt-get install -y --no-install-recommends libhdf5-dev \
    # Required for scCustomize
    && apt-get install -y --no-install-recommends libcairo2-dev

# remotes_2.4.2.1
RUN R --slave -e 'install.packages("remotes")' \
    # parallelly_1.36.0 listenv_0.9.0 globals_0.16.2 digest_0.6.33 future_1.33.0
    && R --slave -e 'remotes::install_version("future.apply", version = "1.11.0", upgrade=FALSE)' \
    # getopt_1.20.4
    && R --slave -e 'remotes::install_version("optparse", version = "1.7.3", upgrade=FALSE)' \
    # rlang_1.1.2 lifecycle_1.0.4 glue_1.6.2 cli_3.6.1 vctrs_0.6.4 utf8_1.2.4 fansi_1.0.5
    # colorspace_2.1-0 pkgconfig_2.0.3 pillar_1.9.0 magrittr_2.0.3 viridisLite_0.4.2 RColorBrewer_1.1-3
    # R6_2.5.1 munsell_0.5.0 labeling_0.4.3 farver_2.1.1 withr_2.5.2 tibble_3.2.1 scales_1.3.0 isoband_0.2.7
    # gtable_0.3.4
    && R --slave -e 'remotes::install_version("ggplot2", version = "3.4.3", upgrade=FALSE)' \
    # Matrix_1.6-4
    && R --slave -e 'remotes::install_version("Matrix", version = "1.6-4", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("tidyverse", version = "2.0.0", upgrade=FALSE)' \
    # Rcpp_1.0.11 RcppEigen_0.3.3.9.4 progressr_0.14.0 sp_2.1-2
    && R --slave -e 'remotes::install_version("SeuratObject", version = "4.1.4", upgrade=FALSE)' \
    # sitmo_2.0.2 BH_1.81.0-1 spatstat.utils_3.0-4 tensor_1.5 abind_1.4-5 polyclip_1.10-6 deldir_2.0-2 
    # spatstat.geom_3.2-7 spatstat.data_3.0-3 fastmap_1.1.1 rappdirs_0.3.3 fs_1.6.3 ellipsis_0.3.2 
    # sass_0.4.7 mime_0.12 memoise_2.0.1 jsonlite_1.8.7 jquerylib_0.1.4 htmltools_0.5.7 cachem_1.0.8 
    # base64enc_0.1-3 promises_1.2.1 later_1.3.1 dotCall64_1.1-1 stringi_1.8.2 stringr_1.5.1 plyr_1.8.9 
    # bitops_1.0-7 caTools_1.18.2 gtools_3.9.5 rprojroot_2.0.4 lazyeval_0.2.2 tidyselect_1.2.0 
    # generics_0.1.3 cpp11_0.4.6 purrr_1.0.2 dplyr_1.1.4 tinytex_0.49 fontawesome_0.5.2 bslib_0.6.1 
    # xfun_0.41 highr_0.10 evaluate_0.23 yaml_2.3.7 rmarkdown_2.25 knitr_1.45 sys_3.4.2 askpass_1.2.0 
    # openssl_2.1.1 curl_5.1.0 commonmark_1.9.0 crayon_1.5.2 sourcetools_0.1.7-1 xtable_1.8-4 httpuv_1.6.12 
    # png_0.1-8 here_1.0.1 RcppTOML_0.2.2 dqrng_0.3.2 RcppProgress_0.4.2 irlba_2.3.5.1 RcppAnnoy_0.0.21 
    # FNN_1.1.3.2 goftest_1.2-3 spatstat.sparse_3.0-3 spatstat.random_3.2-2 spam_2.10-0 RcppArmadillo_0.12.6.6.0 
    # matrixStats_1.1.0 gridExtra_2.3 reshape2_1.4.4 gplots_3.1.3 data.table_1.14.8 crosstalk_1.2.1 tidyr_1.3.0 
    # htmlwidgets_1.6.3 httr_1.4.7 shiny_1.8.0 zoo_1.8-12 igraph_1.5.1 reticulate_1.34.0 uwot_0.1.16 
    # spatstat.explore_3.2-5 sctransform_0.4.1 scattermore_1.2 Rtsne_0.16 ROCR_1.0-11 RANN_2.6.1 plotly_4.10.3 
    # pbapply_1.7-2 patchwork_1.1.3 miniUI_0.1.1.1 lmtest_0.9-40 leiden_0.4.3.1 ica_1.0-3 ggridges_0.5.4 
    # ggrepel_0.9.4 fitdistrplus_1.1-11 cowplot_1.1.1
    && R --slave -e 'remotes::install_version("Seurat", version = "4.4.0", upgrade=FALSE)' \
    # Using commit of Nov 4, 2023: SeuratDisk_0.0.0.9021
    # scPred_af5492e77 SQUAREM_2021.1 numDeriv_2016.8-1.1 shape_1.4.6 lava_1.7.3 diagram_1.6.5
    # timechange_0.2.0 prodlim_2023.08.28 tzdb_0.4.0 timeDate_4022.108 lubridate_1.9.3 ipred_0.9-14
    # hardhat_1.3.0 gower_1.0.1 clock_0.7.0 iterators_1.0.14 proxy_0.4-27 openssl_2.1.1 curl_5.1.0 httr_1.4.7
    # recipes_1.0.8 pROC_1.18.5 ModelMetrics_1.2.2.2 foreach_1.5.2 e1071_1.7-13 plotly_4.10.3 caret_6.0-94
    # data.tree_1.1.0 Seurat_5.0.1 harmony_f054b030b recipes_1.0.8 pROC_1.18.5 ModelMetrics_1.2.2.2
    # foreach_1.5.2 e1071_1.7-13 vipor_0.4.5 beeswarm_0.4.0 RhpcBLASctl_0.23-42 kernlab_0.9-32 caret_6.0-94
    # MLmetrics_1.1.1 ggbeeswarm_0.7.2
    && R --slave -e 'remotes::install_github("mojaveazure/seurat-disk@877d4e18ab38c686f5db54f8cd290274ccdbe295", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("scCustomize", version = "1.1.3", upgrade=FALSE)' \
    # None
    && R --slave -e 'remotes::install_version("viridis", version = "0.6.4", upgrade=FALSE)' \
    # Using commit of Nov 21, 2021: HierscPred_0.1.0
    && R --slave -e 'remotes::install_github("joseah/HierscPred@2ce27feeabdd22ce210f1a9623db91fe7ec6b0e7", upgrade=FALSE)'

####################################
################ CLEAN #############
####################################

RUN apt-get clean \
    && apt-get autoremove -y
