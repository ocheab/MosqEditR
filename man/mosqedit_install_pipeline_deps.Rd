\name{mosqedit_install_pipeline_deps}
\alias{mosqedit_install_pipeline_deps}
\title{Install optional MosqEdit-R pipeline dependencies}
\usage{mosqedit_install_pipeline_deps(ask = FALSE, update = FALSE)}
\arguments{
  \item{ask}{Logical; whether Bioconductor should prompt.}
  \item{update}{Logical; whether Bioconductor may update installed packages.}
}
\value{Invisibly, a list of requested CRAN and Bioconductor packages.}
\description{Installs optional dependencies required by the full manuscript analysis pipeline while keeping the base package installation lightweight.}

