% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebi.eqtls.R
\name{ebi.eqtls}
\alias{ebi.eqtls}
\title{ebi.eqtls}
\description{
A template for building documented, tested R6 classes
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-ebi.eqtls-new}{\code{ebi.eqtls$new()}}
\item \href{#method-ebi.eqtls-getStudyId}{\code{ebi.eqtls$getStudyId()}}
\item \href{#method-ebi.eqtls-getCatalog}{\code{ebi.eqtls$getCatalog()}}
\item \href{#method-ebi.eqtls-get.eQTLsByLocationAndStudyID}{\code{ebi.eqtls$get.eQTLsByLocationAndStudyID()}}
\item \href{#method-ebi.eqtls-fetch.eqtls.in.chunks}{\code{ebi.eqtls$fetch.eqtls.in.chunks()}}
\item \href{#method-ebi.eqtls-create.snpLocs.table}{\code{ebi.eqtls$create.snpLocs.table()}}
\item \href{#method-ebi.eqtls-clone}{\code{ebi.eqtls$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ebi.eqtls-new"></a>}}
\if{latex}{\out{\hypertarget{method-ebi.eqtls-new}{}}}
\subsection{Method \code{new()}}{
Creates a new instance of this [R6][R6::R6Class] class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ebi.eqtls$new(study.id = NA)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{study.id}}{character, the ebi-archived eqtl study name}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new instance of ebi.eqtls
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ebi.eqtls-getStudyId"></a>}}
\if{latex}{\out{\hypertarget{method-ebi.eqtls-getStudyId}{}}}
\subsection{Method \code{getStudyId()}}{
accessor for the object's id field
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ebi.eqtls$getStudyId()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the current value of the id member
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ebi.eqtls-getCatalog"></a>}}
\if{latex}{\out{\hypertarget{method-ebi.eqtls-getCatalog}{}}}
\subsection{Method \code{getCatalog()}}{
retrieves the full ebi eqtl catalog
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ebi.eqtls$getCatalog()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ebi.eqtls-get.eQTLsByLocationAndStudyID"></a>}}
\if{latex}{\out{\hypertarget{method-ebi.eqtls-get.eQTLsByLocationAndStudyID}{}}}
\subsection{Method \code{get.eQTLsByLocationAndStudyID()}}{
gets all eQTLs in the region, or just those for the targetGene,
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ebi.eqtls$get.eQTLsByLocationAndStudyID(
  chrom,
  start,
  end,
  studyIDs,
  targetGene = NA,
  simplify
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{chrom}}{character}

\item{\code{start}}{numeric}

\item{\code{end}}{numeric}

\item{\code{studyIDs}}{character from the catalog's unique_id field}

\item{\code{targetGene}}{character default NA restricts eQTLs if supplied}

\item{\code{simplify}}{logical return the full EBI table, or just a few crucial columsnn}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ebi.eqtls-fetch.eqtls.in.chunks"></a>}}
\if{latex}{\out{\hypertarget{method-ebi.eqtls-fetch.eqtls.in.chunks}{}}}
\subsection{Method \code{fetch.eqtls.in.chunks()}}{
gets all eQTLs in the region, chunking requests as specified
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ebi.eqtls$fetch.eqtls.in.chunks(chrom, start, end, study, simplify, chunk.size)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{chrom}}{character}

\item{\code{start}}{numeric}

\item{\code{end}}{numeric}

\item{\code{study}}{character from the catalog's unique_id field}

\item{\code{simplify}}{logical return the full EBI table, or just a few crucial columsnn}

\item{\code{chunk.size}}{integer often 10000, below the limit at which EBI REST fails}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ebi.eqtls-create.snpLocs.table"></a>}}
\if{latex}{\out{\hypertarget{method-ebi.eqtls-create.snpLocs.table}{}}}
\subsection{Method \code{create.snpLocs.table()}}{
create a table of hg19 and hg38 locs for all rsids in region
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ebi.eqtls$create.snpLocs.table(chrom, start, end)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{chrom}}{character}

\item{\code{start}}{numeric}

\item{\code{end}}{numeric}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ebi.eqtls-clone"></a>}}
\if{latex}{\out{\hypertarget{method-ebi.eqtls-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ebi.eqtls$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
