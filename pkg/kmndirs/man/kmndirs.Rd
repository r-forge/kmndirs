\name{kmndirs}
\alias{kmndirs}
\title{k-Mean-Directions Clustering}
\description{
  Perform the k-Mean-Directions Algorithm of Maitra
  and Ramler (2010, JCGS).
}
\usage{
kmndirs(x, k, nrandom = 1000L, maxiter = 10L)
}
\arguments{
  \item{x}{a numeric matrix of data to be clustered.}
  \item{k}{an integer giving the number of clusters.}
  \item{nrandom}{an integer giving the number of random initial
    clusterings.}
  \item{maxiter}{an integer giving the maximal number of iterations.}
}
\references{
  R. Maitra and Ivan P. Ramler (2010),
  A k-mean-directions Algorithm for Fast Clustering of Data on the
  Sphere.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{19}/2, 377--396.
  \doi{10.1198/jcgs.2009.08155}.
}
