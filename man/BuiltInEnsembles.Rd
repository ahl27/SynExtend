\name{BuiltInEnsembles}
\alias{BuiltInEnsembles}
\docType{data}
\title{
Pretrained EvoWeaver Ensemble Models
}
\description{
EvoWeaver has best performance with an ensemble method combining individual
evidence streams.
This data file provides pretrained models for ease of use.
Two groups of models are provided:
1. Models trained on the KEGG MODULES dataset
2. Models trained on the CORUM dataset

These models are used internally if the user does not provide their own model,
and aren't explicitly designed to be accessed by the user.

See the examples for how to train your own ensemble model.
}
\usage{data("BuiltInEnsembles")}
\format{
  The data contain a named list of objects of class \code{\link[stats]{glm}}.
  This list currently has two entries: \code{"KEGG"} and \code{"CORUM"}
}

\examples{
## Training own ensemble method to avoid using built-ins
## defaults to built-ins when an ensemble isn't provided
set.seed(333L)
exData <- get(data("ExampleStreptomycesData"))
ew <- EvoWeaver(exData$Genes[seq_len(8L)], MySpeciesTree=exData$Tree, NoWarn=TRUE)
datavals <- predict(ew, NoPrediction=TRUE, Verbose=interactive())
datavals <- datavals[datavals[,1] != datavals[,2],]

# Picking random numbers for demonstration purposes
actual_values <- sample(rep(c(1,0), length.out=nrow(datavals)))
datavals[,'y'] <- actual_values
myModel <- glm(y~., datavals[,-c(1,2)], family='binomial')

predictionPW <- EvoWeaver(exData$Genes[9:10], MySpeciesTree=exData$Tree, NoWarn=TRUE)
predict(predictionPW,
          PretrainedModel=myModel, Verbose=interactive())[2,,drop=FALSE]
}
\keyword{datasets}
