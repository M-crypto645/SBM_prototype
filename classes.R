#' @title The ICL (implemented from the paper of Daudin, 2008)
#'
#' @description an abstract class for which the ICL can be calculated
#' @slot ICL the ICL will be stored here after calculation
#' @slot cl the approximation of the hidden variable Z_tilde
#' @slot x adjacency matrix of the sample graph X
#' @slot Q number of clusters the model assumed
#' @param sel_frac fraction of best solutions selected for crossing  (default to 0.75)
#' @return a \code{Genetic-class} object
#' @examples
#' Genetic()
#' Genetic(pop_size = 500)
setClassUnion("adj_matrix_type", c("matrix", "dgCMatrix"))
setClass("ICL_object", 
         slots=list(ICL="numeric", cl="vector",
                    x="adj_matrix_type", Q="numeric"))
setGeneric("ICL", function(object) sim_ICL(object))
setMethod("ICL",
          "ICL_object",
          function(object) {
            if(length(object@ICL)==0){
              res = sim_ICL(object@x, object@cl, object@Q)
              object@ICL = res
            } else {
              res = object@ICL
            }
            res
          }
)

#' @title The (misnamed) Variational EM Algorithm (implemented from the paper of Daudin, 2008)
#'
#' @description calculates a clustering of the nodes in Q classes based on the
#' Variational EM Algorithm (extends ICL_object)
#' @slot ICL an ICL value can be stored here
#' @slot cl the approximation of the hidden variable Z_tilde
#' @slot x adjacency matrix of the sample graph X
#' @slot info some additional information on pi,mu, the e- and m-step
#' @slot Q number of clusters the model assumed
setClass("Var_Bayes_output", 
         slots=list(ICL="numeric", cl="vector", info="list", 
                    x="adj_matrix_type", Q="numeric"), contains="ICL_object")