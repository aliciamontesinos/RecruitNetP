#' Graph of nodes and interactions
#'
#' @description
#' Plant-plant interaction network visualization as a graph. On one hand, it is
#' the visualization of the nodes and interactions, which can be visualized in two
#' formats, either as a graph conducted with this function or as an adjacency matrix,
#' with the function *`RN_heatmap`*. And on the other hand, three functions, one for
#' each interaction type network, to visualize the functional topology of general
#' recruitment networks *`visu_funtopol_rec`* and the structural topology of the
#' recruitment enhancement *`visu_topol_fac`* and depression *`visu_topol_depre`*
#' networks respectively.
#'
#' @inheritParams check_interactions
#' @inheritParams check_cover
#'
#' @param int_type Indicates the type of plant-plant interaction that will be
#' analyzed: general recruitment, recruitment enhancement (i.e. facilitation)
#' or recruitment depression (i.e. competition). See detailed options in
#' [int_significance()].
#' @param weight specifies the metric used to represent interaction strength
#' (i.e., the weight) assigned to each pair of species in the matrix.
#' Explanation of its options (more mathematical information in the description
#' of the function **associndex**):
#' - *Fcr*: **frequency of recruitment** in number of recruits by canopy-recruit pair.
#' - *Dcr*: **density of recruitment** as number of recruits per unit area of canopy
#' species.
#' - *Ns*: The index **Normalized Neighbour Suitability index** (proposed by Mingo, 2014),
#' suitable for comparisons of interaction strength between pairs of species within a
#' local community, which should be preferred in general recruitment networks
#' (Alcantara et al. 2025).
#' - *NIntA*: The index **additive symmetry intensity index** proposed by
#' Diaz-Sierra et al. (2017).
#' - *NIntC*: The index **commutative symmetry intensity index** proposed by
#' Diaz-Sierra et al. (2017).
#' -*RII*: The index **Relative Interaction Index** (Armas et al., 2004).

#' @param mode to be used only for recruitment enhancement("fac") and recruitment
#' depression ("comp") networks. Indicates whether the network should be plotted as a
#' unipartite or a bipartite network. In bipartite networks, canopy species are shown
#' in the upper row and recruits in the lower row of the graph. For general
#' recruitment networks, the network should be considered as unipartite, and it will
#' result in an error if this argument is given the option "bi".
#' @param scale_w is an argument to proportionally increase or decrease the thickness
#' of the links. In some networks, high values can result in the overlapping of links
#' that difficult the visualization.
#'
#' @returns a graph representing a network of plant-plant interactions.
#'
#' @export
#'
#' @examples
#' # Unipartite network representation of a general recruitment network. Link width
#' # corresponds to the scaled frequency of recruitment (*Fcr*):
#' visu_net(Amoladeras_int, Amoladeras_cover, int_type="rec", weight="Fcr", mode="uni", scale_w=0.01)
#'
#' # Unipartite representation of a facilitation network. Link width corresponds to the
#' # scaled *Ns* index:
#' visu_net(Amoladeras_int, Amoladeras_cover, int_type="fac", weight="Ns", mode="uni", scale_w=5)
#'
#' # Bipartite representation of a facilitation network. Link width corresponds to the
#' # scaled *Ns* index. Canopy species are shown in the upper row and recruits in the lower
#' # row of the graph:
#' visu_net(Amoladeras_int, Amoladeras_cover, int_type="fac", weight="Ns", mode="bi", scale_w=5)
#'
#' # Unipartite representation of a recruitment depression (*competition*) network. Link
#' # width corresponds to the scaled *RII* index:
#' visu_net(Amoladeras_int, Amoladeras_cover, int_type="comp", weight="RII", mode="uni", scale_w=5)
#'
#' # Bipartite representation of a recruitment depression (*competition*) network. Link
#' # width corresponds to the scaled *RII* index. Canopy species are shown in the upper row
#' # and recruits in the lower row of the graph:
#' visu_net(Amoladeras_int, Amoladeras_cover, int_type="comp", weight="RII", mode="bi", scale_w=5)
#'
#'
visu_net<-function(int_data,cover_data,int_type=c("rec","fac","comp"),
                   weight = c("Pcr","Fcr","Dcr","Dro","Ns", "NintC", "NintA", "RII"),
                   mode= c("uni","bi"), scale_w=1) {

  int_type <- match.arg(int_type)
  weight <- match.arg(weight)
  mode <- match.arg(mode)

  if(int_type=="rec" && mode=="bi"){
    warning("Here recruitment networks are considered as replacement networks sensu Alcantara et al. 2019. JVS, 30:1239-1249, and thus unipartite by definition. Therefore, the combination of int_type=rec,and mode=uni is not provided")
  }

  if(int_type=="rec" && mode=="uni"){

    mat <- suppressWarnings(t(RN_to_matrix(int_data, cover_data, int_type = "rec", weight = weight)))
    edge_list <- as.data.frame(as.table(mat))
    colnames(edge_list) <- c("from", "to", "weight")
    edge_list <- subset(edge_list, weight > 0)
    RN_igraph <- igraph::graph_from_data_frame(edge_list, directed = TRUE)

    if (weight %in% c("Ns", "NintC", "NintA", "RII")) {
      stop("the index specified in the weight argument uses open as a reference, meanwhile recruitment networks  considered it as a node within the network. This creates an inconsistency as open cannot simultaneously function as a node in the network and as a baseline for weighting interactions.")
    }

    scale<-scale_w # a factor used to adjust the magnitude of the weight ( i.e. width of the arrows)
    # We transpose the adjacency matrix so that arrows point from canopy to recruit,
    #this represents which species will replace a the space ocupied by the canopy in the future.
    plot(RN_igraph,
         edge.arrow.size=.3,
         edge.width = igraph::E(RN_igraph)$weight*scale,
         vertex.color="chartreuse",
         vertex.size=8,
         vertex.frame.color="darkolivegreen",
         vertex.label.color="black",
         vertex.label.cex=0.8,
         vertex.label.dist=2,
         vertex.label.font = 3,
         edge.curved=0.2,
         #layout=layout_with_kk(RN_igraph),
         layout=igraph::layout_in_circle(RN_igraph),
         frame = TRUE)
    graphics::title(main="Recruitment Network")
    return(RN_igraph)


  }


  if(int_type=="fac" && mode=="uni"){


    mat <- RN_to_matrix(int_data, cover_data, int_type = "fac", weight = weight)
    edge_list <- as.data.frame(as.table(mat))
    colnames(edge_list) <- c("from", "to", "weight")
    edge_list <- subset(edge_list, weight > 0)

    RN_igraph <- igraph::graph_from_data_frame(edge_list, directed = TRUE)

    scale<-scale_w # a factor used to adjust the magnitude of the weight ( i.e. width of the arrows)
    # We transpose the adjacency matrix so that arrows point from canopy to recruit,
    #this represents which species enhances its recuitment under another species or itself.
    plot(RN_igraph,
         edge.arrow.size=.3,
         edge.width = igraph::E(RN_igraph)$weight*scale,
         vertex.color="chartreuse",
         vertex.size=8,
         vertex.frame.color="darkolivegreen",
         vertex.label.color="black",
         vertex.label.cex=0.8,
         vertex.label.dist=2,
         vertex.label.font = 3,
         edge.curved=0.2,
         #layout=layout_with_kk(RN_igraph),
         layout=igraph::layout_in_circle(RN_igraph),
         frame = TRUE)
    graphics::title(main="Unipartite Recruitment Enhancement Network")

    return(RN_igraph)
  }


  if(int_type=="comp" && mode=="uni"){


    mat <- RN_to_matrix(int_data, cover_data, int_type = "comp", weight = weight)
    edge_list <- as.data.frame(as.table(mat))
    colnames(edge_list) <- c("from", "to", "weight")

    if (weight %in% c("Ns", "NintC", "NintA", "RII")) {
      edge_list$weight<-abs(edge_list$weight)
    }
    edge_list <- subset(edge_list, weight > 0)
    RN_igraph <- igraph::graph_from_data_frame(edge_list, directed = TRUE)


    scale<-scale_w # a factor used to adjust the magnitude of the weight ( i.e. width of the arrows)
    # We transpose the adjacency matrix so that arrows point from canopy to recruit,
    #this represents which species enhances its recuitment under another species or itself.
    plot(RN_igraph,
         edge.arrow.size=.3,
         edge.width = igraph::E(RN_igraph)$weight*scale,
         vertex.color="chartreuse",
         vertex.size=8,
         vertex.frame.color="darkolivegreen",
         vertex.label.color="black",
         vertex.label.cex=0.8,
         vertex.label.dist=2,
         vertex.label.font = 3,
         edge.curved=0.2,
         #layout=layout_with_kk(RN_igraph),
         layout=igraph::layout_in_circle(RN_igraph),
         frame = TRUE)
    graphics::title(main="Unipartite Recruitment Depression Network")

    return(RN_igraph)

  }


  if(int_type=="fac" && mode=="bi"){


    a <- RN_to_matrix(int_data, cover_data, int_type = "fac", weight = weight)
    scale<-scale_w
    a<-a*scale
    sorted_a <- bipartite::sortweb(a, sort.order = "dec")

    # Graficar la red
    bipartite::plotweb(sorted_a,
            srt = 90,
            higher_italic = TRUE,
            lower_italic = TRUE,
            higher_color = "#E69F00",  # Canopy
            lower_color = "#0072B2"   # Recruit
    )


  }


  if(int_type=="comp" && mode=="bi"){


    a <- RN_to_matrix(int_data, cover_data, int_type = "comp", weight = weight)
    scale<-scale_w

    if(weight%in%c("Ns", "NintC", "NintA", "RII")){a<-a*scale*(-1)}else{a<-a*scale}

    sorted_a <- bipartite::sortweb(a, sort.order = "dec")

    # Graficar la red
    bipartite::plotweb(sorted_a,
            srt = 90,
            higher_italic = TRUE,
            lower_italic = TRUE,
            higher_color = "#E69F00",  # Canopy
            lower_color = "#0072B2"   # Recruit
    )

  }
}

