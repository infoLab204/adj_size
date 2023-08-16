########################################################
# create graph with p-values and significance level alpha
########################################################
cgraph <- function(p_values, n, alpha)
{
	# find number of edges
	pcorr_edges <- (p_values < alpha)
	# length(p_values[pcorr_edges])   # number of edges in graph

	#creating the corresponding network graph
	# see #7.18 and #7.19 in p.125 of the reference
	mat_pcorr <- matrix(0, n, n)
	mat_pcorr[lower.tri(mat_pcorr)] <- as.numeric(pcorr_edges)
	g_pcorr <- graph.adjacency(mat_pcorr, "undirected")

	# find disconnected components. 
	# see p. 57 of “Statistical analysis of netowork data with R”
	#comps <- decompose.graph(g_pcorr)
	#t <- table(sapply(comps, vcount))

	# see #4.14 in p.52-54 of the reference. 
	# t <- table(sapply(cliques(g_pcorr), length))   
	# t <- table(sapply(maximal.cliques(g_pcorr), length))   
	d <- graph.density(g_pcorr) 
	t <- transitivity(g_pcorr)
	return(list(d, t))
	#return(list(g_pcorr, t))
}
#####################################################
# (eg) gres <- cgraph(res[[1]], res[[2]], 0.05)
# ref: Kolaczyk, "Statistical analysis of network data with R"



