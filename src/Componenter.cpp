#include <Rcpp.h>

#include "Array.h"
#include "Componenter.h"

#include "lemon/connectivity.h"

using namespace lemon;

std::vector<int> & Componenter::run ()
{
    Array<double> * kernelArray = kernel->getArray();
    const Neighbourhood &kernelNeighbourhood = kernelArray->getNeighbourhood();
    const Neighbourhood &sourceNeighbourhood = original->getNeighbourhood(kernelArray->getDimensions());
    const size_t neighbourhoodSize = kernelNeighbourhood.size;
    
    const std::vector<int> &dims = original->getDimensions();
    const int nDims = original->getDimensionality();
    const size_t nLabels = original->size();
    labels.resize(nLabels, NA_INTEGER);
    currentLoc.resize(nDims);
    
    connections.clear();
    SmartGraph::NodeMap<size_t> indexMap(connections);
    std::vector<SmartGraph::Node> nodes(original->size(), INVALID);
    
    // Construct the graph
    for (size_t i=0; i<nLabels; i++)
    {
        const double &value = original->at(i);
        if (R_IsNA(value) || value == 0.0)
            continue;
        
        if (nodes[i] == INVALID)
        {
            const SmartGraph::Node node = connections.addNode();
            indexMap[node] = i;
            nodes[i] = node;
        }
        
        original->expandIndex(i, currentLoc);
        
        // We assume the kernel is symmetric (the R code checks this), so we
        // only need to look at half of it
        for (size_t k=(neighbourhoodSize/2)+1; k<neighbourhoodSize; k++)
        {
            const ptrdiff_t loc = i + sourceNeighbourhood.offsets[k];
            
            // Check if we're out of bounds in any dimension
            bool validLoc = true;
            for (int j=0; j<nDims; j++)
            {
                const int currentDimIndex = currentLoc[j] + sourceNeighbourhood.locs(k,j);
                if (currentDimIndex < 0 || currentDimIndex >= dims[j])
                    validLoc = false;
            }
            
            if (!validLoc)
                continue;
            
            const double &neighbourValue = original->at(loc);
            const double &kernelValue = kernelArray->at(k);
            
            // Zero or NA neighbour or kernel value means no connection
            if (R_IsNA(neighbourValue) || neighbourValue == 0.0 || R_IsNA(kernelValue) || kernelValue == 0.0)
                continue;
            
            // Create a node for the neighbour if there isn't already one
            if (nodes[loc] == INVALID)
            {
                const SmartGraph::Node node = connections.addNode();
                indexMap[node] = loc;
                nodes[loc] = node;
            }
            
            // Add an edge
            connections.addEdge(nodes[i], nodes[loc]);
        }
    }
    
    // Find connected components
    SmartGraph::NodeMap<int> componentMemberships(connections);
    connectedComponents(connections, componentMemberships);
    
    for (SmartGraph::NodeIt node(connections); node != INVALID; ++node)
        labels[indexMap[node]] = componentMemberships[node];
    
    return labels;
}
