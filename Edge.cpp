#include "Edge.hpp"
#include "Face.hpp"
#include "Vertex.hpp"
#include "DGP/Vector4.hpp"

Edge *Edge::nextAroundEndpoint(int i)
{
    debugAssertM(i == 0 || i == 1, "Edge: Invalid endpoint index");

    if (faces.size() > 2) // non-manifold
        return NULL;

    // Find which incident face has this endpoint as the origin of the edge when
    // stepping round the face. The required edge
    // is then the predecessor of this edge around the face.
    for (auto fi = faces.begin(); fi != faces.end(); ++fi) {
        Face *face = *fi;
        Edge *prev = face->getPredecessor(this);
        if (prev->hasEndpoint(endpoints[i])) // found it!
            return prev;
    }

    return NULL;
}

double Edge::length() const
{
    return (endpoints[0]->getPosition() - endpoints[1]->getPosition()).length();
}
