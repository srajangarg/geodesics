#ifndef __A2_Edge_hpp__
#define __A2_Edge_hpp__
using namespace std;

#include "Common.hpp"
#include <list>

// Forward declarations
class Vertex;
class Face;

/** Edge of a mesh. */
class Edge
{

public:
    /** Construct from two endpoints. */
    Edge(Vertex *v0 = NULL, Vertex *v1 = NULL)
    {
        endpoints[0] = v0;
        endpoints[1] = v1;

        if (v1 < v0)
            std::swap(endpoints[0], endpoints[1]);
    }

    /** Get an endpoint of the edge. \a i = 0 returns the first endpoint and \a
     * i = 1 the second. */
    Vertex const *getEndpoint(int i) const
    {
        debugAssertM(i == 0 || i == 1, "Edge: Invalid endpoint index");
        return endpoints[i];
    }

    /** Get an endpoint of the edge. \a i = 0 returns the first endpoint and \a
     * i = 1 the second. */
    Vertex *getEndpoint(int i)
    {
        debugAssertM(i == 0 || i == 1, "Edge: Invalid endpoint index");
        return endpoints[i];
    }

    /**
     * Given one endpoint of the edge, get the other one. This function assumes
     * the supplied vertex is indeed a valid endpoint
     * of the edge, and (in release mode) does not check that this is so.
     */
    Vertex const *getOtherEndpoint(Vertex const *endpoint) const
    {
        debugAssertM(hasEndpoint(endpoint),
                     "Edge: Vertex is not an endpoint of the edge");
        return endpoints[0] == endpoint ? endpoints[1] : endpoints[0];
    }

    /**
     * Given one endpoint of the edge, get the other one. This function assumes
     * the supplied vertex is indeed a valid endpoint
     * of the edge, and (in release mode) does not check that this is so.
     */
    Vertex *getOtherEndpoint(Vertex const *endpoint)
    {
        debugAssertM(hasEndpoint(endpoint),
                     "Edge: Vertex is not an endpoint of the edge");
        return endpoints[0] == endpoint ? endpoints[1] : endpoints[0];
    }

    /** Get the index (0 or 1) of an endpoint given a pointer to it, or a
     * negative value if the neither endpoint matches. */
    int getEndpointIndex(Vertex const *endpoint) const
    {
        return endpoints[0] == endpoint ? 0 : (endpoints[1] == endpoint ? 1 : -1);
    }

    /** Check if the edge has a given vertex as an endpoint. */
    bool hasEndpoint(Vertex const *v) const
    {
        return endpoints[0] == v || endpoints[1] == v;
    }

    /** Check if the edge is adjacent to a given face. */
    bool hasIncidentFace(Face const *face) const
    {
        for (auto fi = faces.begin(); fi != faces.end(); ++fi)
            if (*fi == face)
                return true;

        return false;
    }

    /** Check if two edges share the same endpoints. */
    bool isCoincidentTo(Edge const &other) const
    {
        return (endpoints[0] == other.endpoints[0] && endpoints[1] == other.endpoints[1])
               || (endpoints[0] == other.endpoints[1]
                   && endpoints[1] == other.endpoints[0]);
    }

    /** Check if this edge shares an endpoint with another. */
    bool isConnectedTo(Edge const &other) const
    {
        return (endpoints[0] == other.endpoints[0] || endpoints[1] == other.endpoints[1]
                || endpoints[0] == other.endpoints[1]
                || endpoints[1] == other.endpoints[0]);
    }

    Vertex *getCommonVertex(Edge *other) const
    {
        if (endpoints[0] == other->endpoints[0] or endpoints[0] == other->endpoints[1])
            return endpoints[0];
        if (endpoints[1] == other->endpoints[0] or endpoints[1] == other->endpoints[1])
            return endpoints[1];
        return NULL;
    }

    /**
     * Get the next edge when stepping counter-clockwise (when viewed from the
     * "outside" of the mesh) around a specified
     * endpoint, assuming the neighborhood is manifold. This also assumes that
     * face vertices/edges have consistent
     * counter-clockwise winding order when viewed from the outside. On error,
     * returns null.
     */
    Edge const *nextAroundEndpoint(int i) const
    {
        return const_cast<Edge *>(this)->nextAroundEndpoint(i);
    }

    /**
     * Get the next edge when stepping counter-clockwise (when viewed from the
     * "outside" of the mesh) around a specified
     * endpoint, assuming the neighborhood is manifold. This also assumes that
     * face vertices/edges have consistent
     * counter-clockwise winding order when viewed from the outside. On error,
     * returns null.
     */
    Edge *nextAroundEndpoint(int i);

    /**
     * Get the next edge when stepping counter-clockwise (when viewed from the
     * "outside" of the mesh) around a specified
     * endpoint, assuming the neighborhood is manifold. This also assumes that
     * face vertices/edges have consistent
     * counter-clockwise winding order when viewed from the outside. On error,
     * returns null.
     */
    Edge const *nextAroundEndpoint(Vertex const *endpoint) const
    {
        return const_cast<Edge *>(this)->nextAroundEndpoint(endpoint);
    }

    /**
     * Get the next edge when stepping counter-clockwise (when viewed from the
     * "outside" of the mesh) around a specified
     * endpoint, assuming the neighborhood is manifold. This also assumes that
     * face vertices/edges have consistent
     * counter-clockwise winding order when viewed from the outside. On error,
     * returns null.
     */
    Edge *nextAroundEndpoint(Vertex const *endpoint)
    {
        debugAssertM(hasEndpoint(endpoint),
                     "Edge: Vertex is not an endpoint of the edge");
        return nextAroundEndpoint(getEndpointIndex(endpoint));
    }

    /** Check if this is a boundary edge, i.e. if it is adjacent to at most one
     * face. */
    bool isBoundary() const
    {
        return faces.size() <= 1;
    }

    friend class Mesh;

    /** Set an endpoint of the edge. */
    void setEndpoint(int i, Vertex *vertex)
    {
        debugAssertM(i == 0 || i == 1, "Edge: Invalid endpoint index");
        endpoints[i] = vertex;
    }

    /** Replace all references to a vertex with references to another vertex. */
    void replaceVertex(Vertex *old_vertex, Vertex *new_vertex)
    {
        if (endpoints[0] == old_vertex)
            endpoints[0] = new_vertex;
        if (endpoints[1] == old_vertex)
            endpoints[1] = new_vertex;
    }

    /** Add a reference to a face incident at this vertex. */
    void addFace(Face *face)
    {
        faces.push_back(face);
    }

    /** Remove all references to a face incident on this edge. */
    void removeFace(Face *face)
    {
        for (auto fi = faces.begin(); fi != faces.end();) {
            if (*fi == face)
                fi = faces.erase(fi);
            else
                ++fi;
        }
    }

    /** Remove a face incident on this edge. */
    std::list<Face *>::iterator removeFace(std::list<Face *>::iterator loc)
    {
        return faces.erase(loc);
    }

    /** Replace all references to a face with references to another face. */
    void replaceFace(Face *old_face, Face *new_face)
    {
        for (auto fi = faces.begin(); fi != faces.end(); ++fi)
            if (*fi == old_face)
                *fi = new_face;
    }

    /** Is the edge a self-loop (both endpoints same)? */
    bool isSelfLoop() const
    {
        return endpoints[0] == endpoints[1];
    }

    double length() const;

    Vertex *endpoints[2];
    std::list<Face *> faces;

    friend ostream &operator<<(ostream &os, const Edge &e);
}; // class Edge

#endif
