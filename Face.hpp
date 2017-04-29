#ifndef __A2_Face_hpp__
#define __A2_Face_hpp__

#include "Common.hpp"
using namespace std;
#include "DGP/Colors.hpp"
#include "DGP/Vector3.hpp"
#include <list>

// Forward declarations
class Vertex;
class Edge;

/** Face of a mesh. */
class Face
{
public:
    /** Construct with the given normal. */
    Face(Vector3 const &normal_ = Vector3::zero()) : normal(normal_)
    {
    }

    /** Check if the face has a given vertex. */
    bool hasVertex(Vertex const *vertex) const
    {
        for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
            if (*vi == vertex)
                return true;

        return false;
    }

    /** Check if the face has a given edge. */
    bool hasEdge(Edge const *edge) const
    {
        for (auto ei = edges.begin(); ei != edges.end(); ++ei)
            if (*ei == edge)
                return true;

        return false;
    }

    /** Get the predecessor of a vertex around the face. Assumes the iterator
     * points to a valid vertex of the face. */
    Vertex const *getPredecessor(std::list<Vertex *>::iterator vertex) const
    {
        return const_cast<Face *>(this)->getPredecessor(vertex);
    }

    /** Get the predecessor of a vertex around the face. Assumes the iterator
     * points to a valid vertex of the face. */
    Vertex *getPredecessor(std::list<Vertex *>::iterator vertex)
    {
        if (vertices.empty())
            return NULL;
        else if (vertex == vertices.begin())
            return vertices.back();
        else {
            --vertex;
            return *vertex;
        }
    }

    /** Get the predecessor of a vertex around the face. Returns null if the
     * vertex does not belong to the face. */
    Vertex const *getPredecessor(Vertex const *vertex) const
    {
        return const_cast<Face *>(this)->getPredecessor(vertex);
    }

    /** Get the predecessor of a vertex around the face. Returns null if the
     * vertex does not belong to the face. */
    Vertex *getPredecessor(Vertex const *vertex)
    {
        for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
            if (*vi == vertex)
                return getPredecessor(vi);

        return NULL;
    }

    /** Get the successor of a vertex around the face. Assumes the iterator
     * points to a valid vertex of the face. */
    Vertex const *getSuccessor(std::list<Vertex *>::iterator vertex) const
    {
        return const_cast<Face *>(this)->getSuccessor(vertex);
    }

    /** Get the successor of a vertex around the face. Assumes the iterator
     * points to a valid vertex of the face. */
    Vertex *getSuccessor(std::list<Vertex *>::iterator vertex)
    {
        if (vertices.empty())
            return NULL;
        else {
            ++vertex;
            return vertex == vertices.end() ? vertices.front() : *vertex;
        }
    }

    /** Get the successor of a vertex around the face. Returns null if the
     * vertex does not belong to the face. */
    Vertex const *getSuccessor(Vertex const *vertex) const
    {
        return const_cast<Face *>(this)->getSuccessor(vertex);
    }

    /** Get the successor of a vertex around the face. Returns null if the
     * vertex does not belong to the face. */
    Vertex *getSuccessor(Vertex const *vertex)
    {
        for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
            if (*vi == vertex)
                return getSuccessor(vi);

        return NULL;
    }

    /** Get the predecessor of an edge around the face. Assumes the iterator
     * points to a valid edge of the face. */
    Edge const *getPredecessor(std::list<Edge *>::iterator edge) const
    {
        return const_cast<Face *>(this)->getPredecessor(edge);
    }

    /** Get the predecessor of an edge around the face. Assumes the iterator
     * points to a valid edge of the face. */
    Edge *getPredecessor(std::list<Edge *>::iterator edge)
    {
        if (edges.empty())
            return NULL;
        else if (edge == edges.begin())
            return edges.back();
        else {
            --edge;
            return *edge;
        }
    }

    /** Get the predecessor of an edge around the face. Returns null if the edge
     * does not belong to the face. */
    Edge const *getPredecessor(Edge const *edge) const
    {
        return const_cast<Face *>(this)->getPredecessor(edge);
    }

    /** Get the predecessor of an edge around the face. Returns null if the edge
     * does not belong to the face. */
    Edge *getPredecessor(Edge const *edge)
    {
        for (auto ei = edges.begin(); ei != edges.end(); ++ei)
            if (*ei == edge)
                return getPredecessor(ei);

        return NULL;
    }

    /** Get the successor of an edge around the face. Assumes the iterator
     * points to a valid edge of the face. */
    Edge const *getSuccessor(std::list<Edge *>::iterator edge) const
    {
        return const_cast<Face *>(this)->getSuccessor(edge);
    }

    /** Get the successor of an edge around the face. Assumes the iterator
     * points to a valid edge of the face. */
    Edge *getSuccessor(std::list<Edge *>::iterator edge)
    {
        if (edges.empty())
            return NULL;
        else {
            ++edge;
            return edge == edges.end() ? edges.front() : *edge;
        }
    }

    /** Get the successor of an edge around the face. Returns null if the edge
     * does not belong to the face. */
    Edge const *getSuccessor(Edge const *edge) const
    {
        return const_cast<Face *>(this)->getSuccessor(edge);
    }

    /** Get the successor of an edge around the face. Returns null if the edge
     * does not belong to the face. */
    Edge *getSuccessor(Edge const *edge)
    {
        for (auto ei = edges.begin(); ei != edges.end(); ++ei)
            if (*ei == edge)
                return getSuccessor(ei);

        return NULL;
    }

    /** Reverse the order in which vertices and edges wind around the face. The
     * face normal is <b>not</b> modified. */
    void reverseWinding()
    {
        vertices.reverse();
        edges.reverse();
    }

    /** Get the face normal. */
    Vector3 const &getNormal() const
    {
        return normal;
    }

    /** Set the face normal. */
    void setNormal(Vector3 const &normal_)
    {
        normal = normal_;
    }

    /** Update the face normal by recomputing it from vertex data. */
    void updateNormal();

    /** Get the color of the face. */
    ColorRGBA const &getColor() const
    {
        return color;
    }

    /** Set the color of the face. */
    void setColor(ColorRGBA const &color_)
    {
        color = color_;
    }

    bool isTriangle() const
    {
        return vertices.size() == 3;
    }

    bool isQuad() const
    {
        return vertices.size() == 4;
    }

    // get angle on this face at this vertex
    double getAngle(Vertex *v);

    /**
     * Test if the face contains a point (which is assumed to lie on the plane
     * of the face -- for efficiency the function does
     * <b>not</b> explicitly verify that this holds).
     */
    bool contains(Vector3 const &p) const;

public:
    friend ostream &operator<<(ostream &os, const Face &f);

    friend class Mesh;
    friend class Edge;

    /** Add a reference to a vertex of this face. */
    void addVertex(Vertex *vertex)
    {
        vertices.push_back(vertex);
    }

    /** Remove a reference to a vertex. */
    std::list<Vertex *>::iterator removeVertex(std::list<Vertex *>::iterator loc)
    {
        return vertices.erase(loc);
    }

    /** Remove all references to a vertex. */
    void removeVertex(Vertex *vertex)
    {
        for (auto vi = vertices.begin(); vi != vertices.end();) {
            if (*vi == vertex)
                vi = vertices.erase(vi);
            else
                ++vi;
        }
    }

    /** Replace all references to an vertex with references to another vertex.
     */
    void replaceVertex(Vertex *old_vertex, Vertex *new_vertex)
    {
        for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
            if (*vi == old_vertex)
                *vi = new_vertex;
    }

    /** Add a reference to an edge of this face. */
    void addEdge(Edge *edge)
    {
        edges.push_back(edge);
    }

    /** Remove a reference to an edge. */
    std::list<Edge *>::iterator removeEdge(std::list<Edge *>::iterator loc)
    {
        return edges.erase(loc);
    }

    /** Remove all references to an edge. */
    void removeEdge(Edge *edge)
    {
        for (auto ei = edges.begin(); ei != edges.end();) {
            if (*ei == edge)
                ei = edges.erase(ei);
            else
                ++ei;
        }
    }

    /** Replace all references to an edge with references to another edge. */
    void replaceEdge(Edge *old_edge, Edge *new_edge)
    {
        for (auto ei = edges.begin(); ei != edges.end(); ++ei)
            if (*ei == old_edge)
                *ei = new_edge;
    }

    Vector3 normal;
    ColorRGBA color;
    std::list<Vertex *> vertices;
    std::list<Edge *> edges;

}; // class Face

#endif
