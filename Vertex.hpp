//============================================================================
//
// DGP: Digital Geometry Processing toolkit
// Copyright (C) 2016, Siddhartha Chaudhuri
//
// This software is covered by a BSD license. Portions derived from other
// works are covered by their respective licenses. For full licensing
// information see the LICENSE.txt file.
//
//============================================================================

#ifndef ___Vertex_hpp__
#define ___Vertex_hpp__

#include "Common.hpp"
#include "DGP/Colors.hpp"
#include "DGP/Vector3.hpp"
#include <list>
#define M_PI 3.14159265358979323846

// Forward declarations
class Edge;
class Face;

/** Vertex of a mesh. */
class Vertex
{
public:
    /** Default constructor. */
    Vertex()
        : position(Vector3::zero()), normal(Vector3::zero()),
          color(ColorRGBA(1, 1, 1, 1)), has_precomputed_normal(false),
          normal_normalization_factor(0)
    {
    }

    /** Sets the vertex to have a given location. */
    explicit Vertex(Vector3 const &p)
        : position(p), normal(Vector3::zero()), color(ColorRGBA(1, 1, 1, 1)),
          has_precomputed_normal(false), normal_normalization_factor(0)
    {
    }

    /** Sets the vertex to have a location, normal and color. */
    Vertex(Vector3 const &p, Vector3 const &n, ColorRGBA const &c = ColorRGBA(1, 1, 1, 1))
        : position(p), normal(n), color(c), has_precomputed_normal(true),
          normal_normalization_factor(0)
    {
    }

    /**
     * Get the degree of the vertex, i.e. number of edges incident on it.
     * Equivalent to edges.size()().
     *
     * @see edges.size()();
     */
    int degree() const
    {
        return (int)edges.size();
    }

    /** Get the edge from this vertex to another, if it exists, else return
     * null. */
    Edge const *getEdgeTo(Vertex const *v) const
    {
        return const_cast<Vertex *>(this)->getEdgeTo(v);
    }

    /** Get the edge from this vertex to another, if it exists, else return
     * null. */
    Edge *getEdgeTo(Vertex const *v);

    /** Check if the vertex is adjacent to a given edge. */
    bool hasEdgeTo(Vertex const *v) const
    {
        return getEdgeTo(v) != NULL;
    }

    /** Check if the edge is adjacent to a given face. */
    bool hasIncidentEdge(Edge const *edge) const
    {
        for (auto ei = edges.begin(); ei != edges.end(); ++ei)
            if (*ei == edge)
                return true;

        return false;
    }

    /** Check if the edge is adjacent to a given face. */
    bool hasIncidentFace(Face const *face) const
    {
        for (auto fi = faces.begin(); fi != faces.end(); ++fi)
            if (*fi == face)
                return true;

        return false;
    }

    /** Check if the vertex lies on a mesh boundary. */
    bool isBoundary() const;

    /** Get the position of the vertex. */
    Vector3 const &getPosition() const
    {
        return position;
    }

    /** Set the position of the vertex. */
    void setPosition(Vector3 const &position_)
    {
        update_saddle_or_boundary();
        position = position_;
    }

    /** Get the normal at the vertex. */
    Vector3 const &getNormal() const
    {
        return normal;
    }

    /** Set the normal at the vertex. */
    void setNormal(Vector3 const &normal_)
    {
        normal = normal_;
    }

    /** Check if the vertex has a precomputed normal. */
    bool hasPrecomputedNormal() const
    {
        return has_precomputed_normal;
    }

    /**
     * Update the vertex normal by recomputing it from face data. Useful when
     * the faces have been modified externally. Destroys
     * any prior precomputed normal.
     */
    void updateNormal();

    /** Get the color of the vertex. */
    ColorRGBA const &getColor() const
    {
        return color;
    }

    /** Set the color of the vertex. */
    void setColor(ColorRGBA const &color_)
    {
        color = color_;
    }

private:
    friend class Mesh;

    /** Add a reference to an edge incident at this vertex. */
    void addEdge(Edge *edge)
    {
        edges.push_back(edge);
    }

    /** Remove a reference to an edge incident at this vertex. */
    std::list<Edge *>::iterator removeEdge(std::list<Edge *>::iterator loc)
    {
        return edges.erase(loc);
    }

    /** Remove a reference to an edge incident at this vertex. */
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

    /** Add a reference to a face incident at this vertex. */
    void addFace(Face *face, bool update_normal = true);

    /** Remove all references to a face incident at this vertex. */
    void removeFace(Face *face);

    /** Replace all references to a face with references to another face. */
    void replaceFace(Face *old_face, Face *new_face)
    {
        for (auto fi = faces.begin(); fi != faces.end(); ++fi)
            if (*fi == old_face)
                *fi = new_face;
    }

    /**
     * Add (accumulate) normal information from a new face at this vertex,
     * unless the vertex has a precomputed normal.
     *
     * @param n Unit (or weighted) normal of the new face.
     */
    void addFaceNormal(Vector3 const &n)
    {
        if (!has_precomputed_normal) {
            Vector3 sum_normals = normal_normalization_factor * getNormal() + n;
            normal_normalization_factor = sum_normals.length();
            setNormal(normal_normalization_factor < 1e-20f
                          ? Vector3::zero()
                          : sum_normals / normal_normalization_factor);
        }
    }

    /**
     * Remove (subtract) normal information from a new face at this vertex,
     * unless the vertex has a precomputed normal.
     *
     * @param n Unit (or weighted) normal of the face to be removed.
     */
    void removeFaceNormal(Vector3 const &n)
    {
        if (!has_precomputed_normal) {
            Vector3 sum_normals = normal_normalization_factor * getNormal() - n;
            normal_normalization_factor = sum_normals.length();
            setNormal(normal_normalization_factor < 1e-20f
                          ? Vector3::zero()
                          : sum_normals / normal_normalization_factor);
        }
    }

    void update_saddle_or_boundary();

    Vector3 position;
    Vector3 normal;
    ColorRGBA color;
    std::list<Edge *> edges;
    std::list<Face *> faces;
    bool has_precomputed_normal;
    bool saddle_or_boundary;
    float normal_normalization_factor;

}; // class Vertex

#endif
