#ifndef __A2_Mesh_hpp__
#define __A2_Mesh_hpp__

#include "Common.hpp"
#include "DGP/Graphics/RenderSystem.hpp"
#include "DGP/AxisAlignedBox3.hpp"
#include "DGP/Colors.hpp"
#include "DGP/NamedObject.hpp"
#include "DGP/Noncopyable.hpp"
#include "DGP/Vector3.hpp"
#include "Face.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"
#include <list>
#include <type_traits>
#include <vector>
#include <set>

/** A class for storing meshes with arbitrary topologies. */
class Mesh : public virtual NamedObject, private Noncopyable
{
public:
    /** Constructor. */
    Mesh(std::string const &name = "AnonymousMesh") : NamedObject(name)
    {
    }

    /** Deletes all data in the mesh. */
    void clear()
    {
        vertices.clear();
        edges.clear();
        faces.clear();
        bounds = AxisAlignedBox3();
    }

    /** True if and only if the mesh contains no objects. */
    bool isEmpty() const
    {
        return vertices.empty() && faces.empty() && edges.empty();
    }

    /**
     * Add a vertex to the mesh, with a given location.
     *
     * @return A pointer to the newly created vertex on success, null on
     * failure.
     */
    Vertex *addVertex(Vector3 const &point)
    {
        vertices.push_back(Vertex(point));
        bounds.merge(point);
        return &vertices.back();
    }

    /**
     * Add a vertex to the mesh, with precomputed normal and color.
     *
     * @return A pointer to the newly created vertex on success, null on
     * failure.
     */
    Vertex *addVertex(Vector3 const &point, Vector3 const &normal,
                      ColorRGBA const &color = ColorRGBA(1, 1, 1, 1))
    {
        vertices.push_back(Vertex(point, normal, color));
        bounds.merge(point);
        return &vertices.back();
    }

    /**
     * Add a face to the mesh, specified by the sequence of vertices obtained by
     * dereferencing [vbegin, vend).
     * VertexIterator must dereference to a pointer to a Vertex. Unless the
     * mesh is already in an inconsistent state,
     * failure to add the face will not affect the mesh.
     *
     * @return A pointer to the newly created face, or null on error.
     */
    template <typename VertexIterator>
    Face *addFace(VertexIterator vbegin, VertexIterator vend)
    {
        // Check for errors and compute normal
        size_t num_verts = 0;
        Vector3 v[3];
        for (VertexIterator vi = vbegin; vi != vend; ++vi, ++num_verts) {
            debugAssertM(*vi,
                         getNameStr() + ": Null vertex pointer specified for new face");
            if (num_verts < 3)
                v[num_verts] = (*vi)->getPosition();
        }

        if (num_verts < 3) {
            DGP_WARNING << getName() << ": Skipping face -- too few vertices ("
                        << num_verts << ')';
            return NULL;
        }

        // Create the (initially empty) face
        faces.push_back(Face());
        Face *face = &(*faces.rbegin());

        // Add the loop of vertices to the face
        VertexIterator next = vbegin;
        for (VertexIterator vi = next++; vi != vend; ++vi, ++next) {
            if (next == vend)
                next = vbegin;

            face->addVertex(*vi);
            (*vi)->addFace(face, false); // we'll update the normals later

            Edge *edge = (*vi)->getEdgeTo(*next);
            if (!edge) {
                edges.push_back(Edge(*vi, *next));
                edge = &(*edges.rbegin());

                (*vi)->addEdge(edge);
                (*next)->addEdge(edge);
            }

            edge->addFace(face);
            face->addEdge(edge);
        }

        // Update the face and vertex normals;
        face->updateNormal();
        for (auto fvi = face->vertices.begin(); fvi != face->vertices.end(); ++fvi)
            (*fvi)->addFaceNormal(face->getNormal()); // weight by face area?

        return face;
    }

    /**
     * Remove a face of the mesh. This does NOT remove any vertices or edges.
     * Iterators to the face vector remain valid unless the
     * iterator pointed to the removed face.
     *
     * This is a relatively slow operation since the face needs to be looked up
     * in the face vector
     * (linear in number of faces). For speed, use removeFace(auto).
     *
     * @return True if the face was found and removed, else false.
     */
    bool removeFace(Face *face)
    {
        for (auto fi = faces.begin(); fi != faces.end(); ++fi)
            if (&(*fi) == face)
                return removeFace(fi);

        return false;
    }

    /**
     * Remove a face of the mesh. This does NOT remove any vertices or edges.
     * Iterators to the face vector remain valid unless the
     * iterator pointed to the removed face.
     *
     * Use this version in preference to removeFace(Face const *) where
     * possible.
     *
     * @return True if the face was found and removed, else false.
     */
    bool removeFace(std::list<Face>::iterator face)
    {
        Face *fp = &(*face);

        for (auto fvi = face->vertices.begin(); fvi != face->vertices.end(); ++fvi)
            (*fvi)->removeFace(fp);

        for (auto fei = face->edges.begin(); fei != face->edges.end(); ++fei)
            (*fei)->removeFace(fp);

        faces.erase(face);

        return true;
    }

    /** Draw the mesh on a render_system. */
    void draw(Graphics::RenderSystem &render_system, bool draw_edges = false,
              bool use_vertex_data = false, bool send_colors = false) const;

    /** Update the bounding box of the mesh. */
    void updateBounds()
    {
        bounds = AxisAlignedBox3();
        for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
            bounds.merge(vi->getPosition());
    }

    /** Get the bounding box of the mesh. */
    AxisAlignedBox3 const &getAABB() const
    {
        return bounds;
    }

    /** Load the mesh from a disk file. */
    bool load(std::string const &path);

    /** Save the mesh to a disk file. */
    bool save(std::string const &path) const;

private:
    /**
     * Utility function to draw a face. Must be enclosed in the appropriate
     * RenderSystem::beginPrimitive()/RenderSystem::endPrimitive() block.
     */
    void drawFace(Face const &face, Graphics::RenderSystem &render_system,
                  bool use_vertex_data, bool send_colors) const
    {
        if (!use_vertex_data) {
            render_system.setNormal(face.getNormal());
            if (send_colors)
                render_system.setColor(face.getColor());
        }

        for (auto vi = face.vertices.begin(); vi != face.vertices.end(); ++vi) {
            Vertex const &vertex = **vi;

            if (use_vertex_data) {
                render_system.setNormal(vertex.getNormal());
                if (send_colors)
                    render_system.setColor(vertex.getColor());
            }

            render_system.sendVertex(vertex.getPosition());
        }
    }

    /** Load the mesh from an OFF file. */
    bool loadOFF(std::string const &path);

    /** Save the mesh to an OFF file. */
    bool saveOFF(std::string const &path) const;

public:
    std::list<Face> faces;      ///< Set of mesh faces.
    std::list<Vertex> vertices; ///< Set of mesh vertices.
    std::list<Edge> edges;      ///< Set of mesh edges.
    AxisAlignedBox3 bounds;     ///< Mesh bounding box.

    mutable std::vector<Vertex *>
        face_vertices; ///< Internal cache of vertex pointers for a face.

}; // class Mesh

#endif
