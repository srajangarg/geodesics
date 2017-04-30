#include "Mesh.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"
#include "Face.hpp"
#include "DGP/FilePath.hpp"
#include "DGP/Polygon3.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_map>

void Mesh::draw(Graphics::RenderSystem &render_system, bool draw_edges,
                bool use_vertex_data, bool send_colors) const
{
    // Three separate passes over the faces is probably faster than using
    // Primitive::POLYGON for each face

    if (draw_edges) {
        render_system.pushShapeFlags();
        render_system.setPolygonOffset(true, 1);
    }

    // First try to render as much stuff using triangles as possible
    render_system.beginPrimitive(Graphics::RenderSystem::Primitive::TRIANGLES);
    for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        if (fi->isTriangle())
            drawFace(*fi, render_system, use_vertex_data, send_colors);
    render_system.endPrimitive();

    // Now render all quads
    render_system.beginPrimitive(Graphics::RenderSystem::Primitive::QUADS);
    for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        if (fi->isQuad())
            drawFace(*fi, render_system, use_vertex_data, send_colors);
    render_system.endPrimitive();

    // Finish off with all larger polygons
    for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        if (fi->edges.size() > 4) {
            render_system.beginPrimitive(Graphics::RenderSystem::Primitive::POLYGON);
            drawFace(*fi, render_system, use_vertex_data, send_colors);
            render_system.endPrimitive();
        }

    if (draw_edges)
        render_system.popShapeFlags();

    if (draw_edges) {
        render_system.pushShader();
        render_system.pushColorFlags();

        render_system.setShader(NULL);
        render_system.setColor(ColorRGBA(0.2, 0.3, 0.7, 1)); // set default edge color

        render_system.beginPrimitive(Graphics::RenderSystem::Primitive::LINES);
        for (auto ei = edges.begin(); ei != edges.end(); ++ei) {
            render_system.sendVertex(ei->getEndpoint(0)->getPosition());
            render_system.sendVertex(ei->getEndpoint(1)->getPosition());
        }
        render_system.endPrimitive();

        render_system.popColorFlags();
        render_system.popShader();
    }
}

bool Mesh::loadOFF(std::string const &path)
{
    std::ifstream in(path.c_str());
    if (!in) {
        DGP_ERROR << "Could not open '" << path << "' for reading";
        return false;
    }

    clear();

    std::string magic;
    if (!(in >> magic) || magic != "OFF") {
        DGP_ERROR << "Header string OFF not found at beginning of file '" << path << '\'';
        return false;
    }

    long nv, nf, ne;
    if (!(in >> nv >> nf >> ne)) {
        DGP_ERROR << "Could not read element counts from OFF file '" << path << '\'';
        return false;
    }

    if (nv < 0 || nf < 0 || ne < 0) {
        DGP_ERROR << "Negative element count in OFF file '" << path << '\'';
        return false;
    }

    vertices.reserve(nv);
    edges.reserve(max(ne, (nv + nf) * 2 + 10));
    faces.reserve(nf);

    std::vector<Vertex *> indexed_vertices;
    Vector3 p;
    for (long i = 0; i < nv; ++i) {
        if (!(in >> p[0] >> p[1] >> p[2])) {
            DGP_ERROR << "Could not read vertex " << indexed_vertices.size() << " from '"
                      << path << '\'';
            return false;
        }

        Vertex *v = addVertex(p);
        if (!v)
            return false;
        v->index = i;
        indexed_vertices.push_back(v);
    }

    std::vector<Vertex *> face_vertices;
    long num_face_vertices, vertex_index;
    bool triangular = true;

    for (long i = 0; i < nf; ++i) {
        if (!(in >> num_face_vertices) || num_face_vertices < 0) {
            DGP_ERROR << "Could not read valid vertex count of face " << faces.size()
                      << " from '" << path << '\'';
            return false;
        }

        if (num_face_vertices != 3)
          triangular = false;

        face_vertices.resize(num_face_vertices);
        for (size_t j = 0; j < face_vertices.size(); ++j) {
            if (!(in >> vertex_index)) {
                DGP_ERROR << "Could not read vertex " << j << " of face " << faces.size()
                          << " from '" << path << '\'';
                return false;
            }

            if (vertex_index < 0 || vertex_index >= (long)vertices.size()) {
                DGP_ERROR << "Out-of-bounds index " << vertex_index << " of vertex " << j
                          << " of face " << faces.size() << " from '" << path << '\'';
                return false;
            }

            face_vertices[j] = indexed_vertices[(size_t)vertex_index];
        }

        addFace(
            face_vertices.begin(),
            face_vertices.end()); // ok if this fails, just skip the face with a warning
    }

    for (auto &v : vertices)
        v.update_saddle_or_boundary();

    setName(FilePath::objectName(path));

    return triangular;
}

bool Mesh::saveOFF(std::string const &path) const
{
    std::ofstream out(path.c_str(), std::ios::binary);
    if (!out) {
        DGP_ERROR << "Could not open '" << path << "' for writing";
        return false;
    }

    out << "OFF\n";
    out << vertices.size() << ' ' << faces.size() << " 0\n";

    typedef std::unordered_map<Vertex const *, long> VertexIndexMap;
    VertexIndexMap vertex_indices;
    long index = 0;
    for (auto vi = vertices.begin(); vi != vertices.end(); ++vi, ++index) {
        Vector3 const &p = vi->getPosition();
        out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';

        vertex_indices[&(*vi)] = index;
    }

    for (auto fi = faces.begin(); fi != faces.end(); ++fi) {
        out << fi->vertices.size();

        for (auto vi = fi->vertices.begin(); vi != fi->vertices.end(); ++vi) {
            VertexIndexMap::const_iterator existing = vertex_indices.find(*vi);
            if (existing == vertex_indices.end()) {
                DGP_ERROR << "Face references vertex absent from mesh '" << path << '\'';
                return false;
            }

            out << ' ' << existing->second;
        }

        out << '\n';
    }

    return true;
}

bool Mesh::load(std::string const &path)
{
    std::string path_lc = toLower(path);
    bool status = false;
    if (endsWith(path_lc, ".off"))
        status = loadOFF(path);
    else {
        DGP_ERROR << "Unsupported mesh format: " << path;
    }

    if (not status)
    {
      vector<vector<Face>::iterator> to_erase;
      vector<Face> new_faces;

      for (auto fi = faces.begin(); fi != faces.end(); ++fi) {

        if (fi->vertices.size() == 3) 
        {
            new_faces.push_back(*fi);
            continue;
        }
        
        to_erase.push_back(fi);
        vector<Vertex*> vv;
        Polygon3 p;
        
        for (auto fvi = fi->vertices.begin(); fvi != fi->vertices.end(); fvi++)
        {
          p.addVertex((*fvi)->getPosition());
          vv.push_back(*fvi);

        }

        std::vector<long> tri_indx;
        p.triangulate(tri_indx);

        for (int i = 0; i < tri_indx.size(); i += 3)
        {
            Face f;
            f.addVertex(vv[tri_indx[i]]);
            f.addVertex(vv[tri_indx[i+1]]);
            f.addVertex(vv[tri_indx[i+2]]);
            new_faces.push_back(f);
        }
      }

      cout << "Saving new triangular mesh!" << endl;
      faces.swap(new_faces);
      saveOFF(path.substr(0, path.find_last_of(".")) + "_tri.off");
    }

    return status;
}

bool Mesh::save(std::string const &path) const
{
    std::string path_lc = toLower(path);
    if (endsWith(path_lc, ".off"))
        return saveOFF(path);

    DGP_ERROR << "Unsupported mesh format: " << path;
    return false;
}
