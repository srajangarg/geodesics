#include "Vertex.hpp"
#include "Edge.hpp"
#include "Face.hpp"

Edge *Vertex::getEdgeTo(Vertex const *v)
{
    if (v == this)
        return NULL;

    for (auto ei = edges.begin(); ei != edges.end(); ++ei) {
        Edge *e = *ei;
        if (e->hasEndpoint(v))
            return e;
    }

    return NULL;
}

bool Vertex::isBoundary() const
{
    if (edges.empty())
        return true;

    for (auto ei = edges.begin(); ei != edges.end(); ++ei)
        if ((*ei)->isBoundary())
            return true;

    return false;
}

void Vertex::addFace(Face *face, bool update_normal)
{
    faces.push_back(face);
    if (update_normal && !has_precomputed_normal)
        addFaceNormal(face->getNormal());
}

void Vertex::removeFace(Face *face)
{
    for (auto fi = faces.begin(); fi != faces.end();) {
        if (*fi == face) {
            fi = faces.erase(fi);
            if (!has_precomputed_normal)
                removeFaceNormal(face->getNormal());

            // Keep going, just in case the face somehow got added twice
        } else
            ++fi;
    }
}

void Vertex::updateNormal()
{
    if (!faces.empty()) {
        Vector3 sum_normals = Vector3::zero();
        for (auto fi = faces.begin(); fi != faces.end(); ++fi)
            sum_normals += (*fi)->getNormal(); // weight by face area?

        normal_normalization_factor = sum_normals.length();
        setNormal(normal_normalization_factor < 1e-20f
                      ? Vector3::zero()
                      : sum_normals / normal_normalization_factor);
    } else {
        setNormal(Vector3::zero());
        normal_normalization_factor = 0;
    }

    has_precomputed_normal = false;
}

void Vertex::update_saddle_or_boundary()
{
    saddle_or_boundary = false;
    if (isBoundary()) {
        saddle_or_boundary = true;
    } else {
        double angle = 0.0;
        for (auto fi = faces.begin(); fi != faces.end(); ++fi)
            angle += (*fi)->getAngle(this);
        if (angle > 2 * M_PI)
            saddle_or_boundary = true;
    }
}