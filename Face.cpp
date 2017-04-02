#include "Face.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"

void Face::updateNormal()
{
    // Assume the face is planar.
    auto vi2 = vertices.begin();
    auto vi0 = vi2++;
    auto vi1 = vi2++;

    if (vertices.size() > 3) {
        // vi1 might be a concave corner -- we need to add up the cross products
        // at all vertices
        Vector3 sum_cross = Vector3::zero();
        for (; vi0 != vertices.end(); ++vi0, ++vi1, ++vi2) {
            if (vi1 == vertices.end())
                vi1 = vertices.begin();
            if (vi2 == vertices.end())
                vi2 = vertices.begin();

            Vector3 e1 = (*vi0)->getPosition() - (*vi1)->getPosition();
            Vector3 e2 = (*vi2)->getPosition() - (*vi1)->getPosition();
            sum_cross += e2.cross(e1);
        }

        setNormal(sum_cross.unit());
    } else {
        Vector3 e1 = (*vi0)->getPosition() - (*vi1)->getPosition();
        Vector3 e2 = (*vi2)->getPosition() - (*vi1)->getPosition();
        setNormal(e2.cross(e1).unit()); // counter-clockwise
    }
}

bool Face::contains(Vector3 const &p) const
{
    if (vertices.empty())
        return false;

    // Generate a ray for the even-odd test, from p to the midpoint of the first
    // halfedge. Ignore degenerate situations for
    // now.
    auto vi = vertices.begin();
    auto last = vi++;
    Vector3 u = 0.5 * ((*last)->getPosition() + (*vi)->getPosition()) - p;

    long count = 1; // first halfedge is obviously intersected, since we
                    // generated the ray through its midpoint
    for (; last != vertices.end(); ++vi) {
        if (vi == vertices.end())
            vi = vertices.begin();

        Vector3 v0 = (*last)->getPosition() - p;
        Vector3 v1 = (*vi)->getPosition() - p;

        // If winding order is: vector to first vertex, ray, vector to second
        // vertex, then intersects
        Vector3 c0 = v0.cross(u);
        Vector3 c1 = u.cross(v1);
        if (c0.dot(c1) > 0) // intersects, now check forward or reverse
        {
            // Forward if the vector to the point nearest to p on the line
            // containing the edge makes an acute angle with u.
            //
            // The point p' on line v + t * e closest to point p is v + t0 * e,
            // where t0 = e.dot(p - v) / e.dot(e)
            // (see www.geometrictools.com/Documentation/DistancePointLine.pdf).
            //
            // We translate p to the origin for simpler computations.
            Vector3 edge = v1 - v0;
            Real t0 = -edge.dot(v0) / edge.dot(edge);
            Vector3 u0 = v0 + t0 * edge;

            if (u0.dot(u) > 0)
                count++;
        }

        last = vi;
    }

    return (count % 2 == 1);
}

double Face::getAngle(Vertex *v)
{
    assert(hasVertex(v));
    Edge *e1 = NULL;
    Edge *e2 = NULL;
    for (auto it = v->edges.begin(); it != v->edges.end(); ++it) {
        if (hasEdge(*it)) {
            if (e1 == NULL)
                e1 = *it;
            else
                e2 = *it;
        }
    }

    Vector3 edge1, edge2;
    if (e1->getEndpoint(1)->getPosition() == v->getPosition())
        edge1 = e1->getEndpoint(0)->getPosition() - e1->getEndpoint(1)->getPosition();
    else
        edge1 = e1->getEndpoint(1)->getPosition() - e1->getEndpoint(0)->getPosition();

    if (e2->getEndpoint(1)->getPosition() == v->getPosition())
        edge2 = e2->getEndpoint(0)->getPosition() - e2->getEndpoint(1)->getPosition();
    else
        edge2 = e2->getEndpoint(1)->getPosition() - e2->getEndpoint(0)->getPosition();

    double angle = acos((edge1).dot(edge2)
                        / sqrt((edge1).squaredLength() * (edge2).squaredLength()));
    // FILL
    return angle;
}
