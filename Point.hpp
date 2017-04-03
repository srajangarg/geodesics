#ifndef __Point_hpp__
#define __Point_hpp__

#include "Common.hpp"
#include "DGP/Vector3.hpp"
#include <vector>
using namespace std;

// Point on the surface of the mesh
class Point
{
public:
    enum PointType { VERTEX, EDGE, FACE, UNDEFINED };

    void *p;
    PointType ptype;
    Vector3 pos;

    // FILL
    // construct from 1. vertex*, 2. edge* + double(0 to 1) + vertex?,
    // 3. face*, 4. xyz coordinates
    // set p and ptype accordingly

    Point(Vertex *v)
    {
        p = v;
        pos = v->getPosition();
        ptype = VERTEX;
    };

    Point(Face *f)
    {
        p = f;
        for (auto it = f->vertices.begin(); it != f->vertices.end(); ++it)
            pos += (*it)->getPosition();
        pos /= f->vertices.size();
        ptype = FACE;
    };

    Point(Edge *e, double ratio = 0.5)
    {
        assert(ratio <= 1 and ratio >= 0);
        p = e;
        Vertex *v0 = e->getEndpoint(0);
        Vertex *v1 = e->getEndpoint(1);

        pos = ratio * v0->getPosition() + (1 - ratio) * v1->getPosition();
        ptype = EDGE;
    };

    Point(double x, double y, double z)
    {
        pos = Vector3(x, y, z);
        ptype = UNDEFINED;
    };

    vector<Edge *> get_visible_edges()
    {
        assert(ptype != UNDEFINED);

        vector<Edge *> visible;

        switch (ptype) {
            case VERTEX:
                for (auto it = ((Vertex *)p)->edges.begin();
                     it != ((Vertex *)p)->edges.end(); ++it)
                    visible.push_back(*it);
                break;

            case EDGE:
                visible.push_back((Edge *)p);
                break;

            case FACE:
                for (auto it = ((Face *)p)->edges.begin(); it != ((Face *)p)->edges.end();
                     ++it)
                    visible.push_back(*it);
        }
        return visible;
    }
};

#endif