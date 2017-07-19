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
    enum PointType { VERTEX, EDGE, UNDEFINED };

    void *p;
    PointType ptype;
    Vector3 pos;
    double ratio;

    Point()
    {
    }

    Point(Vertex *v)
    {
        p = v;
        pos = v->getPosition();
        ptype = VERTEX;
    }

    Point(Edge *e, double rat = 0.5)
    {
        assert(rat <= 1 and rat >= 0);

        ratio = rat;
        p = e;
        Vertex *v0 = e->getEndpoint(0);
        Vertex *v1 = e->getEndpoint(1);

        pos = (1 - rat) * v0->getPosition() + rat * v1->getPosition();
        ptype = EDGE;
    }

    Point(double x, double y, double z)
    {
        pos = Vector3(x, y, z);
        ptype = UNDEFINED;
    }

    bool operator<(const Point &rhs) const
    {
        if (p != rhs.p)
            return p < rhs.p;
        return pos < rhs.pos;
    }

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
                
            default:
                assert(false);
        }
        return visible;
    }

    friend ostream &operator<<(ostream &os, const Point &e)
    {
        if (e.ptype == VERTEX)
            os << "P(" << *((Vertex *)e.p) << ")";
        else
            os << "P(" << *((Edge *)e.p) << ", " << e.ratio << ")";
        return os;
    }
};

#endif