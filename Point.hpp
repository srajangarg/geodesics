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

    // FILL
    // construct from 1. vertex*, 2. edge* + double(0 to 1) + vertex?,
    // 3. face*, 4. xyz coordinates
    // set p and ptype accordingly

    vector<Edge *> get_visible_edges()
    {
        // FILL
        assert(ptype != UNDEFINED);

        vector<Edge *> visible;

        switch (ptype) {
            case VERTEX:

                break;

            case EDGE:

                break;

            case FACE:

                break;
        }
        return visible;
    }
};

#endif