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
    // construct from vertex*, edge* + double(0 to 1), face*, xyz
    // set p and ptype accordingly

    vector<Edge *> get_visible_edges()
    {
        // FILL
        return {};
    }
};

#endif