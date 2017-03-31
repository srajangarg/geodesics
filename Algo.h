#include "Mesh.hpp"
#include "Point.hpp"
#include <unordered_map>
using namespace std;

class Interval
{
public:
    enum Direction { FROM_SOURCE, FROM_LEFT, FROM_RIGHT, UNDEFINED };

    double x, y, st, end, ps_d;
    double min_d;

    Edge *edge;
    Direction from;

    Interval(double x_, double y_, double st_, double end_, double ps_d_, Direction from_,
             Edge *edge_)
        : x(x_), y(y_), st(st_), end(end_), ps_d(ps_d_)
    {
        edge = edge_;
        from = from_;
        recompute_min_d();
    }

    void recompute_min_d()
    {
        // FILL
    }

    double compute_max_d()
    {
        // FILL
        return 0.0;
    }

    bool operator<(const Interval &rhs)
    {
        if (min_d != rhs.min_d)
            return min_d < rhs.min_d;
        return (edge < rhs.edge);
    }

    // invariants
    // y >= 0
    // recompute_min_d() must be called if any of the 5 parameters change
};

class MMP
{
public:
    Mesh *mesh;
    unordered_map<Edge *, vector<Interval *>> edge_intervals;
    set<Interval> intervals;
    Point *source;
    vector<Point *> destinations;

    MMP(Mesh *m, Point *s, Point *d) : mesh(m), source(s)
    {
        assert(s->ptype != UNDEFINED);
        destinations.push_back(d);
    }

    MMP(Mesh *m, Point *s, const vector<Point *> &dests) : mesh(m), source(s)
    {
        destinations = dests;
    }

    void propagate()
    {
        // FILL
        // removes first element of set, and propagates accordingly once
    }

    void insert_new_interval(Edge *e, const Interval &w)
    {
        // FILL
        // updates edge_intervals[e] and intervals according to algo discussed
        // should be O(edge_intervals[e])
    }

    vector<Interval> get_new_intervals(const Interval &w, Edge *e, Vertex *v)
    {
        // FILL
        // given an interval, propogate it to an edge e, v is the vertex opposite to e
        // guaranteed that e is adjacent to the edge on which w lies
        // return vector of candidate intervals
        return {};
    }

    void initialize()
    {
        // FILL
        // insert edges which correspond to the source vertex
    }

    void algorithm()
    {
        // FILL
    }

    bool check_termination()
    {
        // FILL
        return false;
    }

    // invariants
    // vector of intervals for each edge must be non overlapping and sorted by st
};