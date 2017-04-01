#include "Mesh.hpp"
#include "Point.hpp"
#include "DGP/Vector2.hpp"
#include <unordered_map>
using namespace std;

class Interval
{
public:

    Vector2 pos;
    double st, end, ps_d;
    double min_d;

    Edge* edge;
    Face* from;

    Interval(double x_, double y_, double st_, double end_, double ps_d_, Face* from_,
             Edge *edge_)
        : st(st_), end(end_), ps_d(ps_d_)
    {   
        pos = Vector2(x_, y_);
        edge = edge_;
        from = from_;

        assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
        recompute_min_d();
    }

    Interval(Vector2 pos_, double st_, double end_, double ps_d_, Face* from_,
             Edge *edge_)
        : pos(pos_), st(st_), end(end_), ps_d(ps_d_)
    {   
        edge = edge_;
        from = from_;
        
        assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
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

    bool operator<(const Interval &rhs) const
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
    unordered_map<Edge *, list<set<Interval>::iterator>> edge_intervals;
    set<Interval> intervals_heap;
    Point *source;
    vector<Point *> destinations;

    MMP(Mesh *m, Point *s, Point *d) : mesh(m), source(s)
    {
        // assert(s->ptype != UNDEFINED);
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

        Interval prop_w = *intervals_heap.begin();
        Edge* prop_e = prop_w.edge;

        intervals_heap.erase(intervals_heap.begin());

        for (auto &face : prop_e->faces)
        {
            if (prop_w.from == face)
                continue;

            Edge *e1 = face->getSuccessor(prop_e), *e2 = face->getSuccessor(e1);

            auto candids_1 = get_new_intervals(prop_w, e1);
            auto candids_2 = get_new_intervals(prop_w, e2);

            for (auto &new_w : candids_1)
                insert_new_interval(e1, new_w);
            for (auto &new_w : candids_2)
                insert_new_interval(e2, new_w);
        }
    }

    vector<Interval> source_bisect(double st, double end, Vector2 ps1, Vector2 ps2)
    {
        // FILL
        // takes in  a segment st-end and returns 1or2 interval
        // with the closest source indicated on each interval
    }

    vector<Interval> sanitize_and_merge(vector<Interval>& intervals)
    {
        // FILL

        return {};
    }

    void insert_new_interval(Edge *e, Interval &new_w)
    {
        // FILL
        // updates edge_intervals[e] and intervals according to algo discussed
        // should be O(edge_intervals[e])

        auto & intervals = edge_intervals[e];
        vector<Interval> new_intervals;
        bool new_pushed = false;

        auto it = intervals.begin();
        for (; it != intervals.end(); it++)
        {
            auto w = *it;

            if (new_pushed) break;

            if (w->st >= new_w.end)
            {
                // ----new----
                //               -----w-----
                new_intervals.push_back(new_w);
                new_intervals.push_back(*w);
                new_pushed = true;
            }
            else if (w->end <= new_w.st)
            {
                //             -----new-----
                // ----w----
                new_intervals.push_back(*w);
            }
            else if (w->st < new_w.st and w->end <= new_w.end)
            {
                //       -------new------
                // ------w--------

                Interval w_short(*w); w_short.end = new_w.st;

                new_intervals.push_back(w_short);
                for (auto &interval : source_bisect(new_w.st, w->end, w->pos, new_w.pos))
                    new_intervals.push_back(interval);

                new_w.st = w->end;
            }
            else if (w->st >= new_w.st and w->end > new_w.end)
            {
                // ------new------
                //        -----w--------

                Interval new_w_short(new_w); new_w_short.end = w->st;
                new_intervals.push_back(new_w_short);
                for (auto &interval : source_bisect(w->st, new_w.end, new_w.pos, w->pos))
                    new_intervals.push_back(interval);
                Interval w_short(*w); w_short.st = new_w.end;
                new_intervals.push_back(w_short);

                new_pushed = true;
            }
            else if (w->st >= new_w.st and w->end <= new_w.end)
            {
                // --------new----------
                //      -----w-----

                Interval new_w_short(new_w); new_w_short.end = w->st;
                new_intervals.push_back(new_w_short);
                for (auto &interval : source_bisect(w->st, w->end, new_w.pos, w->pos))
                    new_intervals.push_back(interval);
                new_w.st = w->end;
            }
            else if (w->st < new_w.st and w->end > new_w.end)
            {
                //     ----new----
                //  ----------w-------

                Interval w_short1(*w); w_short1.end = new_w.st;
                new_intervals.push_back(w_short1);
                for (auto &interval : source_bisect(new_w.st, new_w.end, new_w.pos, w->pos))
                    new_intervals.push_back(interval);
                Interval w_short2(*w); w_short2.st = new_w.end;
                new_intervals.push_back(w_short2);
                new_pushed = true;
            }
            else
            {
                // should never be reached
                assert(false);
            }
        }

        if (not new_pushed)
        {
            assert(it == intervals.end());
            new_intervals.push_back(new_w);
        }

        for (; it != intervals.end(); it++)
            new_intervals.push_back(**it);

        for (auto &it : intervals)
            intervals_heap.erase(it);
        intervals.clear();

        auto sanitized_new_intervals = sanitize_and_merge(new_intervals);

        for (auto interval : sanitized_new_intervals)
        {
            auto pp = intervals_heap.insert(interval);
            assert(pp.second);
            intervals.push_back(pp.first);
        }
    }

    vector<Interval> get_new_intervals(const Interval &w, Edge *e)
    {
        // FILL
        // given an interval, propogate it to an edge e, v is the vertex opposite to e
        // guaranteed that e is adjacent to the edge on which w lies
        // return vector of candidate intervals

        Edge* prop_e = w.edge;


        return {};
    }

    void initialize()
    {
        // FILL
        // insert edges which correspond to the source vertex
        // see `list_edges_visible_from_source` in `geodesics_algo_exact.h`
        // check mesh preconditions and sanity
        // all edges have 1/2 faces
        // each face has 3 edges, 3 vertices
    }

    void algorithm()
    {
        // FILL
        initialize();

        while(not intervals_heap.empty() and not terminate())
        {
            propagate();
        }

    }

    bool terminate()
    {
        // FILL
        return false;
    }

    // invariants
    // vector of intervals for each edge must be non overlapping and sorted by st
    // all valid windows MUST be in the set
};