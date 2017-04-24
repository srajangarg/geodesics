#include "Mesh.hpp"
#include "Point.hpp"
#include "DGP/Vector2.hpp"
#include <unordered_map>
#define epsilon 1e-15
using namespace std;

class Interval
{
public:
    Vector2 pos;
    double st, end, ps_d;
    double min_d;

    Edge *edge;
    Face *from;

    Interval(double x_, double y_, double st_, double end_, double ps_d_, Face *from_,
             Edge *edge_, bool invert)
        : ps_d(ps_d_)
    {
        edge = edge_;
        from = from_;
        
        //check if we need to invert the interval
        //INVARIANT - st is always close to endpoints[0] than end

        if(invert)
        {
            st = edge->length() - end_;
            end = edge->length() - st_;
            pos = Vector2(edge->length() - x_, y_);
        }
        else
        {
            st = st_;
            end = end_;
            pos = Vector2(x_, y_);
        }

        assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
        recompute_min_d();
    }

    Interval(Vector2 pos_, double st_, double end_, double ps_d_, Face *from_,
             Edge *edge_, bool invert)
        : pos(pos_), ps_d(ps_d_)
    {
        edge = edge_;
        from = from_;

        if(invert)
        {
            st = edge->length() - end_;
            end = edge->length() - st_;
            pos = Vector2(edge->length() - pos_.x(), pos_.y());
        }
        else
        {
            st = st_;
            end = end_;
        }

        assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
        recompute_min_d();
    }

    void recompute_min_d()
    {
        Vector2 st_v(st, 0);
        Vector2 end_v(end, 0);

        if (pos.x() >= st and pos.x() <= end)
            min_d = pos.y();
        else if (pos.x() < st)
            min_d = (pos - st_v).length();
        else
            min_d = (pos - end_v).length();
        
        min_d += ps_d;
    }

    double compute_max_d()
    {
        Vector2 st_v(st, 0);
        Vector2 end_v(end, 0);

        if (pos.x() >= (st + end) / 2.0)
            return (pos - st_v).length() + ps_d;
        else
            return (pos - end_v).length() + ps_d;
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
        assert(s->ptype != Point::UNDEFINED);
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
        Edge *prop_e = prop_w.edge;

        intervals_heap.erase(intervals_heap.begin());

        for (auto &face : prop_e->faces) {
            if (prop_w.from == face)
                continue;

            auto candidates = get_new_intervals(prop_w, face);

            for (auto &new_w : candidates)
                insert_new_interval(new_w);
        }
    }

    vector<Interval> source_bisect(double st, double end, const Interval &i1,
                                   const Interval &i2)
    {
        // takes in  a segment st-end and returns 1or2 interval
        // with the closest source indicated on each interval

        auto ps1 = i1.pos, ps2 = i2.pos;
        Vector2 st_v(st, 0), end_v(end, 0);

        vector<Interval> b_intervals;
        // FIX  INCOORECT double x = (ps2.squaredLength() -
        // ps1.squaredLength())/(2.0*(ps2.x() -
        // ps1.x()));
        double x = 0.0;
        // (x, 0) is equidistant from ps1 and ps2
        //calculated x directly from research paper
        double alpha = ps2.x() - ps1.x();
        double beta = i2.ps_d - i1.ps_d;
        double gamma = ps1.squaredLength() - ps2.squaredLength() - beta*beta;
        //ax^2 + bx + c = 0
        double a = alpha*alpha - beta*beta;
        //γα + 2s1xβ2
        double b = gamma*alpha + 2.0*ps2.x()*beta*beta;
        //c = 1/4γ2 − s12β2
        double c = 0.25*gamma*gamma - ps2.squaredLength()*beta*beta;

        bool equidist_pt = true;

        if(b*b - 4*a*c > 0)
        {
            double sq_det = sqrt(b*b - 4*a*c); 
            x = (-b + sq_det)/(2*a);
        }
        else
            equidist_pt = false;

        //equidist_pt is false when no such x exists, in that case one pseudo-source is closer

        if (x >= st and x <= end and equidist_pt) {
            // [st, x] closer to ps1 or ps2?
            if ((st_v - ps1).squaredLength() < (st_v - ps2).squaredLength()) {
                b_intervals.push_back(Interval(ps1, st, x, i1.ps_d, i1.from, i1.edge, false));
                b_intervals.push_back(Interval(ps2, x, end, i2.ps_d, i2.from, i2.edge, false));
            } else {
                b_intervals.push_back(Interval(ps2, st, x, i2.ps_d, i2.from, i2.edge, false));
                b_intervals.push_back(Interval(ps1, x, end, i1.ps_d, i1.from, i1.edge, false));
            }
        } else {
            // the complete interval is near to one single source
            if ((ps2 - st_v).squaredLength() < (ps1 - st_v).squaredLength())
                b_intervals.push_back(Interval(ps2, st, end, i2.ps_d, i2.from, i2.edge, false));
            else
                b_intervals.push_back(Interval(ps1, st, end, i1.ps_d, i1.from, i1.edge, false));
        }

        return b_intervals;
    }

    vector<Interval> sanitize_and_merge(vector<Interval> &intervals)
    {
        // FILL
        // remove too small intervals, merge intervals with same pos and ps_d etc.
        vector<Interval> sanitized_intervals;
        Interval last_merged_interval = sanitized_intervals[0];
        auto it1 = intervals.begin();
        auto it2 = it1;
        //last_merged_interval is interval corresponding to [it1, it2)
        for (it2++; it1 != intervals.end() and it2 != intervals.end();)
        {

            auto w = *it2;
            //check if w is too small
            if(w.end - w.st > epsilon)
            {
                //check if last merged interval and w have same pos and everything
                if(last_merged_interval.pos == w.pos and 
                    abs(last_merged_interval.ps_d - w.ps_d) < epsilon)
                {
                    //merge last_merge_interval and w
                    last_merged_interval.end = w.end;
                    last_merged_interval.recompute_min_d();
                    it2++;
                }
                else
                {
                    sanitized_intervals.push_back(last_merged_interval);
                    it1 = it2;
                    it2++;
                    last_merged_interval = *it1;
                }
            }
            else
                it2++;
        }
        sanitized_intervals.push_back(last_merged_interval);
        return sanitized_intervals;
    }

    void insert_new_interval(Interval &new_w)
    {
        // FILL
        // updates edge_intervals[e] and intervals according to algo discussed
        // should be O(edge_intervals[e])

        auto &intervals = edge_intervals[new_w.edge];
        vector<Interval> new_intervals;
        bool new_pushed = false;

        auto it = intervals.begin();
        for (; it != intervals.end(); it++) {
            auto w = *it;

            if (new_pushed)
                break;

            if (w->st >= new_w.end) {
                // ----new----
                //               -----w-----
                new_intervals.push_back(new_w);
                new_intervals.push_back(*w);
                new_pushed = true;
            } else if (w->end <= new_w.st) {
                //             -----new-----
                // ----w----
                new_intervals.push_back(*w);
            } else if (w->st < new_w.st and w->end <= new_w.end) {
                //       -------new------
                // ------w--------

                Interval w_short(*w);
                w_short.end = new_w.st;

                new_intervals.push_back(w_short);
                for (auto &interval : source_bisect(new_w.st, w->end, *w, new_w))
                    new_intervals.push_back(interval);

                new_w.st = w->end;
            } else if (w->st >= new_w.st and w->end > new_w.end) {
                // ------new------
                //        -----w--------

                Interval new_w_short(new_w);
                new_w_short.end = w->st;
                new_intervals.push_back(new_w_short);
                for (auto &interval : source_bisect(w->st, new_w.end, new_w, *w))
                    new_intervals.push_back(interval);
                Interval w_short(*w);
                w_short.st = new_w.end;
                new_intervals.push_back(w_short);

                new_pushed = true;
            } else if (w->st >= new_w.st and w->end <= new_w.end) {
                // --------new----------
                //      -----w-----

                Interval new_w_short(new_w);
                new_w_short.end = w->st;
                new_intervals.push_back(new_w_short);
                for (auto &interval : source_bisect(w->st, w->end, new_w, *w))
                    new_intervals.push_back(interval);
                new_w.st = w->end;
            } else if (w->st < new_w.st and w->end > new_w.end) {
                //     ----new----
                //  ----------w-------

                Interval w_short1(*w);
                w_short1.end = new_w.st;
                new_intervals.push_back(w_short1);
                for (auto &interval : source_bisect(new_w.st, new_w.end, new_w, *w))
                    new_intervals.push_back(interval);
                Interval w_short2(*w);
                w_short2.st = new_w.end;
                new_intervals.push_back(w_short2);
                new_pushed = true;
            } else {
                // should never be reached
                assert(false);
            }
        }

        if (not new_pushed) {
            assert(it == intervals.end());
            new_intervals.push_back(new_w);
        }

        for (; it != intervals.end(); it++)
            new_intervals.push_back(**it);

        for (auto &it : intervals)
            intervals_heap.erase(it);
        intervals.clear();

        auto sanitized_new_intervals = sanitize_and_merge(new_intervals);

        for (auto interval : sanitized_new_intervals) {
            interval.recompute_min_d();
            auto pp = intervals_heap.insert(interval);
            assert(pp.second);
            intervals.push_back(pp.first);
        }
    }

    vector<Interval> get_new_intervals(const Interval &w, Face *face)
    {
        // FILL
        // given an interval, propogate it to an edge e, v is the vertex opposite to e
        // guaranteed that e is adjacent to the edge on which w lies
        // return vector of candidate intervals

        Edge *prop_e = w.edge;

        return {};
    }

    bool check_mesh_sanity()
    {
        // FILL
        // check mesh preconditions and sanity
        // all edges have 1/2 faces
        // each face has 3 edges, 3 vertices
        return false;
    };

    void initialize()
    {
        // FILL
        assert(check_mesh_sanity());

        switch (source->ptype) {
            case Point::VERTEX:

                // FIX should it be 0.0 or edgelength????
                // invert wale pains???

                for (auto &e : ((Vertex *)(source->p))->edges) {
                    //check for invert
                    //Interval(double x_, double y_, double st_, double end_, double ps_d_, Face *from_,
                    // Edge *edge_, bool invert)
                    ////on which side of e does p lie?
                    bool invert = false;
                    if((e->getEndpoint(0) == source->p and e->getEndpoint(1) < e->getEndpoint(0)) or
                        (e->getEndpoint(1) == source->p and e->getEndpoint(0) < e->getEndpoint(1)))
                        invert = true;
                    Interval ii(0, 0, 0, e->length(), 0, NULL, e, invert);
                    insert_new_interval(ii);
                }
                break;

            case Point::EDGE:
            {
                //source->p is edge *
                Edge * e = (Edge *)source->p;
                bool invert = e->getEndpoint(0) > e->getEndpoint(1);
                //x = (pos - endpoint[0]).length()
                //endpoint[0]---------pos------------endpoint[1]
                Interval ii((source->pos - e->getEndpoint(0)->getPosition()).length(), 0, 0, e->length(), 0, NULL, e, invert);
                insert_new_interval(ii);
                break;
            }
            case Point::FACE:
                for (auto &e : ((Face *)(source->p))->edges)
                {
                    Vector3 pos1 = e->getEndpoint(0)->getPosition();
                    Vector3 pos2 = e->getEndpoint(1)->getPosition();
                    //         pos(y)
                    //       /
                    //     /
                    //(0)pos1--(x)-----------pos2
                    double x = (source->pos - pos1).dot((pos2 - pos1).unit());
                    double y = sqrt((source->pos - pos1).squaredLength() - x*x);
                    //pythagoras theoram
                    //Interval(double x_, double y_, double st_, double end_, double ps_d_, Face *from_,
                    // Edge *edge_, bool invert)
                    bool invert = e->getEndpoint(0) > e->getEndpoint(1);
                    Interval ii(x, y, 0, e->length(), 0, (Face *)(source->p), e, invert);
                    insert_new_interval(ii);
                }   
                break;
        }
    }

    void algorithm()
    {
        // FILL
        initialize();

        while (not intervals_heap.empty() and not terminate()) {
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