#include "Mesh.hpp"
#include "Point.hpp"
#include "DGP/Vector2.hpp"
#include <unordered_map>
#include <cmath>
#define EPS 1e-6
using namespace std;

class Interval
{
public:
    Vector2 pos;
    double st, end, ps_d;
    double min_d;

    Edge *edge;
    Face *from;
    Interval *parent;

    Interval()
    {
    }

    void set_st_end_pos(double st_, double end_, bool invert);

    Interval(double x_, double y_, double st_, double end_, double ps_d_, Face *from_,
             Edge *edge_, Interval *par, bool invert = false);

    Interval(Vector2 pos_, double st_, double end_, double ps_d_, Face *from_,
             Edge *edge_, Interval *par, bool invert = false);

    Interval(Vector2 pos_, double st_, double end_, const Interval &i,
             bool invert = false);

    void recompute_min_d();
    double compute_max_d();
    bool operator<(const Interval &rhs) const;
    bool operator==(const Interval &rhs) const;

    struct Info {
        double angle, r;
        double e1, e2;
        Vertex *common;
    };

    Info get_info_edge(Vector2 src, Edge *e, Face *face) const
    {
        Vertex *v = e->getCommonVertex(*edge);
        assert(v != NULL);

        Info i;
        i.angle = face->getAngle(v);
        i.common = v;
        i.e1 = 0;
        i.e2 = std::numeric_limits<double>::infinity();
        double ost, oend;

        if (v == edge->getEndpoint(0))
            i.r = src.x(), ost = st, oend = end;
        else
            i.r = edge->length() - src.x(), ost = edge->length() - st,
            oend = edge->length() - end;

        if (src.y() < EPS)
            return i;

        double x = i.r, y = src.y();
        double angle1 = i.angle, angle2 = atan2(y, x);
        i.r = Vector2(x, y).length();
        i.angle += angle2;

        double num = ost * y;
        double den = sin(angle1) * (x - ost) + cos(angle1) * y;

        if (abs(den) < EPS or num / den < 0)
            i.e1 = std::numeric_limits<double>::infinity();
        else
            i.e1 = num / den;

        num = oend * y;
        den = sin(angle1) * (x - oend) + cos(angle1) * y;

        if (abs(den) < EPS or num / den < 0)
            i.e2 = std::numeric_limits<double>::infinity();
        else
            i.e2 = num / den;

        if (v == edge->getEndpoint(1))
            swap(i.e1, i.e2);

        return i;
    }

    // y >= 0
    // recompute_min_d() must be called if any of the 5 parameters change
    friend ostream &operator<<(ostream &os, const Interval &e);
};

class MMP
{
public:
    Mesh *mesh;
    unordered_map<Edge *, list<Interval *>> edge_intervals;
    Point source;
    vector<Point> destinations;

    struct IntervalPtrComp {
        bool operator()(const Interval *lhs, const Interval *rhs) const
        {
            return *lhs < *rhs;
        }
    };
    set<Interval *, IntervalPtrComp> intervals_heap;

    MMP()
    {
    }

    MMP(Mesh *m, Point s, Point d) : mesh(m), source(s)
    {
        assert(s.ptype != Point::UNDEFINED);
        destinations.push_back(d);
    }

    MMP(Mesh *m, Point s, const vector<Point> &dests) : mesh(m), source(s)
    {
        destinations = dests;
    }

    vector<Interval> prop_thru_interval(Vector2 src, Interval &w, Face *face, double ps_d)
    {
        // face is the face on which the intervals propogate
        // check whether src lies on the interval
        vector<Interval> new_intervals;
        if (src.y() < EPS and (src.x() < w.st or src.x() > w.end))
            return {};

        for (auto &edge : face->edges) {
            if (edge == w.edge)
                continue;

            auto info = w.get_info_edge(src, edge, face);
            bool invert = info.common > edge->getOtherEndpoint(info.common);

            if (info.e1 < edge->length()) {

                double x = info.r * cos(info.angle);
                if (abs(x) < EPS)
                    x = 0;

                double y = abs(info.r * sin(info.angle));
                if (y < EPS)
                    y = 0;

                new_intervals.push_back(Interval(x, y, info.e1,
                                                 min(info.e2, edge->length()), ps_d, face,
                                                 edge, &w, invert));
            }
            // otherwise doesn't intersect this edge
        }

        return new_intervals;
    }

    void propagate()
    {
        Interval prop_w = **intervals_heap.begin();
        Edge *prop_e = prop_w.edge;

        intervals_heap.erase(intervals_heap.begin());

        for (auto &face : prop_e->faces) {
            if (prop_w.from == face)
                continue;

            cout << "Propagating " << prop_w << " on " << *face << endl;
            auto candidates = get_new_intervals(prop_w, face);

            cout << endl;
            for (auto &new_w : candidates)
                insert_new_interval(new_w);
        }
    }

    vector<Point> trace_back(Point destination)
    {
        assert(false);
        // vector<Point> path;
        // assert(destination.ptype != Point::UNDEFINED);

        // path.push_back(destination);

        // Interval best_interval = best_first_interval(destination, path);

        // from the given destination and best_interval, find the next intervals and point
        // accordingly
    }

    vector<Interval> source_bisect(double st, double end, const Interval &i1,
                                   const Interval &i2)
    {
        // takes in  a segment st-end and returns 1or2 interval
        // with the closest source indicated on each interval
        auto ps1 = i1.pos, ps2 = i2.pos;
        Vector2 st_v(st, 0), end_v(end, 0);
        vector<Interval> b_intervals;
        double x, alpha, beta, gamma;

        alpha = ps2.x() - ps1.x();
        beta = i2.ps_d - i1.ps_d;
        gamma = ps1.squaredLength() - ps2.squaredLength() - beta * beta;

        // ax^2 + bx + c = 0
        double a = alpha * alpha - beta * beta;
        double b = gamma * alpha + 2.0 * ps2.x() * beta * beta;
        double c = 0.25 * gamma * gamma - ps2.squaredLength() * beta * beta;
        bool equidist_pt = true;

        if (b * b - 4 * a * c < 0 or (a == 0.0 and b == 0.0 and abs(c) < EPS))
            equidist_pt = false;
        else if (abs(a) < EPS)
            x = -c / b;
        else
            x = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);

        if (x > st and x < end and equidist_pt) {
            // [st, x] closer to ps1 or ps2?
            if ((st_v - ps1).squaredLength() < (st_v - ps2).squaredLength()) {
                b_intervals.push_back(Interval(ps1, st, x, i1));
                b_intervals.push_back(Interval(ps2, x, end, i2));
            } else {
                b_intervals.push_back(Interval(ps2, st, x, i2));
                b_intervals.push_back(Interval(ps1, x, end, i1));
            }
        } else {
            // the complete interval is near to one single source
            if (((ps2 - st_v).length() + i2.ps_d) < ((ps1 - st_v).length() + i1.ps_d))
                b_intervals.push_back(Interval(ps2, st, end, i2));
            else
                b_intervals.push_back(Interval(ps1, st, end, i1));
        }

        return b_intervals;
    }

    vector<Interval> sanitize_and_merge(vector<Interval> &intervals)
    {
        // remove too small intervals, merge intervals with same pos and ps_d etc.
        vector<Interval> sanitized_intervals;
        Interval last_merged_interval = intervals[0];
        auto it1 = intervals.begin();
        auto it2 = it1;
        // last_merged_interval is interval corresponding to [it1, it2)
        for (it2++; it1 != intervals.end() and it2 != intervals.end();) {

            auto w = *it2;
            // check if w is too small
            if (w.end - w.st > EPS) {
                // check if last merged interval and w have same pos and everything
                if (last_merged_interval.pos == w.pos
                    and abs(last_merged_interval.ps_d - w.ps_d) < EPS) {
                    // merge last_merge_interval and w
                    last_merged_interval.end = w.end;
                    last_merged_interval.recompute_min_d();
                    it2++;
                } else {

                    if (abs(last_merged_interval.end - last_merged_interval.st) > EPS)
                        sanitized_intervals.push_back(last_merged_interval);

                    it1 = it2;
                    it2++;
                    last_merged_interval = *it1;
                }
            } else
                it2++;
        }
        sanitized_intervals.push_back(last_merged_interval);
        return sanitized_intervals;
    }

    void insert_new_interval(Interval &new_w)
    {
        // updates edge_intervals[e] and intervals according to algo discussed
        // should be O(edge_intervals[e])
        cout << "Inserting new interval " << new_w << endl;

        auto &intervals = edge_intervals[new_w.edge];
        vector<Interval> new_intervals;
        bool new_fully_pushed = false;

        auto it = intervals.begin();
        for (; it != intervals.end(); it++) {
            auto w = *it;

            if (new_fully_pushed)
                break;

            if (w->st >= new_w.end) {
                // ----new----
                //               -----w-----
                new_intervals.push_back(new_w);
                new_intervals.push_back(*w);
                new_fully_pushed = true;
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

                new_fully_pushed = true;
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
                new_fully_pushed = true;
            } else {
                // should never be reached
                assert(false);
            }
        }

        if (not new_fully_pushed) {
            assert(it == intervals.end());
            new_intervals.push_back(new_w);
        }

        for (; it != intervals.end(); it++)
            new_intervals.push_back(**it);

        for (auto &it : intervals)
            intervals_heap.erase(it);

        set<Interval*> not_to_delete;
        auto old_intervals = intervals;
        intervals.clear();
        
        auto sanitized_new_intervals = sanitize_and_merge(new_intervals);
        cout << "Old intervals cleared, all intervals on "<<*new_w.edge<<" are :"<< endl;

        for (auto &interval : sanitized_new_intervals) {

            Interval* keep_old_itv = NULL;
            for (auto &old_itv : old_intervals)
                if (interval == *old_itv)
                    keep_old_itv = old_itv;

            interval.recompute_min_d();
            cout << interval << endl;

            if (keep_old_itv)
            {
                intervals.push_back(keep_old_itv);
                not_to_delete.insert(keep_old_itv);
            }
            else
            {
                auto added_intv = new Interval(interval);
                auto pp = intervals_heap.insert(added_intv);
                assert(pp.second);
                intervals.push_back(added_intv);
            }
        }
        cout << endl;
    
        for (auto &itv : old_intervals)
            if (not_to_delete.find(itv) == not_to_delete.end())
                delete itv;
    }


    vector<Interval> get_new_intervals(Interval &w, Face *face)
    {
        // propogate this interval to edges on `face`
        Edge *edge = w.edge;
        vector<Interval> new_intervals;
        Vertex *v0 = edge->getEndpoint(0), *v1 = edge->getEndpoint(1);

        if (w.st < EPS and v0->saddle_or_boundary) {
            // treat v0 as a new source

            Vector2 new_src = Vector2(0, 0);
            auto intervals = prop_thru_interval(new_src, w, face,
                                                w.ps_d + (new_src - w.pos).length());

            if (intervals.size())
                cout << intervals.size() << " to be added for saddle " << *v0 << ":"
                     << endl;

            for (auto &ii : intervals)
                cout << ii << endl, new_intervals.push_back(ii);
        }

        if (edge->length() - w.end < EPS and v1->saddle_or_boundary) {
            // treat v1 as a new source

            Vector2 new_src = Vector2(edge->length(), 0);

            auto intervals = prop_thru_interval(new_src, w, face,
                                                w.ps_d + (new_src - w.pos).length());
            if (intervals.size())
                cout << intervals.size() << " to be added for saddle " << *v1 << ":"
                     << endl;

            for (auto &ii : intervals)
                cout << ii << endl, new_intervals.push_back(ii);
        }

        // generate intervals from (usual) current source
        auto intervals = prop_thru_interval(w.pos, w, face, w.ps_d);

        if (intervals.size())
            cout << intervals.size() << " to be added for usual propogate:" << endl;
        for (auto &ii : intervals)
            cout << ii << endl, new_intervals.push_back(ii);
        return new_intervals;
    }

    bool check_mesh_sanity()
    {
        for (auto &edge : mesh->edges)
            if (edge.faces.size() == 0 or edge.faces.size() > 3)
                return false;

        for (auto &face : mesh->faces)
            if (face.edges.size() != 3 or face.vertices.size() != 3)
                return false;

        return true;
    };

    void initialize()
    {
        assert(check_mesh_sanity());
        auto visible = source.get_visible_edges();

        cout << "Making initial intervals ..." << endl;
        switch (source.ptype) {
            case Point::VERTEX: {
                for (auto &e : visible) {
                    Interval ii(0, 0, 0, e->length(), 0, NULL, e, NULL,
                                (source.p == e->getEndpoint(1)));
                    insert_new_interval(ii);
                }
                break;
            }

            case Point::EDGE: {
                auto e = (Edge *)source.p;
                Interval ii(source.ratio * e->length(), 0, 0, e->length(), 0, NULL, e,
                            NULL, false);
                insert_new_interval(ii);
                break;
            }

            case Point::FACE: {
                for (auto &e : visible) {
                    Vector3 pos1 = e->getEndpoint(0)->getPosition();
                    Vector3 pos2 = e->getEndpoint(1)->getPosition();
                    double x = (source.pos - pos1).dot((pos2 - pos1).unit());
                    double y = sqrt((source.pos - pos1).squaredLength() - x * x);

                    Interval ii(x, y, 0, e->length(), 0, (Face *)(source.p), e, NULL,
                                false);
                    insert_new_interval(ii);
                }
                break;
            }

            default:
                assert(false);
        }
    }

    bool best_first_interval(Point dest, Interval &bf)
    {
        switch (dest.ptype) {
            case Point::VERTEX: {

                auto v = (Vertex *)dest.p;
                double mind = std::numeric_limits<double>::infinity();
                for (auto &e : v->edges) {
                    Interval *ii = NULL;
                    double dis;

                    if (edge_intervals[e].empty())
                        continue;

                    if (v == e->getEndpoint(0)) {
                        ii = edge_intervals[e].front();
                        if (ii->st > EPS)
                            continue;
                        dis = ii->pos.length() + ii->ps_d;
                    } else {
                        ii = edge_intervals[e].back();
                        if (e->length() - ii->end > EPS)
                            continue;
                        dis = (ii->pos - Vector2(e->length(), 0)).length() + ii->ps_d;
                    }

                    if (dis < mind) {
                        bf = *ii;
                        mind = dis;
                    }
                }

                if (mind != std::numeric_limits<double>::infinity())
                    return true;
                break;
            }

            case Point::EDGE: {

                auto e = (Edge *)dest.p;
                for (auto &itv : edge_intervals[e])
                    if ((itv->st / e->length() >= dest.ratio)
                        and (itv->end / e->length() <= dest.ratio)) {
                        bf = *itv;
                        return true;
                    }
                break;
            }

            default:
                assert(false);
        }

        return false;
    }

    void algorithm()
    {
        initialize();
        while (not intervals_heap.empty() and not terminate()) {
            propagate();
        }
    }

    bool terminate()
    {
        // Interval ii;
        // for (auto &dest : destinations)
        //     if (not best_first_interval(dest, ii))
        //         return false;
        // return true;
        return false;
    }

    // invariants
    // vector of intervals for each edge must be non overlapping and sorted by st
    // all valid windows MUST be in the set
};