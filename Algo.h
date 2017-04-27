#include "Mesh.hpp"
#include "Point.hpp"
#include "DGP/Vector2.hpp"
#include <unordered_map>
#define EPS 1e-15
using namespace std;

class Interval
{
public:
    Vector2 pos;
    double st, end, ps_d;
    double min_d;

    Edge *edge;
    Face *from;

    Interval()
    {
    }

    void set_st_end_pos(double st_, double end_, bool invert)
    {
        if (invert) {
            st = edge->length() - end_;
            end = edge->length() - st_;
            pos = Vector2(edge->length() - pos.x(), pos.y());
        } else {
            st = st_;
            end = end_;
            pos = Vector2(pos.x(), pos.y());
        }
    }

    Interval(double x_, double y_, double st_, double end_, double ps_d_, Face *from_,
             Edge *edge_, bool invert = false)
    {
        edge = edge_;
        from = from_;
        ps_d = ps_d_;
        pos = Vector2(x_, y_);
        set_st_end_pos(st_, end_, invert);
        assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
        // INVARIANT - st is always close to lower endpoint than higher endpoint pointer
        recompute_min_d();
    }

    Interval(double x_, double y_, double st_, double end_, const Interval &i,
             bool invert = false)
    {
        edge = i.edge;
        from = i.from;
        ps_d = i.ps_d;
        pos = Vector2(x_, y_);
        set_st_end_pos(st_, end_, invert);

        assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
        // INVARIANT - st is always close to lower endpoint than higher endpoint pointer
        recompute_min_d();
    }

    Interval(Vector2 pos_, double st_, double end_, double ps_d_, Face *from_,
             Edge *edge_, bool invert = false)
    {
        edge = edge_;
        from = from_;
        pos = pos_;
        ps_d = ps_d_;
        set_st_end_pos(st_, end_, invert);

        assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
        recompute_min_d();
    }

    Interval(Vector2 pos_, double st_, double end_, const Interval &i,
             bool invert = false)
    {
        edge = i.edge;
        from = i.from;
        ps_d = i.ps_d;
        pos = pos_;
        set_st_end_pos(st_, end_, invert);

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

    struct Info {
        double angle;
        double x, y;
        Vertex *common;
    };

    Info get_info_edge(Edge *e) const
    {
        Vertex *v = e->getCommonVertex(*edge);
        assert(v != NULL);

        Info i;
        i.angle = from->getAngle(v);
        i.y = pos.y();
        i.common = v;

        if (v < edge->getOtherEndpoint(v))
            i.x = pos.x();
        else
            i.x = edge->length() - pos.x();

        return i;
    }

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

    Interval best_first_interval(Point destination, vector<Point> &path)
    {
        if (destination.ptype == Point::EDGE) {
            Edge *e = (Edge *)destination.p;
            auto &intervals = edge_intervals[e];
            double ratio = (destination.pos - e->getEndpoint(0)->getPosition()).length()
                           / e->length();

            for (auto interval : intervals) {
                if ((interval->st / e->length() > ratio)
                    and (interval->end / e->length() < ratio))
                    return *interval;
            }
        } else if (destination.ptype == Point::FACE) {
            Face *f = (Face *)destination.p;
            double distance = std::numeric_limits<double>::infinity();
            Interval min_interval;
            Point next_point(0, 0, 0); // temporarily an undefined point
            for (auto it = f->edges.begin(); it != f->edges.end(); ++it) {
                // pos += (*it)->getPosition();
                auto &intervals = edge_intervals[*it];
                // find the interval which is closest to the destination and source
                for (auto interval : intervals) {
                    if (interval->from == f)
                        continue;

                    // consider only such intervals through which the path can pass
                    // through
                    Vector3 pos = destination.pos;
                    Vector3 pos1 = (*it)->getEndpoint(0)->getPosition();
                    Vector3 pos2 = (*it)->getEndpoint(1)->getPosition();
                    //         pos(source)
                    //       /
                    //     /
                    //(0)pos1--(st,0)------(end,0)---pos2
                    ////////////////////////
                    //------------pos(destination)
                    if ((*it)->getEndpoint(0) > (*it)->getEndpoint(1))
                        swap(pos1, pos2);

                    double dest_x = (pos - pos1).dot((pos2 - pos1).unit());
                    double dest_y = -sqrt((destination.pos - pos1).squaredLength()
                                          - dest_x * dest_x);
                    // pythagoras theoram
                    // get the point where source.pos and dest.pos meet x-axis and
                    // then calculate the miimum distance to source through this interval
                    double axis_x
                        = interval->pos.x()
                          - interval->pos.y() * ((interval->pos.x() - dest_x)
                                                 / (interval->pos.y() - dest_y));

                    double source_dist = interval->ps_d;
                    double ratio;
                    // minimum distance from dest to source through this interval
                    if (axis_x > interval->st and axis_x < interval->end) {
                        source_dist += (Vector2(dest_x, dest_y) - interval->pos).length();
                        ratio = axis_x / (*it)->length();
                    } else if (axis_x < interval->st) {
                        source_dist
                            += ((Vector2(dest_x, dest_y) - Vector2(interval->st, 0.0))
                                    .length()
                                + (interval->pos - Vector2(interval->st, 0.0)).length());
                        ratio = interval->st / (*it)->length();
                    } else {
                        source_dist
                            += ((Vector2(dest_x, dest_y) - Vector2(interval->end, 0.0))
                                    .length()
                                + (interval->pos - Vector2(interval->end, 0.0)).length());
                        ratio = interval->end / (*it)->length();
                    }

                    if (source_dist < distance) {
                        distance = source_dist;
                        min_interval = *interval;
                        next_point = Point(*it, ratio);
                    }
                }
            }

            path.push_back(next_point);
            return min_interval;
        } else if (destination.ptype == Point::VERTEX) {
            Vertex *v = (Vertex *)destination.p;
            double distance = std::numeric_limits<double>::infinity();
            Interval min_interval;
            for (auto it = v->edges.begin(); it != v->edges.end(); ++it) {
                Interval interval;
                double source_dist = interval.ps_d;
                if (v == (*it)->getEndpoint(0)) {
                    interval = *(edge_intervals[*it].front());
                    source_dist += (interval.pos - Vector2(interval.st, 0.0)).length();
                    source_dist += interval.st;
                } else {
                    interval = *(edge_intervals[*it].back());
                    source_dist += (interval.pos - Vector2(interval.end, 0.0)).length();
                    source_dist += ((*it)->length() - interval.end);
                }
                if (source_dist < distance)
                    min_interval = interval;
            }

            return min_interval;
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

        if (b * b - 4 * a * c > 0)
            x = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
        else
            equidist_pt = false;

        // equidist_pt is false when no such x exists, in that case one pseudo-source is
        // closer

        if (x >= st and x <= end and equidist_pt) {
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
            if ((ps2 - st_v).squaredLength() < (ps1 - st_v).squaredLength())
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
        Interval last_merged_interval = sanitized_intervals[0];
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

    vector<Interval> get_intervals_source(Vector2 src, const Interval &w, Face *face,
                                          double ps_d)
    {
        // face is the face on which the intervals propogate
        // check whether src lies on the interval
        vector<Interval> intervals_source;

        if (src.y() < EPS and src.x() >= w.st and src.x() <= w.end) {
            // take two intervals on the other two edges of the face
            for (auto &edge : face->edges) {
                if (edge == w.edge)
                    continue;

                auto info = w.get_info_edge(edge);
                bool invert = info.common > edge->getOtherEndpoint(info.common);
                Interval ii(info.x * cos(info.angle), info.x * sin(info.angle), 0,
                            edge->length(), ps_d, face, edge, invert);
                intervals_source.push_back(ii);
            }
        } else {
            // draw lines and solve this
            for (auto &edge : face->edges) {
                if (edge == w.edge)
                    continue;

                Vertex *v0 = edge->getEndpoint(0) < edge->getEndpoint(1)
                                 ? edge->getEndpoint(0)
                                 : edge->getEndpoint(1);
                Vertex *v1 = edge->getOtherEndpoint(v0);

                double theta1, theta2, r, src_x, src_y;
                //   src
                //  /
                // / <theta2>
                // v0----------v1
                // \<theta1>
                //  \
                //   edge

                if (edge->hasEndpoint(v0)) {
                    theta1 = face->getAngle(v0);
                    theta2 = atan(src.y() / src.x());
                    r = sqrt(src.x() * src.x() + src.y() * src.y());
                    src_x = src.x();
                } else if (edge->hasEndpoint(v1)) {
                    theta1 = face->getAngle(v1);
                    theta2 = atan(src.y() / (edge->length() - src.x()));
                    r = sqrt((edge->length() - src.x()) * (edge->length() - src.x())
                             + src.y() * src.y());
                    src_x = edge->length() - src.x();
                }
                src_y = src.y();

                double theta = theta1 + theta2;

                double e1 = (w.st * src_y)
                            / (sin(theta1) * (src_x - w.st) + cos(theta1) * src_y);
                double e2 = (w.end * src_y)
                            / (sin(theta1) * (src_x - w.end) + cos(theta1) * src_y);

                if (edge->hasEndpoint(v1))
                    swap(e1, e2);

                // Interval(double x_, double y_, double st_, double end_, double ps_d_,
                // Face *from_,
                // Edge *edge_, bool invert)
                bool invert = edge->getEndpoint(0) > edge->getEndpoint(1);

                if (e1 >= 0 and e2 < 0 and e1 <= edge->length()) {
                    e2 = edge->length();
                    Interval ii(r * cos(theta), r * sin(theta), e1, e2, ps_d, face, edge,
                                invert);
                    intervals_source.push_back(ii);
                } else if (e1 >= 0 and e2 >= 0 and e1 <= edge->length()) {
                    if (e2 > edge->length())
                        e2 = edge->length();

                    Interval ii(r * cos(theta), r * sin(theta), e1, e2, ps_d, face, edge,
                                invert);
                    intervals_source.push_back(ii);
                }
                // otherwise doesn't intersect this edge
            }
        }

        return intervals_source;
    }

    vector<Interval> get_new_intervals(const Interval &w, Face *face)
    {
        // propogate to the face opposite to face
        // insert intervals done in propogate, just calculate the intervals here
        // vector<Interval> get_intervals_source(Vector2 src, const Interval &w, Face *
        // face, double ps_d)

        Edge *prop_e = w.edge;

        vector<Interval> new_intervals;

        for (auto &e_face : prop_e->faces) {
            if (e_face != face)
                continue;

            // propogate only on opposite side of source
            // check if any end of the interval is endpoint
            Vertex *v0 = prop_e->getEndpoint(0) < prop_e->getEndpoint(1)
                             ? prop_e->getEndpoint(0)
                             : prop_e->getEndpoint(1);
            Vertex *v1 = prop_e->getOtherEndpoint(v0);
            if (w.st < EPS and v0->saddle_or_boundary) {
                // treat this as a new source
                // pass this to a function
                for (auto ii : get_intervals_source(
                         Vector2(0.0, 0.0), w, e_face,
                         w.ps_d + sqrt(w.pos.x() * w.pos.x() + w.pos.y() * w.pos.y())))
                    new_intervals.push_back(ii);
            }
            if (abs(w.end - prop_e->length()) < EPS and v1->saddle_or_boundary) {
                // treat this as a new source
                // pass this to a function
                double r
                    = sqrt((prop_e->length() - w.pos.x()) * (prop_e->length() - w.pos.x())
                           + w.pos.y() * w.pos.y());
                for (auto ii : get_intervals_source(Vector2(prop_e->length(), 0.0), w,
                                                    e_face, w.ps_d + r))
                    new_intervals.push_back(ii);
            }

            // generate normal intervals
            // pass this to a function
            for (auto ii : get_intervals_source(w.pos, w, e_face, w.ps_d))
                new_intervals.push_back(ii);
        }
        return new_intervals;
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
                    // check for invert
                    // Interval(double x_, double y_, double st_, double end_, double
                    // ps_d_, Face *from_,
                    // Edge *edge_, bool invert)
                    ////on which side of e does p lie?
                    bool invert = false;
                    if ((e->getEndpoint(0) == source->p
                         and e->getEndpoint(1) < e->getEndpoint(0))
                        or (e->getEndpoint(1) == source->p
                            and e->getEndpoint(0) < e->getEndpoint(1)))
                        invert = true;
                    Interval ii(0, 0, 0, e->length(), 0, NULL, e, invert);
                    insert_new_interval(ii);
                }
                break;

            case Point::EDGE: {
                // source->p is edge *
                Edge *e = (Edge *)source->p;
                bool invert = e->getEndpoint(0) > e->getEndpoint(1);
                // x = (pos - endpoint[0]).length()
                // endpoint[0]---------pos------------endpoint[1]
                Interval ii((source->pos - e->getEndpoint(0)->getPosition()).length(), 0,
                            0, e->length(), 0, NULL, e, invert);
                insert_new_interval(ii);
                break;
            }
            case Point::FACE:
                for (auto &e : ((Face *)(source->p))->edges) {
                    Vector3 pos1 = e->getEndpoint(0)->getPosition();
                    Vector3 pos2 = e->getEndpoint(1)->getPosition();
                    //         pos(y)
                    //       /
                    //     /
                    //(0)pos1--(x)-----------pos2
                    double x = (source->pos - pos1).dot((pos2 - pos1).unit());
                    double y = sqrt((source->pos - pos1).squaredLength() - x * x);
                    // pythagoras theoram
                    // Interval(double x_, double y_, double st_, double end_, double
                    // ps_d_, Face *from_,
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