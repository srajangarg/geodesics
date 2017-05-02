#include "Algo.h"

vector<Interval> MMP::prop_thru_interval(Vector2 src, Interval &w, Face *face,
                                         double ps_d)
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

        if (info.possible and info.e1 < edge->length()) {

            double x = info.r * cos(info.angle);
            if (abs(x) < EPS)
                x = 0;

            double y = abs(info.r * sin(info.angle));
            if (y < EPS)
                y = 0;

            Interval ii(x, y, info.e1, min(info.e2, edge->length()),
                        ps_d, face, edge, invert);

            set_parent_iterator(ii, w);
            new_intervals.push_back(ii);
        }
        // otherwise doesn't intersect this edge
    }

    return new_intervals;
}

void MMP::set_parent_iterator(Interval & ii, Interval & parent)
{
    //this is called implies that ii is valid
    for (auto it = edge_intervals[parent.edge].begin();it != edge_intervals[parent.edge].end(); it++)
    {
        if (*it == parent)
        {
            ii.parent = it;
            return;
        }
    }
    //shouldn't reach here
    assert(false);
}

void MMP::propagate()
{
    auto top_ite = *intervals_heap.begin();
    Interval prop_w = *top_ite;
    Edge *prop_e = prop_w.edge;

    intervals_heap.erase(intervals_heap.begin());
    update_not_reached(top_ite);

    // should this be done before or after inserting new intervals?
    top_ite->propagated = true;

    for (auto &face : prop_e->faces) {
        if (prop_w.from == face)
            continue;

        // cout << "Propagating " << prop_w << " on " << *face << endl;
        auto candidates = get_new_intervals(prop_w, face);

        // cout << endl;
        for (auto &new_w : candidates)
            insert_new_interval(new_w);
    }
}

vector<Interval> MMP::get_new_intervals(Interval &w, Face *face)
{
    // propogate this interval to edges on `face`
    Edge *edge = w.edge;
    vector<Interval> new_intervals;
    Vertex *v0 = edge->getEndpoint(0), *v1 = edge->getEndpoint(1);

    if (w.st < EPS and v0->saddle_or_boundary) {
        // treat v0 as a new source

        Vector2 new_src = Vector2(0, 0);
        auto intervals
            = prop_thru_interval(new_src, w, face, w.ps_d + (new_src - w.pos).length());

        // if (intervals.size())
        //     cout << intervals.size() << " to be added for saddle " << *v0 << ":"
        //          << endl;

        for (auto &ii : intervals)
            /*cout << ii << endl,*/ new_intervals.push_back(ii);
    }

    if (edge->length() - w.end < EPS and v1->saddle_or_boundary) {
        // treat v1 as a new source

        Vector2 new_src = Vector2(edge->length(), 0);

        auto intervals
            = prop_thru_interval(new_src, w, face, w.ps_d + (new_src - w.pos).length());
        // if (intervals.size())
        //     cout << intervals.size() << " to be added for saddle " << *v1 << ":"
        //          << endl;

        for (auto &ii : intervals)
            /*cout << ii << endl,*/ new_intervals.push_back(ii);
    }

    // generate intervals from (usual) current source
    auto intervals = prop_thru_interval(w.pos, w, face, w.ps_d);

    // if (intervals.size())
    //     cout << intervals.size() << " to be added for usual propogate:" << endl;
    for (auto &ii : intervals)
        /*cout << ii << endl,*/ new_intervals.push_back(ii);
    return new_intervals;
}