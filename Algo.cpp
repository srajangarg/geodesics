#include "Algo.h"

Interval::Interval(double x_, double y_, double st_, double end_, double ps_d_,
                   Face *from_, Edge *edge_, Interval *par, bool invert)
{
    edge = edge_;
    from = from_;
    ps_d = ps_d_;
    pos = Vector2(x_, y_);
    parent = par;
    set_st_end_pos(st_, end_, invert);
    assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
    // INVARIANT - st is always close to lower endpoint than higher endpoint pointer
    recompute_min_d();
}

Interval::Interval(Vector2 pos_, double st_, double end_, double ps_d_, Face *from_,
                   Edge *edge_, Interval *par, bool invert)
{
    edge = edge_;
    from = from_;
    pos = pos_;
    parent = par;
    ps_d = ps_d_;
    set_st_end_pos(st_, end_, invert);

    assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
    recompute_min_d();
}

Interval::Interval(Vector2 pos_, double st_, double end_, const Interval &i, bool invert)
{
    edge = i.edge;
    from = i.from;
    ps_d = i.ps_d;
    pos = pos_;
    parent = const_cast<Interval *>(&i);
    set_st_end_pos(st_, end_, invert);

    assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
    recompute_min_d();
}

void Interval::recompute_min_d()
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

double Interval::compute_max_d()
{
    Vector2 st_v(st, 0);
    Vector2 end_v(end, 0);

    if (pos.x() >= (st + end) / 2.0)
        return (pos - st_v).length() + ps_d;
    else
        return (pos - end_v).length() + ps_d;
}

bool Interval::operator<(const Interval &rhs) const
{
    if (min_d != rhs.min_d)
        return min_d < rhs.min_d;
    else if (edge != rhs.edge)
        return edge < rhs.edge;
    else
        return st < rhs.st;
}

void Interval::set_st_end_pos(double st_, double end_, bool invert)
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

ostream& operator<<(ostream& os, const Interval& e)
{
    os<<"I("<<e.st<<", "<<e.end<<", "<<e.pos<<", "<<e.ps_d<<", on "<<*e.edge<<")";
    return os;
}
