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
    recompute_min_d();
}

Interval::Interval(Vector2 pos_, double st_, double end_, const Interval &i, bool invert)
{
    edge = i.edge;
    from = i.from;
    ps_d = i.ps_d;
    pos = pos_;
    parent = const_cast<Interval *>(i.parent);
    set_st_end_pos(st_, end_, invert);
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

bool Interval::operator==(const Interval &rhs) const
{
    return (edge == rhs.edge and abs(st - rhs.st) < EPS and abs(end - rhs.end) < EPS
            and abs(ps_d - rhs.ps_d) < EPS and (pos - rhs.pos).length() < EPS);
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

    assert(st < end and st >= 0 and pos.y() >= 0 and end <= edge->length());
}

ostream &operator<<(ostream &os, const Interval &e)
{
    os << "I(st/end: (" << e.st << ", " << e.end << "), pos: " << e.pos
       << ", psd: " << e.ps_d << ", " << *e.edge;

    if (e.from)
        os << ", " << *e.from << ")";
    else
        os << ", NO FACE)";
    return os;
}
