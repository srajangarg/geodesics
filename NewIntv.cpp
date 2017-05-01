#include "Algo.h"

vector<Interval> MMP::source_bisect(double st, double end, const Interval &i1,
                                    const Interval &i2)
{
    // takes in  a segment st-end and returns 1or2 interval
    // with the closest source indicated on each interval
    // cout<<"st : "<<st<<" end : "<<end<<endl;
    // cout<<"i1 : "<<i1<<endl;
    // cout<<"i2 : "<<i2<<endl;
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

    // cout<<"a : "<<a<<" b : "<<b<<" c : "<<c<<endl;

    if ((b * b - 4 * a * c < 0) or (abs(b) < EPS and abs(c) < EPS))
        equidist_pt = false;
    else if (abs(a) < EPS)
        x = -c / b;
    else {
        if (a < 0)
            x = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
        else
            x = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    }

    // cout<<"a "<<a<<" b "<<b<<" c "<<c<<endl;
    // cout<<"x : "<<x<<endl;
    if (equidist_pt and x > st and x < end) {
        // [st, x] closer to ps1 or ps2?
        auto i1val = ((ps1 - st_v).length() + i1.ps_d);
        auto i2val = ((ps2 - st_v).length() + i2.ps_d);

        // cout<<"i1val : "<<i1val<<endl;
        // cout<<"i2val : "<<i2val<<endl;

        if (i1val < i2val) {
            b_intervals.push_back(Interval(ps1, st, x, i1));
            b_intervals.push_back(Interval(ps2, x, end, i2));
        } else {
            b_intervals.push_back(Interval(ps2, st, x, i2));
            b_intervals.push_back(Interval(ps1, x, end, i1));
        }
    } else {
        // the complete interval is near to one single source
        auto i1val = ((ps1 - st_v).length() + i1.ps_d);
        auto i2val = ((ps2 - st_v).length() + i2.ps_d);

        if (abs(i1val - i2val) < EPS) {
            i1val = ((ps1 - end_v).length() + i1.ps_d);
            i2val = ((ps2 - end_v).length() + i2.ps_d);
        }

        if (i2val < i1val)
            b_intervals.push_back(Interval(ps2, st, end, i2));
        else
            b_intervals.push_back(Interval(ps1, st, end, i1));
    }

    // cout<<"final ints : "<<endl;
    // for (auto ii : b_intervals)
    //     cout<<ii<<endl;

    // cout<<"done"<<endl;

    return b_intervals;
}

vector<Interval> MMP::sanitize_and_merge(vector<Interval> &intervals)
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
            if ((last_merged_interval.pos - w.pos).length() < EPS
                and abs(last_merged_interval.ps_d - w.ps_d) < EPS
                and last_merged_interval.from == w.from) {
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

void MMP::insert_new_interval(Interval &new_w)
{
    // updates edge_intervals[e] and intervals according to algo discussed
    // should be O(edge_intervals[e])
    // cout << "Inserting new interval " << new_w << endl;

    auto &intervals = edge_intervals[new_w.edge];
    vector<Interval> new_intervals;
    bool new_fully_pushed = false;

    auto it = intervals.begin();
    for (; it != intervals.end(); it++) {
        // COPY happening here!!!
        auto w = *it;

        if (new_fully_pushed)
            break;

        if (w.st >= new_w.end) {
            // ----new----
            //               -----w-----
            new_intervals.push_back(new_w);
            new_intervals.push_back(w);
            new_fully_pushed = true;
        } else if (w.end <= new_w.st) {
            //             -----new-----
            // ----w----
            new_intervals.push_back(w);
        } else if (w.st < new_w.st and w.end <= new_w.end) {
            //       -------new------
            // ------w--------

            Interval w_short(w);
            w_short.end = new_w.st;

            new_intervals.push_back(w_short);
            for (auto &interval : source_bisect(new_w.st, w.end, w, new_w))
                new_intervals.push_back(interval);

            new_w.st = w.end;
        } else if (w.st >= new_w.st and w.end > new_w.end) {
            // ------new------
            //        -----w--------

            Interval new_w_short(new_w);
            new_w_short.end = w.st;
            new_intervals.push_back(new_w_short);
            for (auto &interval : source_bisect(w.st, new_w.end, new_w, w))
                new_intervals.push_back(interval);
            Interval w_short(w);
            w_short.st = new_w.end;
            new_intervals.push_back(w_short);

            new_fully_pushed = true;
        } else if (w.st >= new_w.st and w.end <= new_w.end) {
            // --------new----------
            //      -----w-----

            Interval new_w_short(new_w);
            new_w_short.end = w.st;
            new_intervals.push_back(new_w_short);

            // cout<<"source bisect- ----"<<endl;
            // cout<<w.st<<" "<<w.end<<endl;
            // cout<<"i1 : "<<new_w<<endl;
            // cout<<"i2 : "<<(w)<<endl;
            for (auto &interval : source_bisect(w.st, w.end, new_w, w))
                new_intervals.push_back(interval);
            // cout<<"----"<<endl;

            new_w.st = w.end;
        } else if (w.st < new_w.st and w.end > new_w.end) {
            //     ----new----
            //  ----------w-------

            Interval w_short1(w);
            w_short1.end = new_w.st;
            new_intervals.push_back(w_short1);
            for (auto &interval : source_bisect(new_w.st, new_w.end, new_w, w))
                new_intervals.push_back(interval);
            Interval w_short2(w);
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
        new_intervals.push_back(*it);

    for (auto it = intervals.begin(); it != intervals.end(); it++)
        intervals_heap.erase(it);

    // set<Interval *> not_to_delete;
    // COPY happening here
    auto old_intervals = intervals;
    intervals.clear();

    auto sanitized_new_intervals = sanitize_and_merge(new_intervals);
    // cout << "Old intervals cleared, all intervals on " << *new_w.edge
    //      << " are :" << endl;

    for (auto &interval : sanitized_new_intervals) {

        Interval *keep_old_itv = NULL;
        for (auto &old_itv : old_intervals)
            if (interval == old_itv)
                keep_old_itv = &old_itv;

        interval.recompute_min_d();
        // cout << interval << endl;

        if (keep_old_itv) {
            intervals.push_back(*keep_old_itv);
            // not_to_delete.insert(keep_old_itv);
            if (not keep_old_itv->propagated) {
                auto pp = intervals_heap.insert(--intervals.end());
                assert(pp.second);
            }
        } else {
            // auto added_intv = new Interval(interval);
            intervals.push_back(interval);
            auto pp = intervals_heap.insert(--intervals.end());
            assert(pp.second);
        }
    }
    // cout << endl;

    // for (auto &itv : old_intervals)
    //     if (not_to_delete.find(itv) == not_to_delete.end())
    //         delete itv;
}