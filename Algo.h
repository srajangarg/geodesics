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
    bool propagated = false;

    Interval()
    {
    }

    void set_st_end_pos(double st_, double end_, bool invert);

    Interval(double x_, double y_, double st_, double end_, double ps_d_, Face *from_,
             Edge *edge_, bool invert = false);

    Interval(Vector2 pos_, double st_, double end_, double ps_d_, Face *from_,
             Edge *edge_, bool invert = false);

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
        bool possible = true;
    };

    Info get_info_edge(Vector2 src, Edge *e, Face *face) const;

    // y >= 0
    // recompute_min_d() must be called if any of the 5 parameters change
    friend ostream &operator<<(ostream &os, const Interval &e);
};

class MMP
{
public:
    Mesh *mesh;
    unordered_map<Edge *, list<Interval>> edge_intervals;
    Point source;
    map<Point, list<Interval>::iterator> best_interval_dest;
    set<Point> not_reached;

    struct IntervalPtrComp {
        bool operator()(const list<Interval>::iterator lhs, const list<Interval>::iterator rhs) const
        {
            return *lhs < *rhs;
        }
    };
    set<list<Interval>::iterator, IntervalPtrComp> intervals_heap;

    MMP()
    {
    }

    MMP(Mesh *m, Point s, Point d) : mesh(m), source(s)
    {
        assert(s.ptype != Point::UNDEFINED);
        not_reached.insert(d);
    }

    MMP(Mesh *m, Point s, const vector<Point> &dests) : mesh(m), source(s)
    {
        assert(s.ptype != Point::UNDEFINED);
        not_reached.insert(dests.begin(), dests.end());
    }

    vector<Interval> prop_thru_interval(Vector2 src, Interval &w, Face *face, double ps_d);

    void update_not_reached(list<Interval>::iterator &w);

    void best_first_saddle(Edge* e, double & cur_x, Interval & cur_itv, int endpoint);

    void propagate();

    void best_first_saddle(Vertex* v, double & cur_x, Interval & cur_itv);

    vector<Point> trace_back(Point destination);

    vector<Interval> source_bisect(double st, double end, const Interval &i1,
                                   const Interval &i2);

    vector<Interval> sanitize_and_merge(vector<Interval> &intervals);

    void insert_new_interval(Interval &new_w);

    vector<Interval> get_new_intervals(Interval &w, Face *face);

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

    void initialize();

    vector<Point> algorithm()
    {
        initialize();
        // int x;
        // cout << "HEAP ---" << endl;
        // for (auto &itv : intervals_heap)
        //     cout << *itv << endl;
        // cout << "----" << endl << endl;

        // cout << "EDGE MAP ---" << endl;
        // for (auto &pp : edge_intervals) {
        //     cout << endl;
        //     cout << *(pp.first) << " : " << endl;
        //     for (auto &w : pp.second)
        //         cout << w << endl;
        // }
        // cout << "----" << endl << endl;
        // cin>>x;

        while (not intervals_heap.empty() and not not_reached.empty()) {
            propagate();

            // cout << "HEAP ---" << endl;
            // for (auto &itv : intervals_heap)
            //     cout << *itv << endl;
            // cout << "----" << endl << endl;

            cout << "EDGE MAP ---" << endl;
            for (auto &pp : edge_intervals) {
                if ((not (pp.first->getEndpoint(0)->index == 33 and pp.first->getEndpoint(1)->index == 52))
                    and (not (pp.first->getEndpoint(0)->index == 33 and pp.first->getEndpoint(1)->index == 42)))
                    continue;

                // if (not (pp.first->getEndpoint(0)->index == 42 and pp.first->getEndpoint(1)->index == 52))
                //     continue;

                cout << endl;
                cout << *(pp.first) << " : " << endl;
                for (auto &w : pp.second)
                    cout << w << endl;
            }
            cout << "----" << endl << endl;
            // cin>>x;
        }

        return trace_back(best_interval_dest.begin()->first);
    }

    // invariants
    // vector of intervals for each edge must be non overlapping and sorted by st
    // all valid windows MUST be in the set
};