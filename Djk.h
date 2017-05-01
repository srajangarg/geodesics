#include "Mesh.hpp"
#include "Point.hpp"
#include "DGP/Vector2.hpp"
#include <unordered_map>
#include <cmath>
#define EPS 1e-6

class DJK
{
public:
    Mesh *mesh;
    Point source;
    set<Point> not_reached;

    struct IntervalPtrComp {
        bool operator()(const Vertex *lhs, const Vertex *rhs) const
        {
            if (lhs->dis != rhs->dis)
                return lhs->dis < rhs->dis;
            return lhs < rhs;
        }
    };

    multiset<Vertex *, IntervalPtrComp> heap;
    vector<Point> destinations;

    DJK()
    {
    }

    DJK(Mesh *m, Point s, Point d) : mesh(m), source(s)
    {
        assert(s.ptype != Point::UNDEFINED);
        not_reached.insert(d);
        destinations.push_back(d);
    }

    DJK(Mesh *m, Point s, const vector<Point> &dests) : mesh(m), source(s)
    {
        assert(s.ptype != Point::UNDEFINED);
        not_reached.insert(dests.begin(), dests.end());
        destinations = dests;
    }

    void update_not_reached(Vertex *v)
    {
        vector<set<Point>::iterator> to_erase;
        for (auto it = not_reached.begin(); it != not_reached.end(); it++) {
            switch (it->ptype) {
                case Point::VERTEX: {

                    if ((Vertex *)it->p == v)
                        to_erase.push_back(it);
                    break;
                }

                case Point::EDGE: {

                    auto e = (Edge *)it->p;
                    if ((e->getEndpoint(0)->visited and e->getEndpoint(1) == v)
                        or (e->getEndpoint(1)->visited and e->getEndpoint(0) == v))
                        to_erase.push_back(it);
                    break;
                }

                default:
                    assert(false);
            }
        }

        for (auto &it : to_erase)
            not_reached.erase(it);
    }

    void propagate()
    {
        auto prop_v = *heap.begin();

        prop_v->visited = true;
        heap.erase(heap.begin());
        update_not_reached(prop_v);

        for (auto &e : prop_v->edges) {
            auto v = e->getOtherEndpoint(prop_v);
            if (v->visited)
                continue;

            if (prop_v->dis + e->length() < v->dis) {
                auto it = heap.find(v);
                if (it != heap.end())
                    heap.erase(it);
                v->dis = prop_v->dis + e->length();
                v->par = prop_v;
                heap.insert(v);
            }
        }
    }

    vector<Point> trace_back(Point destination)
    {
        vector<Point> path;
        path.push_back(destination);

        Vertex *stv;
        switch (destination.ptype) {
            case Point::VERTEX: {
                stv = (Vertex *)(destination.p);
                break;
            }

            case Point::EDGE: {
                auto e = (Edge *)destination.p;
                auto v0 = e->getEndpoint(0), v1 = e->getEndpoint(1);

                if (v0->dis + (destination.pos - v0->getPosition()).length()
                    < v1->dis + (destination.pos - v1->getPosition()).length())
                    stv = v0;
                else
                    stv = v1;

                path.push_back(Point(stv));
                break;
            }

            default:
                assert(false);
        }

        while (stv->par != NULL) {
            stv = stv->par;
            path.push_back(Point(stv));
        }

        return path;
    }

    bool check_mesh_sanity()
    {
        for (auto &edge : mesh->edges)
            if (edge.faces.size() == 0 or edge.faces.size() > 3)
                return false;
        return true;
    };

    void initialize()
    {
        assert(check_mesh_sanity());

        for (auto &v : mesh->vertices)
            v.visited = false, v.par = NULL, v.dis = std::numeric_limits<double>::max();

        switch (source.ptype) {
            case Point::VERTEX: {

                auto v = (Vertex *)(source.p);
                v->dis = 0;
                heap.insert(v);
                break;
            }

            case Point::EDGE: {
                auto e = (Edge *)source.p;

                for (auto &v : e->endpoints) {
                    v->dis = (v->getPosition() - source.pos).length();
                    heap.insert(v);
                }
                break;
            }

            default:
                assert(false);
        }
    }

    vector<Point> algorithm()
    {
        initialize();

        while (not heap.empty() and not not_reached.empty()) {
            propagate();
        }

        return trace_back(destinations.front());
    }

    // invariants
    // vector of intervals for each edge must be non overlapping and sorted by st
    // all valid windows MUST be in the set
};