#include "Algo.h"

void MMP::update_not_reached(list<Interval>::iterator &w)
{
    vector<set<Point>::iterator> to_erase;
    for (auto it = not_reached.begin(); it != not_reached.end(); it++) {
        switch (it->ptype) {
            case Point::VERTEX: {

                auto v = (Vertex *)it->p;
                if (not v->hasIncidentEdge(w->edge))
                    continue;

                if (v == w->edge->getEndpoint(0)) {
                    if (w->st > EPS)
                        continue;
                } else {
                    if (w->edge->length() - w->end > EPS)
                        continue;
                }

                best_interval_dest[*it] = w;
                to_erase.push_back(it);
                break;
            }

            case Point::EDGE: {

                auto e = (Edge *)it->p;
                if (e != w->edge)
                    continue;

                if (w->st / e->length() > it->ratio or w->end / e->length() < it->ratio)
                    continue;

                best_interval_dest[*it] = w;
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

void MMP::initialize()
{
    assert(check_mesh_sanity());
    auto visible = source.get_visible_edges();

    cout << "Making initial intervals ..." << endl;
    switch (source.ptype) {
        case Point::VERTEX: {
            for (auto &e : visible) {
                Interval ii (0, 0, 0, e->length(), 0, NULL, e,
                                       (source.p == e->getEndpoint(1)));
                insert_new_interval(ii);
            }
            break;
        }

        case Point::EDGE: {
            auto e = (Edge *)source.p;
            Interval ii(source.ratio * e->length(), 0, 0, e->length(), 0,
                                   NULL, e, false);
            insert_new_interval(ii);
            break;
        }

        case Point::FACE: {
            for (auto &e : visible) {
                Vector3 pos1 = e->getEndpoint(0)->getPosition();
                Vector3 pos2 = e->getEndpoint(1)->getPosition();
                double x = (source.pos - pos1).dot((pos2 - pos1).unit());
                double y = sqrt((source.pos - pos1).squaredLength() - x * x);

                Interval ii(x, y, 0, e->length(), 0, (Face *)(source.p), e, false);
                insert_new_interval(ii);
            }
            break;
        }

        default:
            assert(false);
    }
}