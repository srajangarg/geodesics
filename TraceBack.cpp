#include "Algo.h"

vector<Point> MMP::trace_back(Point destination)
{   
    cout<<"TRACE BACK ------------"<<endl<<endl;

    vector<Point> path;
    // assert(best_interval_dest.find(destination) != best_interval_dest.end());
    // path.push_back(destination);

    // auto cur_itv = best_interval_dest[destination];
    // double cur_x;

    // switch(destination.ptype)
    // {
    //     case Point::VERTEX:
    //         if (cur_itv.edge->getEndpoint(0) == destination.p)
    //             cur_x = 0;
    //         else
    //             cur_x = cur_itv.edge->length();
    //         break;
    //     case Point::EDGE: 
    //         cur_x = cur_itv.edge->length() * destination.ratio;
    //         break;
    //     default:
    //         assert(false);
    // }

    // // cout<<"cur_x : "<<cur_x<<endl;
    // // cout<<"cur_itv : "<<cur_itv<<endl;

    // //saddle and sources not handled
    // int i = 0;
    // while (i < 2)
    // {
    //     cout<<"cur_x : "<<cur_x<<endl;
    //     cout<<"cur_itv : "<<cur_itv<<endl;

    //     for (auto ee : cur_itv.from->edges)
    //     {
    //         //check whether relative position of src wrt this interval gives a vertex or not 
    //         if (ee == cur_itv.edge)
    //             continue;

    //         if ((Vector2(cur_x, 0) - cur_itv.pos).length() > EPS) {
    //             auto cur_e = cur_itv.edge, par_e = ee;
    //             auto x = cur_itv.pos.x(), y = cur_itv.pos.y();
    //             auto e = cur_e->length();
                
    //             auto common = cur_e->getCommonVertex(par_e);
    //             assert(common != NULL);
    //             auto theta = cur_itv.from->getAngle(common);

    //             if (common == cur_e->getEndpoint(0))
    //                 cur_x = e - cur_x, x = e - x;

    //             auto new_x
    //                 = (y * (e - cur_x)) / (sin(theta) * (x - cur_x) + cos(theta) * y);

    //             if (common == par_e->getEndpoint(1))
    //                 new_x = par_e->length() - new_x;

    //             if (new_x >= -EPS and new_x <= par_e->length() + EPS)
    //             {
    //                 new_x = max(0.0, min(par_e->length(), new_x));

    //                 path.push_back(Point(par_e, new_x / par_e->length()));
                    
    //                 for (auto & ii : edge_intervals[par_e])
    //                     if (ii.st < new_x and new_x < ii.end)
    //                         cur_itv = ii;

    //                 cur_x = new_x;
    //                 break;
    //             }
    //             else
    //             {
    //                 //revert the changes
    //                 if (common == cur_e->getEndpoint(0))
    //                     cur_x = e - cur_x;
    //                 continue;
    //             }

    //             //find which interval new_x belongs to

    //         }
    //     }
    //     i++;

    // }

    // while (cur_itv.parent != NULL) {

    //     cout<<"cur_it :"<<*cur_itv<<endl;
    //     cout<<"cur_x  :"<<cur_x<<endl;
    //     cout<<"par    :"<<(cur_itv->parent)<<endl;

    //     if ((Vector2(cur_x, 0) - cur_itv->pos).length() > EPS) {
    //         auto cur_e = cur_itv->edge, par_e = cur_itv->parent->edge;
    //         auto x = cur_itv->pos.x(), y = cur_itv->pos.y();
    //         auto e = cur_e->length();
            
    //         auto common = cur_e->getCommonVertex(par_e);
    //         assert(common != NULL);
    //         auto theta = cur_itv->from->getAngle(common);

    //         if (common == cur_e->getEndpoint(0))
    //             cur_x = e - cur_x, x = e - x;

    //         auto new_x
    //             = (y * (e - cur_x)) / (sin(theta) * (x - cur_x) + cos(theta) * y);

    //         if (common == par_e->getEndpoint(1))
    //             new_x = par_e->length() - new_x;

    //         assert(new_x >= -EPS and new_x <= par_e->length() + EPS);
    //         new_x = max(0.0, min(par_e->length(), new_x));

    //         path.push_back(Point(par_e, new_x / par_e->length()));
    //         cur_x = new_x;
    //         cur_itv = cur_itv->parent;
    //     }
    // }

    return path;
}

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