#include "Algo.h"

// vector<Point> MMP::trace_back(Point destination)
// {   
//     cout<<"TRACE BACK ------------"<<endl<<endl;

//     vector<Point> path;
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

//     return path;
// }


void MMP::best_first_saddle(Vertex* v, double & cur_x, Interval & cur_itv)
{   
    // Vertex * v = e->getEndpoint(endpoint);

    // double distance = std::numeric_limits<double>::infinity();
    // Interval min_interval; 
    // for (auto it = v->edges.begin(); it != v->edges.end(); ++it)
    // {
    //     Interval interval;
    //     double source_dist = interval.ps_d;
    //     if (v == (*it)->getEndpoint(0))
    //     {
    //         interval = *(edge_intervals[*it].front());
    //         if (abs(interval.st) < EPS)
    //         source_dist += (interval.pos - Vector2(interval.st, 0.0)).length();
    //         source_dist += interval.st;
    //     }
    //     else
    //     {
    //         interval = *(edge_intervals[*it].back());
    //         source_dist += (interval.pos - Vector2(interval.end, 0.0)).length();
    //         source_dist += ((*it)->length() - interval.end);
    //     }
    //     if (source_dist < distance)
    //         min_interval = interval;
    // }

    // auto v = ed->getEndpoint(endpoint);
    double best_x;
    double mind = std::numeric_limits<double>::infinity();
    cout<<"bfsd"<<endl;

    Interval bf;
    for (auto &e : v->edges) {
        Interval ii;
        double temp_x;
        double dis;

        if(e == cur_itv.edge)
            continue;

        if (edge_intervals[e].empty())
            continue;

        if (v == e->getEndpoint(0)) {
            ii = edge_intervals[e].front();
            if (ii.st > EPS)
                continue;
            dis = ii.pos.length() + ii.ps_d;
            temp_x = 0.0;
        } else {
            ii = edge_intervals[e].back();
            cout<<ii<<endl;
            if (e->length() - ii.end > EPS)
                continue;
            dis = (ii.pos - Vector2(e->length(), 0)).length() + ii.ps_d;
            temp_x = e->length();
        }

        if (dis < mind) {
            bf = ii;
            best_x = temp_x;
            mind = dis;
        }
    }

    assert(mind != std::numeric_limits<double>::infinity());
    cur_x = best_x;
    cur_itv = bf;
}

vector<Point> MMP::trace_back(Point destination)
{   
    cout<<"TRACE BACK ------------"<<endl<<endl;

    vector<Point> path;
    assert(best_interval_dest.find(destination) != best_interval_dest.end());
    path.push_back(destination);

    auto cur_itv = *best_interval_dest[destination];
    double cur_x;

    switch(destination.ptype)
    {
        case Point::VERTEX:
            if (cur_itv.edge->getEndpoint(0) == destination.p)
                cur_x = 0;
            else
                cur_x = cur_itv.edge->length();
            break;
        case Point::EDGE: 
            cur_x = cur_itv.edge->length() * destination.ratio;
            break;
        default:
            assert(false);
    }

    cout<<"initial : "<<endl;
    cout<<"cur_x : "<<cur_x<<endl;
    cout<<"cur_itv : "<<cur_itv<<endl;

    //saddle and sources not handled
    // int i = 0;
    while (true)
    {

        //check if cur_x and cur_itv 
        cout<<"cur_x : "<<cur_x<<endl;
        cout<<"cur_itv : "<<cur_itv<<endl;

        //check whether cur_x and cur_itv are close to source or saddle
        if ((abs(cur_x) < EPS and source.p == cur_itv.edge->getEndpoint(0)) or 
            (abs(cur_x - cur_itv.edge->length()) < EPS and source.p == cur_itv.edge->getEndpoint(1)) or 
            (source.p == cur_itv.edge and abs(cur_x/cur_itv.edge->length() - source.ratio) < EPS))
        {
            //source reached
            path[path.size() - 1] = source;
            return path;
        }
        if((Vector2(cur_x, 0) - cur_itv.pos).length() < EPS and abs(cur_x) < EPS)
        {
            //get the closest interval to this point
            cout<<"saddle "<<*cur_itv.edge->getEndpoint(0)<<endl;
            cout<<"finding"<<endl;
            best_first_saddle(cur_itv.edge->getEndpoint(0), cur_x, cur_itv);
            cout<<"after saddle 1"<<endl;
            cout<<"cur_sadd : "<<cur_x<<endl;
            cout<<"cur_itv_sad : "<<cur_itv<<endl;
            continue;
        }
        else if((Vector2(cur_x, 0) - cur_itv.pos).length() < EPS and abs(cur_x - cur_itv.edge->length()) < EPS)
        {
            cout<<"saddle "<<*cur_itv.edge->getEndpoint(1)<<endl;
            //get the closest interval from this point
            cout<<"finding"<<endl;
            best_first_saddle(cur_itv.edge->getEndpoint(1), cur_x, cur_itv);
            cout<<"after saddle 2"<<endl;
            cout<<"cur_sadd : "<<cur_x<<endl;
            cout<<"cur_itv_sad : "<<cur_itv<<endl;
            continue;
        }



        //figure out which edge to propogate
        double angle1 = atan2(cur_itv.pos.y(), cur_itv.pos.x() - cur_x);
        double angle2;

        cout<<"a1 : "<<angle1<<endl;
        cout<<"a2 : "<<angle2<<endl;
        int endpoint;
        //next interval propogates through this
        //calculate angle2 = angle 
        for (auto ee : cur_itv.from->edges)
        {
            if (ee == cur_itv.edge)
                continue;

            if (cur_itv.edge->getCommonVertex(ee) == cur_itv.edge->getEndpoint(0))
            {
                double theta = cur_itv.from->getAngle(cur_itv.edge->getEndpoint(0));
                double x = ee->length()*cos(theta);
                double y = ee->length()*sin(theta);
                angle2 = atan2(y, x - cur_x);
            }
        }

        if (angle1 > angle2)
        {
            //propogate it on endpoint-0
            endpoint = 0;
        }
        else
            endpoint = 1;

        cout<<"endp : "<<endpoint<<endl;

        for (auto ee : cur_itv.from->edges)
        {
            //check whether relative position of src wrt this interval gives a vertex or not 
            // cout<<"Iterating on : "<<*ee<<endl;
            // cout<<"cur_x : "<<cur_x<<endl;
            // cout<<"cur_itv : "<<cur_itv<<endl;
            // cout<<ee<<" "<<cur_itv.edge<<endl;
            if (ee == cur_itv.edge)
                continue;

            if (cur_itv.edge->getCommonVertex(ee) != cur_itv.edge->getEndpoint(endpoint))
                continue;

            if ((Vector2(cur_x, 0) - cur_itv.pos).length() > EPS) {
                auto cur_e = cur_itv.edge, par_e = ee;
                auto x = cur_itv.pos.x(), y = cur_itv.pos.y();
                auto e = cur_e->length();
                
                auto common = cur_e->getCommonVertex(par_e);
                assert(common != NULL);
                auto theta = cur_itv.from->getAngle(common);
                cout<<"before : "<<endl;
                cout<<"cur_x : "<<cur_x<<" y : "<<y<<" e : "<<e<<endl;
                if (common == cur_e->getEndpoint(0))
                    cur_x = e - cur_x, x = e - x;

                cout<<"form : after : "<<endl;
                cout<<"cur_x : "<<cur_x<<" y : "<<y<<" e : "<<e<<endl;
                auto new_x
                    = (y * (e - cur_x)) / (sin(theta) * (x - cur_x) + cos(theta) * y);

                if (common == par_e->getEndpoint(1))
                    new_x = par_e->length() - new_x;

                cout<<"edge : "<<*ee<<endl;
                cout<<"cur_x : "<<cur_x<<" new_x : "<<new_x<<endl;

                if (new_x >= -EPS and new_x <= par_e->length() + EPS)
                {
                    cout<<"inside, new_x : "<<new_x<<endl;
                    new_x = max(0.0, min(par_e->length(), new_x));

                    path.push_back(Point(par_e, new_x / par_e->length()));
                    cout<<"par_e : "<<*par_e<<endl;
                    if(abs(new_x) < EPS)
                        cur_itv = edge_intervals[par_e].front();
                    else if(abs(new_x - par_e->length()) < EPS)
                        cur_itv = edge_intervals[par_e].back();
                    else
                    {
                        for (auto & ii : edge_intervals[par_e])
                            if (ii.st <= new_x and new_x <= ii.end/* and ii.from != cur_itv.from*/)
                            {
                                cout<<"new_itv found"<<endl;
                                cur_itv = ii;
                            }
                    }

                    cout<<"new_itv : "<<cur_itv<<endl;

                    cur_x = new_x;

                    // cout<<"cur_x : "<<cur_x<<endl;
                    // cout<<"cur_itv : "<<cur_itv<<endl;
                    break;
                }
                //find which interval new_x belongs to 

            }
        }
        // i++;
    }

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