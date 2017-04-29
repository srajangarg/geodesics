#include "Mesh.hpp"
#include "Viewer.hpp"
#include "Algo.h"
#include <cstdlib>

int usage(int argc, char *argv[])
{
    DGP_CONSOLE << "";
    DGP_CONSOLE << "Usage: " << argv[0] << " <mesh-in>";
    DGP_CONSOLE << "";

    return -1;
}

void test(Mesh & mesh)
{   

    for (auto &f : mesh.faces)
    {
        cout<<&f<<endl;
        for (auto v : f.vertices)
            cout<<v<<" ";
        cout<<endl;
    }
    // cout<<endl<<"Done\n";

    for (auto &v : mesh.vertices)
        v.update_saddle_or_boundary();
    // cout<<"oul"<<endl;

    Point p(&(mesh.vertices.front()));

    Point d(&(mesh.vertices.back()));
    
    MMP mmp(&mesh, p, d);
    // std::cout<<"src : edge "<<mesh.edges.front().getEndpoint(0)->getPosition()<<" "
    // <<mesh.edges.front().getEndpoint(1)->getPosition()<<std::endl;

    cout<<"src : vertex "<<mesh.vertices.front().getPosition()<<endl;


    mmp.source = p;

    mmp.algorithm();

    for (auto &pp : mmp.edge_intervals)
    {   
        cout<<endl;
        cout<<*(pp.first)<<" : "<<endl;
        for (auto& w: pp.second)
            cout<<*w<<endl;
    }
    cout<<endl;

    return;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
        return usage(argc, argv);

    std::string in_path = argv[1];

    Mesh mesh;
    if (!mesh.load(in_path))
        return -1;

    DGP_CONSOLE << "Read mesh '" << mesh.getName() << "' with " << mesh.vertices.size()
                << " vertices, " << mesh.edges.size() << " edges and "
                << mesh.faces.size() << " faces from " << in_path;

    test(mesh);
    // Viewer viewer;
    // viewer.setObject(&mesh);
    // viewer.launch(argc, argv);

    return 0;
}
