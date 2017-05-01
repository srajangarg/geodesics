#include "Mesh.hpp"
#include "Viewer.hpp"
#include "Algo.h"
#include "Djk.h"
#include <cstdlib>
#include <chrono>

int usage(int argc, char *argv[])
{
    DGP_CONSOLE << "";
    DGP_CONSOLE << "Usage: " << argv[0] << " <mesh-in>";
    DGP_CONSOLE << "";

    return -1;
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

    Point p(&(mesh.vertices.front()));
    Point d(&(mesh.vertices.back()));

    MMP mmp(&mesh, p, d);
    DJK djk(&mesh, p, d);

    double len1 = 0, len2 = 0;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto path1 = mmp.algorithm();
    auto t2 = std::chrono::high_resolution_clock::now();
    auto tt1 = chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    for (uint i = 1; i < path1.size(); i++)
        len1 += (path1[i-1].pos - path1[i].pos).length();

    t1 = std::chrono::high_resolution_clock::now();
    auto path2 = djk.algorithm();
    t2 = std::chrono::high_resolution_clock::now();
    auto tt2 = chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    for (uint i = 1; i < path2.size(); i++)
        len2 += (path2[i-1].pos - path2[i].pos).length();

    cout<<endl;
    cout<<"Path length MMP : "<<len1<<endl;
    cout<<"Path length DJK : "<<len2<<endl;
    cout<<"Ratio : "<<(len1/len2)<<endl<<endl;

    cout<<endl;
    cout<<"Time taken MMP : "<<tt1<<" us"<<endl;
    cout<<"Time taken DJK : "<<tt2<<" us"<<endl;
    cout<<"Ratio : "<<(tt1/tt2)<<endl<<endl;

    Viewer viewer;
    viewer.setObject(&mesh, path1, path2);
    viewer.launch(argc, argv);

    return 0;
}
