#include "Mesh.hpp"
#include "Viewer.hpp"
#include "Algo.h"
#include "Djk.h"
#include <cstdlib>

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

    DJK mmp(&mesh, p, d);
    cout << "src : vertex " << mesh.vertices.front().index << endl;
    cout << "dest : vertex " << mesh.vertices.back().index << endl;
    auto path = mmp.algorithm();

    cout<<endl;
    // for (auto &v : mesh.vertices)
    //     cout<<v<<" "<<v.dis<<endl;
    if (path.size() > 0)
    {
        cout<<path[0];
        for (int i = 1; i < path.size(); i++)
            cout<<" -- "<<path[i];
        cout<<endl;
    }    
    cout << endl;
    
    Viewer viewer;
    viewer.setObject(&mesh, path);
    viewer.launch(argc, argv);

    return 0;
}
