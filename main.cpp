#include "Mesh.hpp"
#include "Viewer.hpp"
#include "Algo.h"
#include <cstdlib>

int usage(int argc, char *argv[])
{
    DGP_CONSOLE << "";
    DGP_CONSOLE << "Usage: " << argv[0] << " <mesh-in> [<target-num-faces> [<mesh-out>]]";
    DGP_CONSOLE << "";

    return -1;
}

void check_test()
{
    MMP mmp;

    std::cout << "Checking" << std::endl;

    Vertex *v1 = new Vertex(Vector3(0, 0, 0));
    Vertex *v2 = new Vertex(Vector3(0, 1, 0));

    Edge *e = new Edge(v1, v2);

    Interval ii(0.25, 1, 0, 0.7, 0.5, 0, e, 0);
    Interval ii2(0.75, 1, 0.3, 1, 0, 0, e, 0);

    auto ret = mmp.intervals_heap.insert(ii);

    // // Interval(double x_, double y_, double st_, double end_, double ps_d_, Face
    // *from_,
    // //          Edge *edge_, bool invert = false)

    mmp.edge_intervals[e].push_back(ret.first);

    mmp.insert_new_interval(ii2);

    // for (auto interval : mmp.edge_intervals[e])
    // {
    //     Interval i = *interval;
    //     i.print();
    // }

    return;
}

int main(int argc, char *argv[])
{
    // if (argc < 2)
    //     return usage(argc, argv);

    check_test();

    return 0;

    std::string in_path = argv[1];

    long target_num_faces = -1;
    std::string out_path;
    if (argc >= 3) {
        target_num_faces = std::atoi(argv[2]);

        if (argc >= 4)
            out_path = argv[3];
    }

    Mesh mesh;
    if (!mesh.load(in_path))
        return -1;

    DGP_CONSOLE << "Read mesh '" << mesh.getName() << "' with " << mesh.vertices.size()
                << " vertices, " << mesh.edges.size() << " edges and "
                << mesh.faces.size() << " faces from " << in_path;

    if (target_num_faces >= 0 && (int)mesh.faces.size() > target_num_faces) {
        mesh.updateBounds();
    }

    if (!out_path.empty()) {
        if (!mesh.save(out_path))
            return -1;

        DGP_CONSOLE << "Saved mesh to " << out_path;
    }

    Viewer viewer;
    viewer.setObject(&mesh);
    viewer.launch(argc, argv);

    return 0;
}
