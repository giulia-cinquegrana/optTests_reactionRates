#include <mpi.h>
#include <iostream>
#include <cmath>

const double LN2 = std::log(2.0);

enum Preset : int { SNeII = 0, UNKNOWN = -1};

Preset parse_preset(const std::string& s){
    if (s == "SNeII") return SNeII;
    return UNKNOWN;
}

double calc_decay(Preset preset){

    double hl;

    switch (preset) {
        case SNeII:
            std::cout << "following SNeII chains (Ni56 -> Co56 -> Fe56 and Ni57 -> Co57 -> Fe57)";
            hl = 10.0;
            break;
        default:
            throw std::runtime_error("only SNeII available now, come back later");
    }

    return LN2 / hl;
}

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int preset_id = UNKNOWN; 

    if (rank == 0) {
        std::cout << "which chains? (SNeII, X or X)";
        std::string input;
        std::cin >> input; 
        preset_id = static_cast<int>(parse_preset(input));
    }

    MPI_Bcast(&preset_id, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    try {
        double decay_const = calc_decay(static_cast<Preset>(preset_id));
        std::cout << "Rank " << rank
                  << " — decay constant: " << decay_const << " s^-1\n";
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();
    return 0;
}