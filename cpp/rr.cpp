#include <mpi.h>
#include <iostream>
#include <cmath>
#include <array>
#include <stdexcept>

const double LN2 = std::log(2.0);

enum Preset : int { SNeII = 0, UNKNOWN = -1};

Preset parse_preset(const std::string& s){
    if (s == "SNeII") return SNeII;
    return UNKNOWN;
}

std::array<std::array<double, 2>, 2> calc_decay(Preset preset){

    std::array<std::array<double, 2>, 2> hl;
    std::array<std::array<double, 3>, 2> Y;
    std::array<std::array<double, 3>, 2> dY;
    std::array<std::array<double, 2>, 2> decay_const;

    switch (preset) {
        case SNeII:
            // std::cout << "following SNeII chains (Ni56 -> Co56 -> Fe56 and Ni57 -> Co57 -> Fe57)\n";
            hl = {{
                {{5.251e5, 6.672e6}}, // Ni56 -> Co56 -> Fe56
                {{1.28e5, 2.35e7}} // Ni57 -> Co57 -> Fe57
            }};
            Y = {{
                {{1,0,0}},
                {{0.3,0,0}},
            }};
            dY = {{
                {{0,0,0}},
                {{0,0,0}},
            }};
            break;

        default:
            throw std::runtime_error("only SNeII available now, come back later");
    }

    for (size_t i = 0; i < 2; ++i){
        for (size_t j = 0; j < 2; ++j){

            decay_const[i][j] = LN2 / hl[i][j];
        }
    }

    // for (double val : decay_const[0]) {
    //     std::cout << val << " ";
    // }
    // std::cout << "\n";
    // return decay_const;

    for (std::size_t j = 0; j < Y[0].size(); ++j) {

        Y[0][j] * decay_const[0][j] ;

    }

    return decay_const;
}

            // reaction flux

        // rf1 = decay_const[:,0]* self.Y[:,0]   
        // rf2 = decay_const[:,1]* self.Y[:,1]

        // dY = np.zeros_like(self.Y)

        // dY[:,0] = -rf1
        // dY[:,1] = rf1 - rf2
        // dY[:,2] = rf2

        // return dY

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int preset_id = UNKNOWN; 

    if (rank == 0) {
        std::cout << "which chains? SNeII, X or X\n";
        std::string input;
        std::cin >> input; 
        preset_id = static_cast<int>(parse_preset(input));
    }

    // preset_id is an int
    MPI_Bcast(&preset_id, 1, MPI_INT, 0, MPI_COMM_WORLD);

    try {

        auto decay_consts = calc_decay(static_cast<Preset>(preset_id));

        for (size_t i = 0; i < decay_consts.size(); ++i) {
            for (size_t j = 0; j < decay_consts[i].size(); ++j) {
                std::cout << "Rank " << rank
                          << " — decay_const[" << i << "][" << j << "] = "
                          << decay_consts[i][j] << " s^-1\n";
            }
        }
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();
    return 0;
}