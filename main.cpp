#include "experiments/experiment_cli.h"

#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "RESULT status=failed reason=\"missing experiment; expected exp0--exp9\"\n";
        return 2;
    }
    const std::string command = argv[1];

    if (command == "exp0") return tacc::run_exp0(argc, argv);

    if (command == "exp1") return tacc::run_exp1(argc, argv);

    if (command == "exp2") return tacc::run_exp2(argc, argv);

    if (command == "exp3") return tacc::run_exp3(argc, argv);

    if (command == "exp4") return tacc::run_exp4(argc, argv);

    if (command == "exp5") return tacc::run_exp5(argc, argv);

    if (command == "exp6") return tacc::run_exp6(argc, argv);

    if (command == "exp7") return tacc::run_exp7(argc, argv);

    if (command == "exp8") return tacc::run_exp8(argc, argv);

    if (command == "exp9") return tacc::run_exp9(argc, argv);

    std::cerr << "RESULT status=failed reason=\"unknown experiment: " << command << "\"\n";
    return 2;
}
