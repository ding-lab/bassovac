#include "BassovacApp.hpp"

#include "bvprob/Bassovac.hpp"
#include <iostream>
#include <stdexcept>

using namespace std;

int main(int argc, char** argv) {
    try {
        BassovacApp app(argc, argv);
        app.run();
    } catch (const exception& e) {
        cerr << e.what() << endl;
        return 1;
    }

    return 0;
}
