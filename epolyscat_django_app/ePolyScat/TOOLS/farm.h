#ifndef FARM_H
#define FARM_H

#include <memory>
#include <vector>
#include <map>
#include <set>
#include "mpiWrapper.h"
#include "timer.h"

class InputGenerator;
class Farm
{

    MPI_Comm _comFarm;  // overall farm communicator
    MPI_Comm _comFlock; // flock communicator of present thread
    std::vector<int>_flock; // ranks in _comFarm that belong to current flock
    int _maxFlockSize=INT_MAX;

    std::shared_ptr<InputGenerator> _inputGenerator;
    mutable std::string _inp; // path to current inpc-file

    std::vector<std::string> assigned;
    Timer* _statusTimer;
public:
    Farm();

    Farm(std::shared_ptr<InputGenerator> Generator);

    /// generate from command line input
    Farm(int argc, char* argv[]);

    /// advance Farm to next input, false if none
    bool next();

    /// true if current input has completed
    bool runCompleted() const;

    /// return farm communicator
    MPI_Comm communicator() const {return _comFarm;}

    /// current input
    std::string inp();

    /// fork into the flocks, set communicator to this thread's flock
    void fork();

    /// join all flocks, return to overall farm communicator
    void join(){MPIwrapper::Barrier(_comFarm);MPIwrapper::setCommunicator(_comFarm);}

    /// return status info for runs in Assigne
    std::string status();
};

#endif // FARM_H
