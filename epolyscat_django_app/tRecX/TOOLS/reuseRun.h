#ifndef REUSERUN_H
#define REUSERUN_H

#include <string>
#include <vector>

class ReadInput;

/// decide whether earlier run with identical input is to be re-used
class ReuseRun
{
    std::vector<std::string> _disregard;
    std::string _reuse;
public:
    /// construct from input that contains re-use directives (if any)
    ReuseRun(ReadInput& Inp, std::vector<std::string> Disregard={});
    /// return name of file with matching InputLines ("" if no match)
    std::string match(const std::vector<std::string> & InputLines) const;
};

#endif // REUSERUN_H
