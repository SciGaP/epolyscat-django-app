#ifndef MULTISTREAM_H
#define MULTISTREAM_H

#include <map>
#include <iostream>
#include <fstream>
#include "str.h"

/// a simple multi-stream with limited functionality for small amounts of output
class MultiStream {
    std::map<std::string,std::ostream*> _os;
public:
    MultiStream():_os({{"Screen",&std::cout}}){}
    /// add additional stream Os, unique name
    void add(std::string Name, std::ostream * Os);
    /// replace stream Name with Os
    void replace(std::string Name, std::ostream * Os);
    /// remove stream Name
    void remove(std::string Name);

    /// insert data into all current streems
    template<class T>
    MultiStream& operator<<(const T& x) {
        for(auto o: _os)*(o.second)<<x<<std::flush;
        return *this;
    }
    /// return all current streams
    std::map<std::string,std::ostream*> Streams(){return _os;}
};



#endif // MULTISTREAM_H
