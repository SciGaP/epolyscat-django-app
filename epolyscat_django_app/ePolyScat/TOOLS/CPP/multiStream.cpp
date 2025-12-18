#include "../multiStream.h"

#include "abort.h"

void MultiStream::add(std::string Name, std::ostream * Os){
    if(_os.count(Name)){
        if(Os!=_os[Name])DEVABORT("cannot add "+Name+" to MulitStream, already exists with different pointer");
    }
    replace(Name,Os);
}

void MultiStream::replace(std::string Name, std::ostream* Os){
    if(Os==&std::cout and Name!="Screen")DEVABORT("std::cout must be named \"Screen\"");
    if(_os.count(Name)){
        std::ofstream*  f=dynamic_cast<std::ofstream*>(_os[Name]);
        if(Os!=_os[Name] and f)f->close();
        else DEVABORT("cannot replace "+Name+" with new output stream");
    }
    _os[Name]=Os;

}

void MultiStream::remove(std::string Name){
    if(_os.count(Name)){
        if(_os[Name]!=&std::cout){
            *_os[Name]<<std::flush;
            delete _os[Name];
        }
        _os.erase(Name);
    }
}
