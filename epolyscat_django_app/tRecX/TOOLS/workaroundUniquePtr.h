// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef WORKAROUNDUNIQUEPTR_H
#define WORKAROUNDUNIQUEPTR_H

template<class T>
class WorkaroundUniquePtr{
    // GCC 4.8.5 is inconsistent as shiped with SuSe Leap 42.3
    // this can be used to replace std::unique_ptr
    // Warning: does not enforce uniqueness
    T* _ptr;
public:
    WorkaroundUniquePtr(){}
    WorkaroundUniquePtr(T* Ptr):_ptr(Ptr){}
    T & operator*(){return *_ptr;}
    void operator=(T* Ptr){_ptr=Ptr;}
    T* operator->(){return _ptr;}
    ~WorkaroundUniquePtr(){delete _ptr;}
};


#endif // WORKAROUNDUNIQUEPTR_H
