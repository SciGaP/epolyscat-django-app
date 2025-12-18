// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef CIRCULARQUEUE_H
#define CIRCULARQUEUE_H

#include <vector>

/*!
 * \brief Template to handle buffered data is a circular queue container.
 * Current index data member points to the oldest element or the beginning of the queue.
 * It has-a vector, provides operator [] and size(), front(), back() and push_back() member functions
 * also admits operator << to fill the buffer.
 * you should have in mind that back() is the element which the internal counter points to,
 * front is the one it pointed to before.
 *
 * Circular Queue features are emulated using vector instead of using a deque because of contiguous memory in the case of vector.
 */
template<typename bufferedType>
class CircularQueue
{
public:
    explicit CircularQueue(int bufferSize) : current_index(0), buffer(std::vector<bufferedType>(bufferSize)) {
        if (bufferSize==0) {
            buffer.resize(1);
        }
    }
    ~CircularQueue() {} // Attention!!! buffer of pointers will lead to memory leaks!

    bufferedType& operator[] (int i) { ///< allow access via arbitrary i's
        while (i<0) {i+=buffer.size();} // When trying to access negative indices (behaviour of % and / is undefined in this case)
        return buffer.at(i%buffer.size());
    }

    CircularQueue<bufferedType>& operator<< (const bufferedType& rhs) { ///< add a new value and purge the oldest one that was added via <<
        buffer.at(current_index) = rhs;
        ++current_index%=buffer.size();
        return *this;
    }

    bufferedType& back() { ///< reference to oldest element (next to be purged)
        return buffer.at(current_index);
    }

    bufferedType& pop_back() { ///< reference to latest element and Incrementation of the internal counter
        ++current_index%=buffer.size();
        return this->operator [](current_index-1);
    }

    bufferedType& front() { ///< reference to newest element
        return this->operator [](current_index-1);
    }

    int size() const {return buffer.size();} ///< size of buffer

    // data
private:
    int current_index; // DEMAND: 0<=current_index<buffer.size(), current index always points to the oldest element
    std::vector<bufferedType> buffer; // access data only by operator[] !
};

#endif // CIRCULARQUEUE_H
