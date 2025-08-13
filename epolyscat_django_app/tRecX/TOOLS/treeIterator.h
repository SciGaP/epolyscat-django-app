// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef TREE_ITERATOR_H
#define TREE_ITERATOR_H

#include <cassert>
#include <iterator>
#include <limits>

#include "tree.h"

template <typename T>
class tree_iterator { };

// given a level depth, the iterator will run on this level only, skipping subtrees that won't reach this level (relative to the root)
template <typename T>
class strict_level_iterator : public tree_iterator<T> { };

// given a level depth, the iterator will run on the level or at the deepest possible level for a certain subtree (relative to the root)
template <typename T>
class level_iterator : public tree_iterator<T> {
    const Tree<T>* _root;
    int _nominalLevel;
    int _actualLevel;
    const Tree<T>* _curr;
    std::function<bool (const Tree<T>*)> _crit;
public:
    static constexpr int custom_level = std::numeric_limits<int>::max();

    level_iterator(const Tree<T>* root, int level)
        : _root{ root }, _nominalLevel{ level }, _actualLevel{ 0 }, _curr{ root }, _crit{ nullptr }
    {
        _descent();
    }

    level_iterator(const Tree<T>* root, std::function<bool (const Tree<T>*)> criterium)
        : _root{ root }, _nominalLevel{ custom_level }, _actualLevel{ 0 }, _curr{ root }, _crit{ criterium }
    {
        _descent();
    }

    level_iterator(const level_iterator&) = default;
    level_iterator(level_iterator&&) = default;
    ~level_iterator() = default;

    level_iterator()
        : _root{ nullptr }, _nominalLevel{ }, _actualLevel{ }, _curr{ nullptr }
    { }

    // copy-and-swap idiom
    level_iterator<T>& operator =(level_iterator<T> it) {
        swap(it, *this);
        return *this;
    }

    friend void swap(level_iterator<T>& lhs, level_iterator<T>& rhs) {
        using std::swap;
        swap(lhs._root, rhs._root);
        swap(lhs._nominalLevel, rhs._nominalLevel);
        swap(lhs._actualLevel, rhs._actualLevel);
        swap(lhs._curr, rhs._curr);
        swap(lhs._crit, rhs._crit);
    }

private:
    // try to descent as deep as possible from the current node downwards on the left-most branch up to nominal level
    void _descent() {
        assert(_curr != nullptr); // do not advance beyond end-iterator
        while(_actualLevel < _nominalLevel) {
            if(_curr->isLeaf() || (_nominalLevel == custom_level && _crit(_curr))) break;
            _curr = _curr->child(0);
            ++_actualLevel;
        } 
    }
    
    void _advance() {
        assert(_curr != nullptr); // do not advance beyond end-iterator
        if(_curr->parent() == nullptr) {
            _curr = nullptr;
            return; // advanced on the absolute root
        }

        auto as_tree = [](const T* t) {
            return static_cast<const Tree<T>*>(t);
        };

        const Tree<T>* node = _curr->rightSibling();
        if(node != nullptr) {
            _curr = node;
            _descent();
            return;
        }

        --_actualLevel;
        for(node = as_tree(_curr->Tree<T>::parent());
            node->Tree<T>::parent() != nullptr && node->nSibling() == node->Tree<T>::parent()->childSize() - 1;
            node = as_tree(node->Tree<T>::parent()), --_actualLevel)
        {
            if(node == _root) {
                _curr = nullptr;
                return;
            }
            assert(node != nullptr); // node has to be in subtree
        }

        if(node == _root || node->rightSibling() == nullptr) {
            _curr = nullptr;
            return;
        }

        _curr = node->rightSibling();
        _descent();
    }

public:
    template <typename U>
    friend bool operator ==(const level_iterator<U>& lhs, const level_iterator<U>& rhs) {
        return lhs._curr == rhs._curr;
    }

    level_iterator<T>& operator ++() {
        _advance();
        return *this;
    }

    level_iterator<T> operator ++(int) {
        level_iterator<T> ret{ *this };
        _advance();
        return ret;
    }

    const Tree<T>* operator ->() {
        return _curr;
    }

    const Tree<T>* operator *() {
        return _curr;
    }

};

template <typename U>
bool operator !=(const level_iterator<U>& lhs, const level_iterator<U>& rhs) {
    return !(lhs == rhs);
}

namespace std {
    template <typename T>
    struct iterator_traits<level_iterator<T>> {
        using value_type = T*;
        using pointer = value_type*;
        using reference = const value_type&;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag; // can be generalized to bi-directional
    };
}

#endif
