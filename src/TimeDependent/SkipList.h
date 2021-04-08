//
// Created by Pascal Baehr on 07.04.21.
//

#ifndef MOLFLOW_PROJ_SKIPLIST_H
#define MOLFLOW_PROJ_SKIPLIST_H

#include "../TimeMoments.h"

using Skip_t = Moment;

//==============================================================================
struct Skip_Node {
    int key;
    Skip_t value;

    // pointers to successor nodes
    std::vector<Skip_Node*> forward;

    Skip_Node (int k, const Skip_t& v, int level);
};

//==============================================================================
class Skip_list {
public:
    Skip_list ();
    ~Skip_list ();

    // non-modifying member functions
    void print ();
    Skip_Node* find (int searchKey);

    // modifying member functions
    void insert (int searchKey, Skip_t newValue);
    void erase (int searchKey);
private:
    // pointer to first node
    Skip_Node* head;
    // last node
    Skip_Node* NIL;

    // implicitly used member functions
    int randomLevel ();
    int nodeLevel(const std::vector<Skip_Node*>& v);
    Skip_Node* makeNode (int key, Skip_t val, int level);

    // data members
    float probability;
    int maxLevel;
};

#endif //MOLFLOW_PROJ_SKIPLIST_H
