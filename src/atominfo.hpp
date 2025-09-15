#pragma once

#include <string>
#include <cstddef>

class AtomInfo
{
public:
    std::size_t     indx;
    std::string     name;
    double          mass;

    AtomInfo(std::size_t indx, 
             std::string name, double mass);

    std::string to_string() const;
};

const std::size_t MAX_ATOM_INDX = 118;
extern AtomInfo elem_table[MAX_ATOM_INDX + 1];

std::size_t get_atom_indx(std::string name);
std::string get_atom_name(std::size_t indx);
