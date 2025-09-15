#include "atominfo.hpp"
#include <cassert>
#include <string>
#include <iostream>
#include <cstdlib>

#define BUFF_LEN 128

AtomInfo::AtomInfo(std::size_t indx, std::string name, double mass)
:indx(indx), name(name), mass(mass) {}

std::string AtomInfo::to_string() const {
    char buff[BUFF_LEN];
    std::sprintf(buff, "%3ld  %-2s  %12.8lf", 
                        indx, name.data(), mass);
    return std::string(buff);
}

std::size_t  get_atom_indx(std::string name) {
    for (std::size_t i = 1; i <= MAX_ATOM_INDX; ++i) {
        if (elem_table[i].name == name) return i;
    }

    std::cerr << "No such atom: \"" + name + "\"" << std::endl;
    std::abort();

    return 0;
}

std::string  get_atom_name(std::size_t indx){
    if (indx == 0 || indx > MAX_ATOM_INDX) {
        std::cerr << "Invaild atom index: "
                  << indx << std::endl;
        std::abort();
    }

    return elem_table[indx].name;
}

AtomInfo elem_table[MAX_ATOM_INDX + 1] = {
    AtomInfo(  0,  "X" ,  0.0000000000),
    AtomInfo(  1,  "H" ,  1.0080000000),
    AtomInfo(  2,  "He",  4.0026020000),
    AtomInfo(  3,  "Li",  6.9400000000),
    AtomInfo(  4,  "Be",  9.0121831000),
    AtomInfo(  5,  "B" ,  10.810000000),
    AtomInfo(  6,  "C" ,  12.011000000),
    AtomInfo(  7,  "N" ,  14.007000000),
    AtomInfo(  8,  "O" ,  15.999000000),
    AtomInfo(  9,  "F" ,  18.998403163),
    AtomInfo( 10,  "Ne",  20.179700000),
    AtomInfo( 11,  "Na",  22.989769280),
    AtomInfo( 12,  "Mg",  24.305000000),
    AtomInfo( 13,  "Al",  26.981538500),
    AtomInfo( 14,  "Si",  28.085000000),
    AtomInfo( 15,  "P" ,  30.973761998),
    AtomInfo( 16,  "S" ,  32.060000000),
    AtomInfo( 17,  "Cl",  35.450000000),
    AtomInfo( 18,  "Ar",  39.948000000),
    AtomInfo( 19,  "K" ,  39.098300000),
    AtomInfo( 20,  "Ca",  40.078000000),
    AtomInfo( 21,  "Sc",  44.955908000),
    AtomInfo( 22,  "Ti",  47.867000000),
    AtomInfo( 23,  "V" ,  50.941500000),
    AtomInfo( 24,  "Cr",  51.996100000),
    AtomInfo( 25,  "Mn",  54.938044000),
    AtomInfo( 26,  "Fe",  55.845000000),
    AtomInfo( 27,  "Co",  58.933194000),
    AtomInfo( 28,  "Ni",  58.693400000),
    AtomInfo( 29,  "Cu",  63.546000000),
    AtomInfo( 30,  "Zn",  65.380000000),
    AtomInfo( 31,  "Ga",  69.723000000),
    AtomInfo( 32,  "Ge",  72.630000000),
    AtomInfo( 33,  "As",  74.921595000),
    AtomInfo( 34,  "Se",  78.971000000),
    AtomInfo( 35,  "Br",  79.904000000),
    AtomInfo( 36,  "Kr",  83.798000000),
    AtomInfo( 37,  "Rb",  85.467800000),
    AtomInfo( 38,  "Sr",  87.620000000),
    AtomInfo( 39,  "Y" ,  88.905840000),
    AtomInfo( 40,  "Zr",  91.224000000),
    AtomInfo( 41,  "Nb",  92.906370000),
    AtomInfo( 42,  "Mo",  95.950000000),
    AtomInfo( 43,  "Tc",  97.907210000),
    AtomInfo( 44,  "Ru",  101.07000000),
    AtomInfo( 45,  "Rh",  102.90550000),
    AtomInfo( 46,  "Pd",  106.42000000),
    AtomInfo( 47,  "Ag",  107.86820000),
    AtomInfo( 48,  "Cd",  112.41400000),
    AtomInfo( 49,  "In",  114.81800000),
    AtomInfo( 50,  "Sn",  118.71000000),
    AtomInfo( 51,  "Sb",  121.76000000),
    AtomInfo( 52,  "Te",  127.60000000),
    AtomInfo( 53,  "I" ,  126.90447000),
    AtomInfo( 54,  "Xe",  131.29300000),
    AtomInfo( 55,  "Cs",  132.90545196),
    AtomInfo( 56,  "Ba",  137.32700000),
    AtomInfo( 57,  "La",  138.90547000),
    AtomInfo( 58,  "Ce",  140.11600000),
    AtomInfo( 59,  "Pr",  140.90766000),
    AtomInfo( 60,  "Nd",  144.24200000),
    AtomInfo( 61,  "Pm",  144.91276000),
    AtomInfo( 62,  "Sm",  150.36000000),
    AtomInfo( 63,  "Eu",  151.96400000),
    AtomInfo( 64,  "Gd",  157.25000000),
    AtomInfo( 65,  "Tb",  158.92535000),
    AtomInfo( 66,  "Dy",  162.50000000),
    AtomInfo( 67,  "Ho",  164.93033000),
    AtomInfo( 68,  "Er",  167.25900000),
    AtomInfo( 69,  "Tm",  168.93422000),
    AtomInfo( 70,  "Yb",  173.05400000),
    AtomInfo( 71,  "Lu",  174.96680000),
    AtomInfo( 72,  "Hf",  178.49000000),
    AtomInfo( 73,  "Ta",  180.94788000),
    AtomInfo( 74,  "W" ,  183.84000000),
    AtomInfo( 75,  "Re",  186.20700000),
    AtomInfo( 76,  "Os",  190.23000000),
    AtomInfo( 77,  "Ir",  192.21700000),
    AtomInfo( 78,  "Pt",  195.08400000),
    AtomInfo( 79,  "Au",  196.96656900),
    AtomInfo( 80,  "Hg",  200.59200000),
    AtomInfo( 81,  "Tl",  204.38000000),
    AtomInfo( 82,  "Pb",  207.20000000),
    AtomInfo( 83,  "Bi",  208.98040000),
    AtomInfo( 84,  "Po",  208.98243000),
    AtomInfo( 85,  "At",  209.98715000),
    AtomInfo( 86,  "Rn",  222.01758000),
    AtomInfo( 87,  "Fr",  223.01974000),
    AtomInfo( 88,  "Ra",  226.02541000),
    AtomInfo( 89,  "Ac",  227.02775000),
    AtomInfo( 90,  "Th",  232.03770000),
    AtomInfo( 91,  "Pa",  231.03588000),
    AtomInfo( 92,  "U" ,  238.02891000),
    AtomInfo( 93,  "Np",  237.04817000),
    AtomInfo( 94,  "Pu",  244.06421000),
    AtomInfo( 95,  "Am",  243.06138000),
    AtomInfo( 96,  "Cm",  247.07035000),
    AtomInfo( 97,  "Bk",  247.07031000),
    AtomInfo( 98,  "Cf",  251.07959000),
    AtomInfo( 99,  "Es",  252.08300000),
    AtomInfo(100,  "Fm",  257.09511000),
    AtomInfo(101,  "Md",  258.09843000),
    AtomInfo(102,  "No",  259.10100000),
    AtomInfo(103,  "Lr",  262.11000000),
    AtomInfo(104,  "Rf",  267.12200000),
    AtomInfo(105,  "Db",  268.12600000),
    AtomInfo(106,  "Sg",  271.13400000),
    AtomInfo(107,  "Bh",  270.13300000),
    AtomInfo(108,  "Hs",  269.13380000),
    AtomInfo(109,  "Mt",  278.15600000),
    AtomInfo(110,  "Ds",  281.16500000),
    AtomInfo(111,  "Rg",  281.16600000),
    AtomInfo(112,  "Cn",  285.17700000),
    AtomInfo(113,  "Nh",  286.18200000),
    AtomInfo(114,  "Fl",  289.19000000),
    AtomInfo(115,  "Mc",  289.19400000),
    AtomInfo(116,  "Lv",  293.20400000),
    AtomInfo(117,  "Ts",  293.20800000),
    AtomInfo(118,  "Og",  294.21400000),
};


// Non-relativistic spin-restricted spherically averaged Hartree-Fock
// configurations for use in atomic SAD calculations. Reference
// configurations from Phys. Rev. A 101, 012516 (2020).
int nrsrhf_config[MAX_ATOM_INDX + 1][4] = {
    { 0, 0, 0, 0},     //   0  GHOST
    { 1, 0, 0, 0},     //   1  H
    { 2, 0, 0, 0},     //   2  He
    { 3, 0, 0, 0},     //   3  Li
    { 4, 0, 0, 0},     //   4  Be
    { 4, 1, 0, 0},     //   5  B
    { 4, 2, 0, 0},     //   6  C
    { 4, 3, 0, 0},     //   7  N
    { 4, 4, 0, 0},     //   8  O
    { 4, 5, 0, 0},     //   9  F
    { 4, 6, 0, 0},     //  10  Ne
    { 5, 6, 0, 0},     //  11  Na
    { 6, 6, 0, 0},     //  12  Mg
    { 6, 7, 0, 0},     //  13  Al
    { 6, 8, 0, 0},     //  14  Si
    { 6, 9, 0, 0},     //  15  P
    { 6,10, 0, 0},     //  16  S
    { 6,11, 0, 0},     //  17  Cl
    { 6,12, 0, 0},     //  18  Ar
    { 7,12, 0, 0},     //  19  K
    { 8,12, 0, 0},     //  20  Ca
    { 8,13, 0, 0},     //  21  Sc
    { 8,12, 2, 0},     //  22  Ti
    { 8,12, 3, 0},     //  23  V
    { 8,12, 4, 0},     //  24  Cr
    { 6,12, 7, 0},     //  25  Mn
    { 6,12, 8, 0},     //  26  Fe
    { 6,12, 9, 0},     //  27  Co
    { 6,12,10, 0},     //  28  Ni
    { 7,12,10, 0},     //  29  Cu
    { 8,12,10, 0},     //  30  Zn
    { 8,13,10, 0},     //  31  Ga
    { 8,14,10, 0},     //  32  Ge
    { 8,15,10, 0},     //  33  As
    { 8,16,10, 0},     //  34  Se
    { 8,17,10, 0},     //  35  Br
    { 8,18,10, 0},     //  36  Kr
    { 9,18,10, 0},     //  37  Rb
    {10,18,10, 0},     //  38  Sr
    {10,19,10, 0},     //  39  Y
    {10,18,12, 0},     //  40  Zr
    {10,18,13, 0},     //  41  Nb
    { 8,18,16, 0},     //  42  Mo
    { 8,18,17, 0},     //  43  Tc
    { 8,18,18, 0},     //  44  Ru
    { 8,18,19, 0},     //  45  Rh
    { 8,18,20, 0},     //  46  Pd
    { 9,18,20, 0},     //  47  Ag
    {10,18,20, 0},     //  48  Cd
    {10,19,20, 0},     //  49  In
    {10,20,20, 0},     //  50  Sn
    {10,21,20, 0},     //  51  Sb
    {10,22,20, 0},     //  52  Te
    {10,23,20, 0},     //  53  I
    {10,24,20, 0},     //  54  Xe
    {11,24,20, 0},     //  55  Cs
    {12,24,20, 0},     //  56  Ba
    {12,24,21, 0},     //  57  La
    {12,24,22, 0},     //  58  Ce
    {12,24,21, 2},     //  59  Pr
    {12,24,20, 4},     //  60  Nd
    {12,24,20, 5},     //  61  Pm
    {12,24,20, 6},     //  62  Sm
    {12,24,20, 7},     //  63  Eu
    {11,24,20, 9},     //  64  Gd
    {10,24,20,11},     //  65  Tb
    {10,24,20,12},     //  66  Dy
    {10,24,20,13},     //  67  Ho
    {10,24,20,14},     //  68  Er
    {11,24,20,14},     //  69  Tm
    {12,24,20,14},     //  70  Yb
    {12,25,20,14},     //  71  Lu
    {12,24,22,14},     //  72  Hf
    {12,24,23,14},     //  73  Ta
    {10,24,26,14},     //  74  W
    {10,24,27,14},     //  75  Re
    {10,24,28,14},     //  76  Os
    {10,24,29,14},     //  77  Ir
    {10,24,30,14},     //  78  Pt
    {11,24,30,14},     //  79  Au
    {12,24,30,14},     //  80  Hg
    {12,25,30,14},     //  81  Tl
    {12,26,30,14},     //  82  Pb
    {12,27,30,14},     //  83  Bi
    {12,28,30,14},     //  84  Po
    {12,29,30,14},     //  85  At
    {12,30,30,14},     //  86  Rn
    {13,30,30,14},     //  87  Fr
    {14,30,30,14},     //  88  Ra
    {14,30,31,14},     //  89  Ac
    {14,30,32,14},     //  90  Th
    {14,30,30,17},     //  91  Pa
    {14,30,30,18},     //  92  U
    {14,30,30,19},     //  93  Np
    {13,30,30,21},     //  94  Pu
    {12,30,30,23},     //  95  Am
    {12,30,30,24},     //  96  Cm
    {12,30,30,25},     //  97  Bk
    {12,30,30,26},     //  98  Cf
    {12,30,30,27},     //  99  Es
    {12,30,30,28},     // 100  Fm
    {13,30,30,28},     // 101  Md
    {14,30,30,28},     // 102  No
    {14,30,31,28},     // 103  Lr
    {14,30,32,28},     // 104  Rf
    {14,30,33,28},     // 105  Db
    {12,30,36,28},     // 106  Sg
    {12,30,37,28},     // 107  Bh
    {12,30,38,28},     // 108  Hs
    {12,30,39,28},     // 109  Mt
    {12,30,40,28},     // 110  Ds
    {13,30,40,28},     // 111  Rg
    {14,30,40,28},     // 112  Cn
    {14,31,40,28},     // 113  Nh
    {14,32,40,28},     // 114  Fl
    {14,33,40,28},     // 115  Mc
    {14,34,40,28},     // 116  Lv
    {14,35,40,28},     // 117  Ts
    {14,36,40,28},     // 118  Og
};

