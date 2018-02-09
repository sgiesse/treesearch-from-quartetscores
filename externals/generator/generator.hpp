#ifndef GENERATOR_HPP
#define GENERATOR_HPP

// generator/continuation for C++
// author: Andrew Fedoniouk @ terrainformatica.com
// idea borrowed from: "coroutines in C" Simon Tatham,
//                     http://www.chiark.greenend.org.uk/~sgtatham/coroutines.html
// modified by: Sebastian Gieße

//++ coroutine, generator, continuation for C++

struct _generator
{
    int _line;
    _generator():_line(-1) {}
};

#define GENERATOR(NAME) struct NAME : public _generator

#define EMIT(T) bool operator()(T& _rv) {       \
    if(_line < 0) { _line=0;}                   \
        switch(_line) { case 0:;

#define STOP  } _line = 0; return false; }

#define YIELD(V)                                \
    do {                                        \
        _line=__LINE__;                         \
        _rv = (V); return true; case __LINE__:; \
    } while (0)

//-- coroutine, generator, continuation for C++

#endif
