#ifndef EOS_HPP
#define EOS_HPP

#include <string>
/* 
Maybe in the future try to implement EOS as an abstract
class to give flexibility and use cached virtual member 
functions to improve performance
*/

class EOS {
    public:
        EOS(std::string, std::string); // general constructor
    public:
        std::string eos_name;
        std::string type; // may be tabular, polytropic or something special like SIDM EOS

};



#endif