#include <iostream>
#include <map>
#include <string>

#include "libdynamix_input_parser.h"

int main () {
 // make map of strings to bools
 std::map<std::string, bool> outputs;

 std::cout << "\n";
 std::cout << "testing getting outputs to create from file ins/parameters.in\n";

 assignOutputs("ins/outputs.in", outputs);

 // print out outputs to be made
 for (std::map<std::string, bool>::iterator it = outputs.begin(); it != outputs.end(); it++) {
  std::cout << "Output file: " << it->first << "\n";
 }

 return 0;
}
