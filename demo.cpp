#include <string>
#include <sstream>
#include <iostream>
#include "serializeable_counting_bloom_filter.hpp"

#define NUM_BITS 100 

using namespace std;
using namespace boost::bloom_filters;

int main(int argc, char** argv) {

   serializeable_counting_bloom_filter<string, NUM_BITS> bf;
   const string sentence("See Spot Run");
   stringstream ss(sentence);

   string buf;
   while( ss >> buf ) { 
      cout << buf << '\t' << bf.probably_contains(buf) << endl << bf.str() << endl; 
      bf.insert(buf); 
   }
   
}

