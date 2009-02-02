#include <string>
#include <iostream>
using namespace std;

// Auxiliary stuff. Perhaps move to separate file.
int lookup_format(const string *supported_formats, const string& path)
{
  size_t suffix_start = path.rfind('.')+1;
  string suffix = path.substr(suffix_start);
    
  cout << "Suffix is " << suffix << endl;    

  for(size_t i=0;!supported_formats[i].empty(); i++)
    if(supported_formats[i] == suffix) return i;

  return -1;
}
