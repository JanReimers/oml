// File: UnitTest.C  Support routines.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include <string>
#include <cstring>

double eps=1e-14;

void StartTest(char* Member, char* Test)
{
  int nspace=78-strlen(Member)-8-strlen(Test)-6;

  std::cout << "-------------------------------------------------------------------------------"
	    << std::endl << std::endl;
  std::cout << "Member: " << Member;
  for (int i=0; i<nspace;i++) std::cout << " ";
  std::cout << "Test: " << Test << std::endl << std::endl;
}


void StartClass(const char* Class)
{
  size_t l=strlen(Class);
  size_t i=34-l/2;
  char* cline=new char[80];
  strcpy(cline,"                    |                          | ");
  strncpy(&(cline[i]),Class,l);

  std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;
  std::cout << "                    /--------------------------\\ " << std::endl;
  std::cout << "                    |                          | " << std::endl;
  std::cout << "                    |  Begin Diagnostics For:  | " << std::endl;
  std::cout << cline << std::endl;
  std::cout << "                    |                          | " << std::endl;
  std::cout << "                    \\--------------------------/ " << std::endl;
  std::cout << std::endl << std::endl;
  delete [] cline;
}

bool fancy_compare(std::stringstream& result,std::stringstream& expected)
{
   bool ret=true;
   int i=0;
//   std::cout << "fancy comparing '" << result << "' and '" << expected << "'" << std::endl;
   std::string next_res("");
   std::string next_exp("");
   do
   {
     result   >> next_res;
     expected >> next_exp;
      if (next_exp!="*") ret=ret&&(next_res==next_exp);
//      std::cout << "comparing '" << next_res << "' and '" << next_exp << "' " << ret << std::endl;
      i++;
   } while (!result.eof() && !expected.eof());
   return ret;
}
