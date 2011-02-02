// File: diag.h  Header file for all diagnostic routines.
#ifndef _diag_h_
#define _diag_h_

// Copyright (1994-2005), Jan N. Reimers

#include <iostream>
#include <sstream>
#include <string>

extern double eps;


void StartTest (char* Member, char* Test);
void StartClass(const char* Class);

std::string Extract(const std::string& s,char delimiter,int occurance);
bool fancy_compare(std::stringstream& result,std::stringstream& expected);

#define EXPECT(Test,Expected) {pass=pass&&DoTest(Test,Expected,#Test);}
#define EXPECT1(Test,Expected,TestString) {pass = pass && DoTest(Test,Expected,TestString);}
#define EXPECT2(Test,Function,Expected) {Test.Function;pass = pass && DoTest(Test,Expected,#Function);}
#define EXPECT3(Test,Function,Expected) {Function;pass = pass && DoTest(Test,Expected,#Function);}




template <class Ob, class T> bool DoTest(const Ob& result,const T& expected,const char* name)
{
//  cout << (void*)&result << " " << result << endl;
  bool pass;
  {
    std::stringstream res_stream, exp_stream;
    res_stream << result << std::ends;
    exp_stream << expected << std::ends;
    pass=fancy_compare(res_stream,exp_stream);
  }
  if (!pass)
  {
    std::stringstream res_stream, exp_stream;
    res_stream << result << std::ends;
    exp_stream << expected << std::ends;    
    std::string res_string,exp_string;
    getline(res_stream,res_string,(char)0);
    getline(exp_stream,exp_string,(char)0);
    std::cerr << "Test " << name << " FAILED:" << std::endl;
    std::cerr << "  Expected :'" << exp_string << "'" << std::endl;
    std::cerr << "  Got      :'" << res_string << "'" << std::endl;
  }
  else 
  {
    std::cout << "Passed: " << name << std::endl;
  }
  return pass;
   
}

#endif //_diag_h_
