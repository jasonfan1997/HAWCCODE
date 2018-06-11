#include <iostream>

//#include "TMath.h"

#include <fstream>

#include <string>
#include <vector>
#include <xcdf/XCDF.h>
#include <xcdf/XCDFField.h>
#include <xcdf/XCDFFile.h>
#include <xcdf/utility/XCDFUtility.h>

#include <sstream>

using namespace std;


template <class Type>
Type stringToNum(const string& str)
{
    istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}

int main(int argc, char** argv) {
  //char* name=getFileName(argv[0]);
  XCDFFile rec;
  rec.Open(argv[1], "r");
  XCDFFile recFile(argv[2],"w");
  FieldCopyBuffer buffer(recFile);
  XCDFFloatingPointField newfield = recFile.AllocateFloatingPointField("rec.proba", 0.000000001);
  std::set<std::string> fields;

  GetFieldNamesVisitor getFieldNamesVisitor(fields);
  rec.ApplyFieldVisitor(getFieldNamesVisitor);

  SelectFieldVisitor selectFieldVisitor(rec, fields, buffer);
  rec.ApplyFieldVisitor(selectFieldVisitor);

  ifstream file ( argv[3]);
  string value;
  while ( rec.Read() )
  {
     getline ( file, value);
     newfield <<stringToNum<double>(value);
     buffer.CopyData();

     recFile.Write();
   }
   rec.Close();
    recFile.Close();

  return 0;
}
