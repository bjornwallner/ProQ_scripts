// Mode is: -*- C++ -*- 
// --- Cpp.h --- Created: Thu Sep 21 12:15:12 1995
// 
// Time-stamp: <October 24, 1996 17:36:15 paolo@dsi>
// Update Count: 2
// 
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   $Log: Cpp.h,v $
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   Revision 1.1  2003/11/19 23:42:16  jcheng
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   *** empty log message ***
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   Revision 1.1  2003/07/21 21:02:13  jcheng
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   *** empty log message ***
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $
// * Revision 1.3  1996/10/24  14:54:52  paolo
// * *** empty log message ***
// *
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   Revision 1.2  1996/03/20 10:12:53  paolo
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   *** empty log message ***
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   Revision 1.1  1996/02/01 15:49:18  paolo
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   Initial revision
// RCS: $Id: Cpp.h,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $

// = FILENAME
//    Cpp.hh
// = LIBRARY
// 
// = AUTHOR
// Paolo Frasconi (paolo@mcculloch.ing.unifi.it)
// = COPYRIGHT
// Copyright (C) 1995 Paolo Frasconi


#ifndef Cpp_h
#define Cpp_h 1

#include "General.h"


class Cpp_istreambase {
protected:
  filebuf inbuf;
  char  *gzip_tmp_fname;
  char  *cpp_tmp_fname;
  bool  Unzipped;
  Cpp_istreambase(const char* fname);
};

class Cpp_istream : public istream, public Cpp_istreambase {
private:
public:
  Cpp_istream(const char *fname);
  ~Cpp_istream();
  void close();
};


#endif // Cpp_h
