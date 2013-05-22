// Mode is: -*- C++ -*- 
// --- General.cc --- Created: Sun Sep 21 16:12:10 1997
// 
// Time-stamp: <September 21, 1997 20:44:31 paolo@tequila.dsi.unifi.it>
// Update Count: 4
// 
// RCS: $Id: General.cxx,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   $Log: General.cxx,v $
// RCS: $Id: General.cxx,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   Revision 1.1  2003/11/19 23:42:16  jcheng
// RCS: $Id: General.cxx,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   *** empty log message ***
// RCS: $Id: General.cxx,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $
// RCS: $Id: General.cxx,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   Revision 1.1  2003/07/21 21:02:13  jcheng
// RCS: $Id: General.cxx,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $   *** empty log message ***
// RCS: $Id: General.cxx,v 1.1 2003/11/19 23:42:16 jcheng Exp $   $Locker:  $

// = FILENAME
//    General.cc
// = LIBRARY
// 
// = AUTHOR
// Paolo Frasconi (paolo@mcculloch.ing.unifi.it)
// = COPYRIGHT
// Copyright (C) 1997 Paolo Frasconi

#include "General.h"
#include <math.h>

int debugLevel=3;


int 
argmax(const Float* x, int n)
{
  Float xmax=x[0];
  int amax=0;
  for (int i=1; i<n; i++) {
    if (x[i]>xmax) {
      xmax=x[i];
      amax=i;
    }
  }
  return amax;
}

void
meanvar(Float* v, int n, Float& mean, Float& std)
{
  mean = 0.0;
  int i;
  for (i=0; i<n; i++) {
    mean += v[i];
  }
  mean /= Float(n);
  std = 0.0;
  for (i=0; i<n; i++) {
    std += (mean - v[i])*(mean - v[i]);
  }
  std /= Float(n-1);
  std = (Float)sqrt(std);
}

void
printvec(const Float* x, int n)
{
  for (int i=0; i<n; i++) {
    cout << x[i] << " ";
    if (i%8==7)
      cout << "\n";
  }
  cout << "\n";
}

#ifdef NEVER
#include <math.h>
int
matherr(struct exception *x) 
{
  char msg[1024];
  
  switch (x->type) {
  case DOMAIN:
    sprintf(msg, "DOMAIN");
    break;
  case SING:
    sprintf(msg, "SING");
    break;
  case OVERFLOW:
    sprintf(msg, "OVERFLOW");
    break;
  case UNDERFLOW:
    sprintf(msg, "UNDERFLOW");
    break;
  case TLOSS:
    sprintf(msg, "TLOSS");
    break;
  case PLOSS:
    sprintf(msg, "PLOSS");
    break;
  }
  if (!strcmp(x->name, "log")) {
    explain = "You probably need to increase IOHMM.scalingFactor";
  }
  FAULT("Type=" << msg << "\nname = " << x->name 
	<< " args= " << x->arg1 << " " << x->arg2
	<< "\nHint: " << explain);
  return 0;
}

#endif
