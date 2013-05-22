#!/usr/bin/perl -w
while(<>)
{
    if(/ATOM/)
    {
	substr($_,21,1)=" ";
    }
    print;
}
