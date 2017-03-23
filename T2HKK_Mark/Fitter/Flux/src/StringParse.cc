#include "MakeFlux.h"
#include "StringParse.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
using namespace std;


/**
 * It will take string and get number information inside of string
 */
StringParse::StringParse()
{
}

char StringParse::NumberExtract(string str)
{
	if(str.at(0) == 'a')
	{
		return HeaderExtract(str); 		// starts with 'a' is header of each contents including angle value
	}
	else
	{
		return ContentsExtract(str);	// return energy and flux 
	}
	
}

StringParse::~StringParse()
{

}

double * StringParse::NumberGet()
{
	return value;
}

char StringParse::HeaderExtract(string str)
{
	char c;
	int strLength = str.length();
	int i = 0, j = 0;
	int dotLocation = 0;
	int IntCheck;
	double tmp = 0;
	string tmpStr;

	// To get doecimal value, find "."	
	for(int i = 0; i < 2; i++)
	{
		dotLocation = str.find(".", dotLocation + 1);
		if(str.at(dotLocation - 2) == '-')
		{
			tmpStr = str.substr(dotLocation-2, 5);	
		}
		else
		{
			tmpStr = str.substr(dotLocation-1, 4);
		}
		
		value[i] = atof(tmpStr.c_str());
	}
	i = 2;
	j = dotLocation + 4;
	while(j < (strLength))
	{
		c = str.at(j);
		if(c >= 0x30 && c <= 0x39)
		{
			IntCheck = 1;
			c -= 0x30;
			tmp = tmp * 10 + c;
		}
		else
		{
			if(IntCheck == 1)
			{
				// cout << tmp << endl;
				IntCheck = 0;
				value[i++] = tmp;
				tmp = 0;
			}
		}
		j++;
	}
	return 1;
}

char StringParse::ContentsExtract(string str)
{
	int length = 0;
	double tmpValue;
	stringstream linestream(str);
	while(linestream >> tmpValue)
	{
		if(tmpValue != 0)
		{
			value[length] = tmpValue;
			cout << value[length] << "\t\t";
			length++;
		}
	}
	if(length == 0)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}
