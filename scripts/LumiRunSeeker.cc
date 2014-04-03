#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>
//#include <boost/regex.hpp>

using namespace std;

int wildcmp(const char *, const char *);

int main(int argc, char **argv) {

  string file;
  file.assign(argv[1]);
  string trigger;
  trigger.assign(argv[2]);
  cout << "Will search for runs with trigger "<<trigger<<" available" <<endl;
  map<string,vector<int> > RunNbPerTrigger;
  map<int,vector<string> > TriggerPerRunNb;
  
  ifstream TriggerSummaryFile(file.c_str(),ifstream::in);
  string dummy, triggername;
  int runnb, triggerbit, nfired, npassed, idummy;

  if ( TriggerSummaryFile )
  {
	while ( !TriggerSummaryFile.eof() ) { 
		TriggerSummaryFile >> dummy >> runnb >> triggerbit >> nfired >> npassed >> idummy >> triggername ;
		//cout << runnb << " " << triggername << endl;
		RunNbPerTrigger[triggername].push_back(runnb);
	}
  }
  else
  {
	cout << "File " << file << " does not exist!" << endl;
  }

  //size_t found;
  for(map<string,vector<int> >::iterator it = RunNbPerTrigger.begin(); it != RunNbPerTrigger.end(); it++)
  {
	//boost::regex xtrigger(trigger.c_str());
	//boost::regex xtrigger(argv[2]);
	//if (boost::regex_search(it->first, xtrigger))

	//it->first.find(trigger);
	//if (found!=string::npos)
	if(wildcmp(trigger.c_str(), it->first.c_str()))
	{
		sort(it->second.begin(),it->second.end());
		cout<<"Trigger "<<it->first<<" available for runs "<<it->second.front()<<"-"<<it->second.back()<<endl;
	}
  }
  
}

int wildcmp(const char *wild, const char *string) {
  // Written by Jack Handy - jakkhandy@hotmail.com

  const char *cp = NULL, *mp = NULL;

  while ((*string) && (*wild != '*')) {
    if ((*wild != *string) && (*wild != '?')) {
      return 0;
    }
    wild++;
    string++;
  }

  while (*string) {
    if (*wild == '*') {
      if (!*++wild) {
        return 1;
      }
      mp = wild;
      cp = string+1;
    } else if ((*wild == *string) || (*wild == '?')) {
      wild++;
      string++;
    } else {
      wild = mp;
      string = cp++;
    }
  }

  while (*wild == '*') {
    wild++;
  }
  return !*wild;
}
