#include <string>
using namespace std;

class StringParse
{
	public:
		StringParse();
		~StringParse();
		double * NumberGet();
		char NumberExtract(string str);


	private:
		char HeaderExtract(string str);
		char ContentsExtract(string str);
		// string str;
		double value[4];
};