#include <iostream>
#include <math.h>
#include <TRandom1.h>

using namespace std;

int main(void)
{
	TRandom1 * myrand = new TRandom1();
	for(int i = 0; i < 10; ++i)
	{
		cout << myrand->Gaus(5,1) << endl;
	}
	cout << sqrt(4.2) << endl;
	return 0;
}