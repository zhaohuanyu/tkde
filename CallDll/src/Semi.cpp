#include <iostream>
using namespace std;
void semiOpt(string fileName, string resultName, string connectedIndex,string connectedValue, double sigma);
int main() {
	string strfile="OliveOil";
	string strresult="result.txt";
	string strconnect="value.txt";
	string strindex="index.txt";
	semiOpt(strfile, strresult, strindex , strconnect, 1.24088);
	return 1;
}
