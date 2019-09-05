#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]){
	//char filepath[256] = "/home/szhang/work/ampt-initial-fluctuation/analysis-flow/ampt-initial-fluctuation-DATA/most-central/c12-Au197-200GeV/data-Chain/";
	//char filepath[256] = "/home/szhang/work/ampt-initial-fluctuation/analysis-flow/ampt-initial-fluctuation-DATA/most-central/c12-Au197-200GeV/data-nofluc/";
	//char filepath[256] = "/home/szhang/work/ampt-initial-fluctuation/analysis-flow/ampt-initial-fluctuation-DATA/most-central/c12-Au197-200GeV/data-Triangle/";
	//char filepath[256] = "/home/hejunjie/161208AMPT-HBTRadius/analysis-HBTRadius/ampt-initial-fluctuation-DATA/most-central/c12-Au197-200GeV/data-Chain/";
	//char filepath[256] = "/home/hejunjie/161208AMPT-HBTRadius/analysis-HBTRadius/ampt-initial-fluctuation-DATA/most-central/c12-Au197-200GeV/data-nofluc/";
	//char filepath[256] = "/home/hejunjie/161208AMPT-HBTRadius/analysis-HBTRadius/ampt-initial-fluctuation-DATA/most-central/c12-Au197-200GeV/data-Triangle/";
	int n = 0;
	int m=0;

	ofstream listout;
	listout.open("rawdatapath.list");
	if(listout.fail()){
		cout<<"Sorry, couldn¡¯t open file"<<endl;
		exit(1);
	}

	for(int i1 = 1; i1 < 8; i1++){
		for(int j1 = 0; j1 < 1000; j1++){
			listout<<argv[1]<<"run"<<i1<<"/data-"<<j1<<".root"<<endl;
//			listout<<argv[1]<<"data-"<<j1<<".root"<<endl;
			n++;
		}
	}
	listout.close();

if (argc >=3){
	ofstream listout;
	listout.open("inputpath.list");
	if(listout.fail()){
		cout<<"Sorry, couldn¡¯t open file"<<endl;
		exit(1);
	}

	for(int i1 = 0; i1 < 10000; i1++){
		listout<<argv[2]<<"data-"<<i1<<".dat"<<endl;
		m++;
	}
	listout.close();
}

if (argc ==4){
	ofstream listout;
	listout.open("outputdirpath.list");
	if(listout.fail()){
		cerr<<"Sorry, couldn¡¯t open file"<<endl;
		exit(1);
	}

	for(int i1 = 0; i1 < 10000; i1++){
		listout<<argv[3]<<"data-"<<i1<<endl;
	}
	listout.close();
}

	//pi
	/*//most-central
	for(int i1 = 1; i1 < 2; i1++){
		for(int j1 = 0; j1 < 400; j1++){
			listout<<filepath<<"run"<<i1<<"/data-"<<j1<<".root"<<endl;
			n++;
		}
	}

	//minibias
	for(int i1 = 1; i1 < 5; i1++){
		for(int j1 = 0; j1 < 200; j1++){
			listout<<filepath<<"run"<<i1<<"/data-"<<j1<<".root"<<endl;
			n++;
		}
	}
	for(int i2 = 5; i2 < 6; i2++){
		for(int j2 = 0; j2 < 200; j2++){
			listout<<filepath<<"run"<<i2<<"/data-"<<j2<<".root"<<endl;
			n++;
		}
	}
	for(int i3 = 5; i3 < 6; i3++){
		for(int j3 = 0; j3 < 1000; j3++){
			listout<<filepath<<"run"<<i3<<"/data-"<<j3<<".root"<<endl;
			n++;
		}
	}*/


	cout<<"there are "<<n<<" .root files"<<endl;
	if(argc>=3)cout<<"there are "<<m<<" .dat files"<<endl;
	if(argc==4)cout<<"there are "<<m<<" save-directions"<<endl;
	cout<<"//////////////////////////////////////////////"<<endl;
	cout<<endl;
	cout<<"this program finished at "<<endl;
	system("date");
	cout<<endl;
	cout<<"//////////////////////////////////////////////"<<endl;

	return 0;
}
