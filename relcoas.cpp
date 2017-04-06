
// g++ -std=c++0x -O3 -I/home/christoph/dlib-19.1/ -o optimize optimization_relatedness_realdata.cpp
//time ./optimize_release -indiv 1,2,46,47 -in_file Software_Release/release_ws-0.txt -meta_file Software_Release/release_metadata_ws-0.txt  -max_e 0.001 -max_c 0.25

#include <dlib/optimization.h>			// optimization library


using namespace std;

typedef dlib::matrix<double,0,1> column_vector;					// column vector for lbfgs matrix of 0 rows and 1 column

struct parameters {
	int lrij,N,L,lcontam,si,ei;
	vector <double> pk,qk,fk; 
	vector < vector <int>> haps; 
	vector <int> individuals;
	vector <double> sll_vector,p_temp;
	int ws;
}sp;



vector <string> splitstr2vec(string s,string delim) {
	size_t pos = 0;
	string token;
	vector <string> ret;
	while ((pos = s.find(delim)) != string::npos) {
		token = s.substr(0, pos);
		ret.push_back(token);
		s.erase(0, pos + delim.length());
	}
	ret.push_back(s);
	return ret;	
}

vector <int> splitstr2vec2 (string s) {
	vector <int> ret;
	for(string::iterator it = s.begin(); it != s.end(); ++it) {
		ret.push_back(stoi(string(1 , *it )));
	}
	return ret;
}

void fill_cv_with_randomnumbers(int l, column_vector & pv){
	random_device rnd_device;
	mt19937 mersenne_engine(rnd_device());
	uniform_real_distribution<double> dist(0,1);
   auto gen = bind(dist, mersenne_engine);
   vector<double> rands(l);
   generate(begin(pv), end(pv), gen);
}

double mn(double a, double b) {
	return (a+b)/2.0;
}
void finalize_results(column_vector & pv){
		
			cout<<"###	Estimated relatedness coefficients r_ij"<<endl;
			cout<<"\t";
			for(int i=0;i<sp.N;i++){
				if(sp.individuals.size()==1 && sp.N==sp.individuals[0]){
					cout<<(i+1)<<"\t";
				}else{
					cout<<sp.individuals[i]<<"\t";
				}				
			}cout<<endl;
			int pc=0;
			double c1=0.0,c2=0.0;
			for(int rc=0;rc<sp.N;rc++){
				if(sp.individuals.size()==1 && sp.N==sp.individuals[0]){
					c1=pv(sp.lrij+rc);
				}else{
					c1=pv(sp.lrij+sp.individuals[rc]-1);
				}
				for(int cc=0;cc<sp.N;cc++){
					if(sp.individuals.size()==1 && sp.N==sp.individuals[0]){
						c2=pv(sp.lrij+cc);
					}else{
						c2=pv(sp.lrij+sp.individuals[cc]-1);
					}
					if(cc==0){
						if(sp.individuals.size()==1 && sp.N==sp.individuals[0]){
							cout<<(rc+1)<<"\t";
						}else{
							cout<<sp.individuals[rc]<<"\t";
						}	
					}	
					if(cc<rc){
						cout<<"-\t";
					}else if(cc==rc){	
						cout<<"-\t";
					}else{
						printf ("%1.2f\t", (double)pv(pc)/(1-mn(c1,c2)));	// 
						pc++;
					}
					if(cc==(sp.N-1))
						cout<<endl;
				}
			}
			
			cout<<"\n###	Contamination rate estimates C_i"<<endl;
			for(int rc=pc;rc<pc+sp.lcontam;rc++){
				printf ("%1.2f ", pv(rc));
			}cout<<endl;
			pc+=sp.lcontam;
			
			
			
			cout<<"\n###	Sequencing error estimates e "<<endl;
			cout<<pv(pc)<<endl;
			pc+=1;
			

}


double llhood_ij (int lrij,const vector <double> & pk,const vector <double> & qk,const vector <double> & fk,int i1,int j1,int i2,int j2,const vector < vector <int>> & haps,const column_vector& param_vector,int N, int L,int contamsize) {
	
	double ci,cj,r,e,C,i3,j3,llhoods;	
	
	i3=i1+1;
	j3=j1+1;
	r = param_vector((i3-1)*N+j3-((i3*i3+i3)/2)-1);
	
	ci = param_vector(lrij+i2);
	cj = param_vector(lrij+j2);
	e = param_vector(lrij+contamsize+1-1);
	double sum = accumulate(param_vector.begin()+lrij, param_vector.begin()+lrij+contamsize, 0.0);
	C = sum / (double)contamsize;
	
	double P11=0.0,P10=0.0,P01=0.0,P00=0.0,Q11=0.0,Q10=0.0,Q01=0.0,Q00=0.0,R11=0.0,R10=0.0,R01=0.0,R00=0.0;
	double r1=(1-(r/2));
	double r2=r/2;
	double c1=(1-ci)*(1-cj);
	double c2=(ci*(1-cj)+cj*(1-ci));
	double c3=ci*cj;
	double c4=ci*(1-cj);
	double c5=cj*(1-ci);
	double ee=e*e;
	double ee2=(1-e)*(1-e);
	double ee3=(1-e)*e;
	
	
	
	double p=0.0,f=0.0,q=0.0,tmp=0.0,p2=0.0,f2=0.0;
	for(int site=0; site<L;site++) {
		f=fk[site],f2=1-fk[site],q=qk[site],tmp=0.0;
		p=(q-C*f)/(1-C);
		if(p>1.0)
			p=0.9999;
		if(p<0)
			p=0.001;
		
		p2=1-p;
		P11=r1*p*p + r2*p;
		P10=P01=r1*p*p2;
		P00=r1*p2*p2 + r2*p2;
		
		Q11=c1*P11 + c2*p*f + c3*f*f;
		Q10=c1*P10 + c4*p2*f + c5*p*f2 + c3*f*f2;
		Q01=c1*P01 + c4*p*f2 + c5*f*p2 + c3*f*f2;
		Q00=c1*P00 + c2*p2*f2 + c3*f2*f2;

		if(haps[i2][site]==1 && haps[j2][site]==1){
			tmp=ee2*Q11 + ee3*(Q10+Q01) + ee*Q00;
			if(tmp>0)
				R11+=log(tmp);
		}else if(haps[i2][site]==1 && haps[j2][site]==0){
			tmp=ee2*Q10 + ee3*(Q11+Q00) + ee*Q01;
			if(tmp>0)
				R10 +=log(tmp); 
		}else if(haps[i2][site]==0 && haps[j2][site]==1){
			tmp=ee2*Q01 + ee3*(Q11+Q00) + ee*Q10;
			if(tmp>0)
				R01 +=log(tmp); 
		} else if(haps[i2][site]==0 && haps[j2][site]==0){	
			tmp=ee2*Q00 + ee3*(Q10+Q01) + ee*Q11;
			if(tmp>0)
				R00 += log(tmp); 
		}

	}

	llhoods=R11+R10+R01+R00;
	return llhoods;
}


double llhood (const column_vector& param_vector) {
	double sll =0.0;
	if(sp.individuals.size()==1 && sp.N==sp.individuals[0]){
		for ( int j=sp.si;j<sp.N;j++) {
			for (unsigned int i=0;i<=(j+sp.ei);i++) {
				sll = sll + llhood_ij(sp.lrij,sp.pk,sp.qk,sp.fk,i,j,i,j,sp.haps,param_vector,sp.N,sp.L,sp.lcontam);
			}
		}
	}else{
		for (unsigned int j=sp.si;j<sp.individuals.size();j++) {	
			for (unsigned int i=0;i<=(j+sp.ei);i++) {
				sll=sll + llhood_ij(sp.lrij,sp.pk,sp.qk,sp.fk,i,j,sp.individuals[i]-1,sp.individuals[j]-1,sp.haps,param_vector,sp.N,sp.L,sp.lcontam);
			}
		}
	}
  return -sll;
}


int main(int argc, char* argv[]) {

string in_file,meta_file,indis,line;
vector <string> linevec;
vector <int> linevecint;
int ic=0;
double max_c,max_e;

column_vector param_vector;
column_vector low;
column_vector up;

//READ CMD LINE ARGUMENTS
if(argc !=11) {
	cout <<to_string(argc) +" There should be 5 parameters";
	exit(0);
}
for(int i = 1; i < argc; i+=2) {
	string a=argv[i];
    if(a == "-in_file") {
       in_file = argv[i+1];
    }else if(a == "-meta_file") {
       meta_file = argv[i+1];
    } else if (a == "-indiv") {
        linevec = splitstr2vec(argv[i + 1],",");
    }else if (a == "-max_e") {
         max_e = stod(argv[i+1]);
    }else if (a == "-max_c") {
         max_c = stod(argv[i+1]);
    }else{
		 cout <<argv[i]<< " Unknown parameter!";
		exit(0);
	 }
}

// FIGURE OUT NUMBER OF INDIVIDUALS
if(linevec.size()==1){
	try{
		sp.N=stoi(linevec[0]);
	}catch (exception const & e) {
       cout <<" ERROR:	-indiv contains invalid format"<< endl;
		 exit(0);
    }
    sp.individuals.push_back(sp.N);
}else{	
	sp.N=linevec.size();
	for(unsigned int i=0;i<linevec.size();i++){
		try{
			sp.individuals.push_back(stoi(linevec[i]));
		}catch (exception const & e) {
       cout <<" ERROR:	-indiv contains invalid format"<< endl;
		 exit(0);
    }
	}
	sort(sp.individuals.begin(), sp.individuals.end());	// sort if user did not sort them (e.g. 1,5,2,10,3)
}


//READ IN HAPLOTYPES
ifstream in( in_file );
if (in) {
	while (getline( in, line )) {
		if(line[0]=='0' || line[0]=='1' || line[0]=='2') {
			linevecint=splitstr2vec2(line);
			sp.haps.push_back(linevecint);
			ic+=1;
			linevecint.clear();
		}
   }
   in.close();
}else{
	cout<<"In_File not found"<<endl; 
	exit(0);
}

for(unsigned int i=0;i<sp.individuals.size();i++){
	if(sp.individuals.size()>1){
		if(sp.individuals[i]>ic) {
			cout<<"ERROR: Individual "<<sp.individuals[i]<<" not found in dataset"<<endl;
			exit(0);
		}
	}else{
		if(sp.individuals[i]!=ic){
			cout<<"ERROR: Dataset does not have "<<sp.individuals[i]<<" individuals"<<endl;
			exit(0);
		}
	}	
}


//READ IN METADATA
ifstream meta( meta_file );
if (meta) {
	while (getline( meta, line )) {
		if(isdigit(line[0])) {
			linevec=splitstr2vec(line,"\t");						//site	qk	fk
			try{
				sp.qk.push_back(stod(linevec[1]));
				sp.fk.push_back(stod(linevec[2]));
			}catch (exception const & e) {
				cout <<" ERROR:	Invalid format in "<< meta_file<<": "<<line<< endl;
				exit(0);
			}
		}
   }
   meta.close();
}else{
	cout<<"Meta_File not found"<<endl;
	exit(0);
}

//NUMBER OF SITES
sp.L=sp.qk.size();

//FILL THE PARAMETER VECTORS
int ltotal=0,lfreq=0,lerr=1;			

sp.lcontam=ic;
sp.lrij=(sp.N*(sp.N-1)/2);
sp.si=1;
sp.ei=-1;
sp.sll_vector.resize(sp.lrij);

ltotal=sp.lcontam+sp.lrij+lerr;
param_vector.set_size(ltotal);
low.set_size(ltotal);
up.set_size(ltotal);

fill_cv_with_randomnumbers(ltotal,param_vector);

fill(low.begin(),low.begin()+sp.lrij,0.0);
fill(low.begin()+sp.lrij,low.begin()+sp.lrij+sp.lcontam,0.0001);
fill(low.begin()+sp.lrij+sp.lcontam,low.begin()+sp.lrij+sp.lcontam+lerr,0.0001);

fill(up.begin(),up.begin()+sp.lrij,1);
fill(up.begin()+sp.lrij,up.begin()+sp.lrij+sp.lcontam,max_c);
fill(up.begin()+sp.lrij+sp.lcontam,up.begin()+sp.lrij+sp.lcontam+lerr,max_e);

cout<<"N="<<sp.N<<" lrij="<<sp.lrij<<" lcontam="<<sp.lcontam<<" lfreq="<<lfreq<<" L="<<sp.L<<endl;


// DLIB optimization
try {
	dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(100),dlib::objective_delta_stop_strategy(1e-10),llhood,derivative(llhood), param_vector,low,up);
	finalize_results(param_vector);
}catch (exception const & e)
    {
        cout << e.what() << endl;
    }
return 0;

}

