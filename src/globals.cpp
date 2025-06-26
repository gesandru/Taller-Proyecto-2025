#include "../include/globals.hpp"

using namespace std;

Param AuxParam;

void AParam(){
	AuxParam.Mjd_UTC = 0;
	AuxParam.n = 0;
	AuxParam.m = 0;
	AuxParam.sun = 0;
	AuxParam.moon = 0;
	AuxParam.planets = 0;
	AuxParam.Mjd_TT = 0;
}

    Matrix eopdata;

    void eop19620101(int c){
        eopdata = zeros(13, c);

        FILE *fid = fopen("../data/eop19620101.txt","r");
        if(fid==NULL){
            printf("Error al abrir el archivo eop19620101");
            exit(EXIT_FAILURE);
        }

for(int i=1; i<=c; i++){
fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
&((eopdata)(1,i)),
&((eopdata)(2,i)),
&((eopdata)(3,i)),
&((eopdata)(4,i)),
&((eopdata)(5,i)),
&((eopdata)(6,i)),
&((eopdata)(7,i)),
&((eopdata)(8,i)),
&((eopdata)(9,i)),
&((eopdata)(10,i)),
&((eopdata)(11,i)),
&((eopdata)(12,i)),
&((eopdata)(13,i)));
}

fclose(fid);
}

Matrix Cnm;
Matrix Snm;

void GGMO3S(){
    Cnm = zeros(181, 181);
    Snm = zeros(181,181);

    FILE *fid = fopen("../data/GGM03S.txt","r");
    if(fid==NULL){
        printf("Error al abrir el archivo GGM03S.txt");
        exit(EXIT_FAILURE);
    }

    Matrix temp = zeros(6);

    for(int n=0; n<=180; n++){
        for (int m=0; m<=n; m++){
                fscanf(fid, "%lf %lf %lf %lf %lf %lf",
				&((temp)(1)),
				&((temp)(2)),
				&((temp)(3)),
				&((temp)(4)),
				&((temp)(5)),
				&((temp)(6)));
                Cnm(n+1, m+1) = temp(3);
                Snm(n+1, m+1) = temp(4);
        }
    }

    fclose(fid);

}

Matrix PC;

void DE430Coeff(){
    PC = zeros(2285, 1020);

    FILE *fid = fopen("../data/DE430Coeff.txt","r");
    if(fid==NULL){
        printf("Error al abrir el archivo DE430Coeff.txt");
        exit(EXIT_FAILURE);
    }

    Matrix temp = zeros(6);

	for(int i=1; i<=2285; i++){
		for(int j=1; j<=1020; j++){
			fscanf(fid,"%lf", &((PC)(i,j)));
		}
	}
	
fclose(fid);        
}

Matrix obs;

void GEOS3(int f){
    obs = zeros(f, 4);

ifstream fid("../data/GEOS3.txt");
        if(!fid){
            printf("Error al abrir el archivo GEOS3");
            exit(EXIT_FAILURE);
        }
		
string tline;
int Y, M, D, h, m, s;
double az, el, Dist; 

for(int i=1; i<=f; i++){
    getline(fid, tline);
    Y = stoi(tline.substr(0,4));
    M = stoi(tline.substr(5,2));
    D = stoi(tline.substr(8,2));
    h = stoi(tline.substr(12,2));
    m = stoi(tline.substr(15,2));
    s = stoi(tline.substr(18,6));
    az = stod(tline.substr(25,8));
    el = stod(tline.substr(35,7));
    Dist = stod(tline.substr(44,10));
    obs(i,1) = Mjday(Y,M,D,h,m,s);
    obs(i,2) = Const::Rad*az;
    obs(i,3) = Const::Rad*el;
    obs(i,4) = 1e3*Dist;
}

fid.close();      
}