#include "../include/globals.hpp"

using namespace std;

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

/*
Matrix *Global::obs;

void Global::GEOS3(){
    Global::obs = new Matrix(46, 4);

    FILE* fid = fopen("../data/GEOS3.txt","r");
    if(fid==nullptr){
        printf("Error al abrir el archivo");
        exit(EXIT_FAILURE);
    }

    char* auxi;
    const char *tline;
    double Y, M, D, h, m, s, az, el, Dist;

    for (int i=0; i<=45; i++){
            fgets (auxi, 100 , fid);
        //no entiendo porqué esto no se castea a const char* aquí
            tline = reinterpret_cast<const char *>(auxi);
            Y = atof(reinterpret_cast<const char *>(tline[0, 1, 2, 3]));
            M = atof(reinterpret_cast<const char *>(tline[5, 6]));
            D = atof(reinterpret_cast<const char *>(tline[8, 9]));
            h = atof(reinterpret_cast<const char *>(tline[12, 13]));
            m = atof(reinterpret_cast<const char *>(tline[15, 16]));
            s = atof(reinterpret_cast<const char *>(tline[18, 19, 20, 21, 22, 23]));
            az = atof(reinterpret_cast<const char *>(tline[25, 26, 27, 28, 29, 30, 31, 32]));
            el = atof(reinterpret_cast<const char *>(tline[35, 36, 37, 38, 39, 40, 41]));
            Dist = atof(reinterpret_cast<const char *>(tline[44, 45, 46, 47, 48, 49, 50, 51, 52, 53]));
            Global::obs[i, 0] = Mjday::mjday(Y, M, D, h, m, s);
            Global::obs[i, 1] = Const::Rad*az;
            Global::obs[i, 2] = Const::Rad*el;
        Global::obs[i, 3] = 1e3*Dist;
    }

    fclose(fid);
}
*/

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

/*
void Global::GGMO3S(int c){
    Global::GGMO3Sdata = new Matrix(6, c);

    FILE* fid = fopen("../data/GGMO3S.txt","r");
    if(fid==nullptr){
        printf("Error al abrir el archivo");
        exit(EXIT_FAILURE);
    }

    for(int i=1; i<=c; i++){
        fscanf(fid,"%lf %lf %lf %lf %lf %lf",
               &((*Global::GGMO3Sdata)(1,i)),
               &((*Global::GGMO3Sdata)(2,i)),
               &((*Global::GGMO3Sdata)(3,i)),
               &((*Global::GGMO3Sdata)(4,i)),
               &((*Global::GGMO3Sdata)(5,i)),
               &((*Global::GGMO3Sdata)(6,i)));
    }

    fclose(fid);

}
 */
/*
Param* Global::getParam(){
    return new Param();
}

Matrix *Global::PC;
*/