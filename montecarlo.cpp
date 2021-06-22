
//Usage: ./montecarlo 
#include <iostream>
#include<fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

//// Función que hace la distancia en un medio periódico ////
double distmic(double r1[3], double r2[3], double lbox){

  double dx=r1[0]-r2[0];
  dx=dx-floor(dx/lbox+0.5)*lbox;
  double dy=r1[1]-r2[1];
  dy=dy-floor(dy/lbox+0.5)*lbox;
  double dz=r1[2]-r2[2];
  dz=dz-floor(dz/lbox+0.5)*lbox;
  double d=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
  
  //  std::cout<<"d ="<<d<<"\n";
  return d;
}


//// Fórmula de Lennard Jones ////
  double enerLJ(double d){
  double sigma=1;
  double epsilon=1;
  double rcut=2.5*sigma;
  double  E=4*epsilon*(pow(sigma/d,12)-pow(sigma/d,6))-4*epsilon*(pow(sigma/rcut,12)-pow(sigma/rcut,6));
  //   
   if (d>rcut){
 double E=0;
   }
  
  return E;

  }

////Formula de LJ derivada

double enerLJderivative(double d){
  double sigma=1;
  double epsilon=1;
  double rcut=2.5*sigma;

  //    double E=0;
  // if (d<rcut){
  double dE=4*epsilon*(1/d)*(6*pow(sigma/d,6)-12*pow(sigma/d,12));
  // }
  
  return dE;
}

////Energía de la partícula i en el sistema de todas las particulas
double energypart(int Npart, int ipart, double ri[3],  std::vector<double> &pos, double lbox){


  double ener=0;
  double dist; //hay que definirla porque sino no va a ser capaz de crearla a partir de la funcion
  double rj[3];
  

  for(int j=0; j<Npart; j++){
    if(j!=ipart){
      rj[0]=pos[3*j];
      rj[1]=pos[3*j+1];
      rj[2]=pos[3*j+2];
      dist=distmic(rj, ri,lbox);// ri y rj tienen 3 argumentos cada uno. lbox es solo un argumento

      ener=ener+enerLJ(dist);
    }
    
    
      //   std::cout<<"dist "<<dist<<"\n";     
      // std::cout<<"EnerLJ(dist) "<<enerLJ(dist)<<"\n";   
   }
  //   std::cout<<"Ener "<<ener<<"\n";   
  return ener;
}

////Radial distribution function ////

void radial_distr(std::vector<double> &pos, double lbox,  double rho, int sw, double Temp ){
 //si sw==0, inicializo todo.
 //Si sw==1 va sumando g(r) para calcular histograma
  //Si sw==2, normaliza e imprime para dibujar

  int Npart=pos.size()/3;
  double ri[3];
  double rj[3];
  static  int counter=0;
  double lower=0;
   int  upper=round(sqrt(3*pow(lbox/2,2))); //segun creo, este es el valor maximo que puede tomar dis
  
  double dr=0.05;
  int nbin=(upper-lower)/dr;
  static std::vector<double> histogram(nbin,0); //el static sirve para que se quede guardado el valor para las siguientes veces que llamemos a la funcion

  if (sw==1){

    counter++;
  for (int i=0; i< Npart-1; i++){
      ri[0]=pos[3*i];
      ri[1]=pos[3*i+1];
      ri[2]=pos[3*i+2];
      for (int j=i+1; j<Npart; j++){
      rj[0]=pos[3*j];
      rj[1]=pos[3*j+1];
      rj[2]=pos[3*j+2];
      double dist=distmic(ri,rj,lbox);
     
      int bin = nbin * (dist - lower)/(upper-lower);
      if (bin < nbin and bin>=0){
	histogram[bin] +=2;

      }
      
    }
  }

  
  }
  else if (sw==2){
    std::ofstream out2("rdf.res");
    
    double integral=0;
    double integral2=0;   
    // double dx = (upper - lower)/nbin;
    for (int i =0; i<nbin; i++){
      double x= ((i+0.5)/nbin)*(upper-lower) + lower;//centro de cada bin del histograma
      double dV=4*M_PI/3*(pow((i+1)*dr,3)-pow(i*dr,3));
      out2<<x<<" "<<histogram[i]/(counter*dV*Npart*rho)<<"\n";
      integral=integral+(histogram[i]/(counter*dV*Npart*rho))*pow(x,2)*enerLJ(x)*dr;
      integral2=integral2+(histogram[i]/(counter*dV*Npart*rho))*pow(x,3)*enerLJderivative(x)*dr;
      //      std::cout<<integral2<<"\n";
    }
    double e=Temp*3/2+2*M_PI*rho*integral;
    double u= 2*M_PI*rho*integral;
    double pres=rho*Temp-integral2*2*M_PI*pow(rho,2)/3;
    double pres_excess=-integral2*2*M_PI*pow(rho,2)/3;
    std::cout<<"Energia por particula= "<<e<<"\n";
    std::cout<<"Presion="<<pres<<"\n";
    std::cout<<"Exceso de presion= "<<pres_excess<<"\n";
    std::cout<<"Exceso energia por particula= "<<u<<"\n";
  }
}




 
//// Aqui empezamos con el programa ////
int main(){

  std::vector<double> pos; //vector de posiciones
  int Nstep=100000;
  double rho=0.6;
  double Temp=5;
  std::cout<<"temperature= " <<Temp<<"\n";
  double lbox=10;
  double Npart=rho*pow(lbox,3);
  //  std::cout<<"Npart="<<Npart<<"\n";
  int Nsample;
  int Nvideo;
  double djump=0.1; //tamaño del salto en unidades de sigma (sigma es el diametro de la particula)
  int nyes=0;
  int ntot=0;
  double rold[3];
  double rnew[3];
  double Enew;
  double Eold;
  double Esample;  
  int ipart; //indice de la particula
  double accratio=0.5;//preguntar a alguien
  double rr;
  std::ofstream out("energyTmuyalta.res");//para guardar los datos de las energias en un archivo
  std::ofstream posout("posiciones.res");//para guardar los datos de las posiciones en un archivo

  //generador de numeros aleatorios
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis11(-1.0,1.0);
  std::uniform_real_distribution<double> dis01(0,1.0);
  std::uniform_int_distribution<int> disN(0,Npart-1);
  
  //lectura de datos de un archivo (pos)

   // std::ifstream input("mc_dens05.rst");//montecarlo restart
  // input>>Npart>>lbox>>useless1>>useless2>>djump;
   // pos.resize(3*Npart);

   //   if (not input.good()){
   



  //// Generamos la primera configuracion de forma aleatoria //// 
    for(int i=0; i<3*Npart; i++){
      pos.resize(3*Npart);
      pos[i]=dis11(gen)*lbox*0.5; //coordenada desde el centro de la caja
     
     }
    
  //   //    }
   
  //    //  for(int i=0;i<3*Npart; i++){
  //   //   input>>pos[i];
    
  //   // }
  //   // input.close();
  
  // //lectura de parametros (Nstep, Nsample, Temp)

   Nsample=100;
   Nvideo=10000;
  

  
   for (int i=0; i<Nstep; i++){


    
     ipart = disN(gen);   
     rold[0]=pos[ipart*3];
     rold[1]=pos[ipart*3+1];
     rold[2]=pos[ipart*3+2];
     rnew[0]=rold[0]+dis11(gen)*djump;
     rnew[1]=rold[1]+dis11(gen)*djump;
     rnew[2]=rold[2]+dis11(gen)*djump;
     //     std::cout<<"pos"<<rold[0]<<" "<<rold[1]<<" "<<rnew[0]<<" "<<rnew[1]<<"\n";
    
     Eold=energypart(Npart,ipart, rold, pos,lbox);
     Enew=energypart(Npart,ipart, rnew, pos,lbox);


     
     // std::cout<<"Eold y Enew  "<<Eold<<" "<<Enew<<"\n";
    double P = exp(-(Enew-Eold)/Temp);//probabilidad entre el estado anterior y el intento
    
    //        std::cout<<"P  "<<P<<"\n";
    //   
    double r=dis01(gen);
    //  std::cout<<"r"<<r<<"\n";
      if(P>r){
	pos[ipart*3] = 	rnew[0];
	pos[ipart*3+1]=	rnew[1];
	pos[ipart*3+2]=	rnew[2];
	
	nyes+=1;
	ntot+=1;
	Esample=Enew;
	
      }
      else {
	    ntot+=1;
	    Esample=Eold;
	    

      }
      // std::cout<<"Esample"<<Esample<<"\n";  
       
    //Escribimos los resultados (evolucion de la energia)

    if (i%Nsample==0){
      out<<Esample*0.5<<std::endl;//los guarda en out, definido anteriormente
      //      std::cout<<Esample<<"\n";


      //tambien dentro de este bucle vamos a ir creando el histograma
     
    
      
    }
      if (i%Nvideo==0){
	 posout<<"\n";	     
	 for (int j=0 ;j<Npart; j++){
	      
	   posout<<" "<<pos[3*j]<<pos[3*j+1]<<pos[3*j+2]<<"\n";//<<std::endl; //guarda las posiciones en posout
	   	 }
	  }	 
     //control del salto
      rr=nyes/double(ntot);

      if (i>1000){
       if (rr > accratio){
        djump=djump*1.005;
      }
      if (rr < accratio){
        djump=djump/1.005;
      }
 
      }
      //         std::cout<<" "<<rr<<"\n";   
  }
  radial_distr(pos,lbox,rho,1,Temp);
  radial_distr(pos,lbox,rho,2,Temp);

  //energia media calculada "manualmente"
  double eprima;
  double eprima_cuad;
  int N=1;
  double  eprimamedia=0;
  double eprima_cuad_media=0;
  double var;
  double cv;
    for(int i=0; i<Npart-2; i++){
      double R[3];
       R[0]=pos[i*3];
       R[1]=pos[i*3+1];
       R[2]=pos[i*3+2];
       eprima=energypart(Npart, i, R,pos,lbox)/2;//hay que dividir por dos porque esta seria la energia de la pareja entera
       eprima_cuad=pow(eprima,2); 
       eprimamedia=eprimamedia*(N-1)/N+eprima/N;
       eprima_cuad_media=eprima_cuad_media*(N-1)/N+eprima_cuad/N;
       N=N+1;
     }
    std::cout<<"energia media manual  "<<eprimamedia<<"\n";
    var=eprima_cuad_media-pow(eprimamedia,2);
    std::cout<<"Var="<<var<<"\n";
    cv=var*Npart/pow(Temp,2);
    std::cout<<"Cv="<<cv<<"\n";
  





  
   return 0;
}
	

      


















	 
