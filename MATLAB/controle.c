//-------------------------------------------------------------------------//
//                Federal University of Santa Maria - UFSM                 //
//          Power Electronics and Control Research Group - GEPOC           //
//                 97105-900 - Santa Maria, RS - BRAZIL                    //
//                       http://www.ufsm.br/gepoc                          //
//                                                                         //
// Desenvolvido por: Caio Ruviaro Dantas Osorio                            //
//=========================================================================//
//                         Descrição do Programa                           //
//=========================================================================//
//
//-------------------------------------------------------------------------//

#include <math.h>
#include <stdlib.h>
#include "controle.h"

#define pontosig 167 // Esse valor vai ser gerado quando fizer a senoide no .m, tem que mudar no vetor tbm
#define pi 3.14159265

// Ganhos do Ressonante fs=10kHz
#define Sd11 0.03125
#define Sd21 0
#define Rd11 1.998579261776922
#define Rd21 1
#define Rd12 -0.999999980007103
#define Rd22 0

// Ganhos da realimentação

// Deadbeat Robusto para r=0.67  variac 10% fs=10kHz ----------------------

/*
#define k1 -129.5821
#define k2 -1.6971
#define k3 1714.2756
#define k4 -1458.2573
*/

// Deadbeat Robusto para r=0.92  fs=10kHz ----------------------

/*
#define k1 -43.077497268237991
#define k2 -1.161581703009778
#define k3 201.6441834456933
#define k4 -192.7823981606899
*/

// Deadbeat Robusto para r=0.95 e fs=10kHz----------------------
/**/
#define k1 -22.626994881450905
#define k2 -0.561421978901323
#define k3 79.163345802203452
#define k4 -76.791082664183065


// Deadbeat Robusto para r=0.95 e fs=20kHz----------------------
/*
#define k1 -41.6226785465594
#define k2 -0.5159384350689
#define k3  269.7946254981793
#define k4 -261.1675197016376
*/

// dlqr -----------------------------------

//#define k1  -9.648000338852821
//#define k2 -0.177855230132305
//#define k3 25.526886588975405
//#define k4 -24.736778873466072

// Deadbeat Convencional ------------------
/*
#define k1 -268.331177193895
#define k2 -2.786579261777
#define k3 5406.371471524585
#define k4 -4136.987822288876
*/

/*
#define k1 -299.245259442553
#define k2 -2.996579261777
#define k3 6377.287689278131
#define k4 -4790.910441172113
*/

/*
#define k1 -2.3236
#define k2 0
#define k3 0.8943
#define k4 -0.9592
*/


// Variáveis
static double t_ks = 0;	                  		// Constantes
static int n = 0;
static int pontos_am = 0;
static double modul_ma_k = 0, modul_ma_nxt1 = 0, v_abf_k = 0;
static double ig_k = 0, ig_ref_ks=0, ig_ref_k=0, Vcc = 200;
static double ig_ref_nxt2 = 0, ig_cont_nxt1 = 0, ig_ref_amost_ks = 0;
static double v_tri = 0, Sa = 0, Sb = 0;

static double u_ks=0, x_ks=0, theta_ks=0, rho1_ks=0, rho2_ks=0;

static double Rf = 0.1, Lf = 0.005, Ig_pk = 10, theta = 0, w = 376.99;
static double u_k=0, u_nxt1=0, cont_aux=1, TPER=0;

static double teste=0,upwm=0,cont=1, rho1_ks_nxt=0, rho2_ks_nxt=0;

// Definicoes de frequencia
//No typhoon estamos utilizando dT=5e-7
static double Fs = 10000, Ts = 1.0e-4;		// Freq. e Periodo de amostragem
static double Fsw = 10000, Tsw = 1.0e-4;    // Freq. e Periodo de comutação    

#define Pontos_fs 200         // fs =10kHz - Pontos_fs = 1/(dT*fs)  
#define Pontos_fsw 200        // fsw=10kHz - Pontos_fsw = 1/(dT*fsw)  

static int t_sample = 200;     // Deve ser inicializada com o valor de Pontos_fs para entrar a primeira vez no loop. 


//=========================================================================//
// Programa principal
//=========================================================================//
//, double t, double delt
double control(double tensao, double corrente) 
{
    int t = 1;
    
    ig_ref_k = Ig_pk*sin(w*t);
     
	// Controlador
	if(t_sample==Pontos_fs) // Amostragem Ts (aqui vai o Pontos_fs do matlab, ou seja Ts/dT) -> Cuidar para inicializar essa variável já com esse valor
	{
		ig_k = corrente;
		v_abf_k = tensao;
		
		rho1_ks = rho1_ks_nxt;
		rho2_ks = rho2_ks_nxt;
		
        ig_ref_ks = ig_ref_k;
        //ig_ref_ks = Ig_pk*sin(w*t_ks);
        
        modul_ma_k = modul_ma_nxt1;
        u_k = u_nxt1;
       
        //u_ks = k1*x_ks + k2*theta_ks + k3*rho1_ks + k4*rho2_ks;
        modul_ma_nxt1 = k1*ig_k + k2*modul_ma_k + k3*rho1_ks + k4*rho2_ks;
        u_nxt1 = (modul_ma_nxt1/(2*Vcc))+0.5;
        
        rho1_ks_nxt = -Sd11*ig_k + 0*modul_ma_k + Rd11*rho1_ks + Rd12*rho2_ks + 0*modul_ma_nxt1 + Sd11 * ig_ref_ks;
        rho2_ks_nxt = -Sd21*ig_k + 0*modul_ma_k + Rd21*rho1_ks + Rd22*rho2_ks + 0*modul_ma_nxt1 + Sd21 * ig_ref_ks;
        
        x_ks = ig_k;

        t_ks = t_ks + Ts;      
		t_sample = 0;		
	}

    //Modulador DSP

    TPER = Pontos_fsw/2;   // PWM simétrico (triangular)
       
    if(cont_aux<Pontos_fsw){
       cont_aux=cont_aux+1;}
    else
       cont_aux=1;    
       
    if(cont_aux<TPER){
       cont=cont+1;}
    else if(cont_aux>=TPER && cont_aux<Pontos_fsw){
         cont=cont-1;}
    else 
       cont=1;
       
        
    if(u_k*TPER>=cont){
       Sa=1;}
    else
       Sa=0;
       
   // Sb=1-Sa;

    t_sample = t_sample + 1;
	
	// PWM
	//out[0]=modul_ma_k;

	// Sinais
	//out[1]=u_k;
	//out[2]=ig_k;
	//out[3]=cont;
	//out[4]=ig_ref_ks;
	//out[5]=ig_ref_k;
    
    return Sa;

}











