// Scritto da Pietro Squilla e Francesco Valerio Servilio
// Simula il moto di un pendolo generale (smorzato forzato)
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

// usiamo le macro invece di funzioni normali per velocizzare la chiamata

/*** equazioni differenziali definite con le macro ***/
#define F1(sp,t) (sp.v)
#define F2(sp,t) (-w2*sin(sp.x)-gammma*sp.v+fe*cos(we*t))
/*** equazioni differenziali definite con le macro ***/

// #define NOBACINI 1           // per evitare di fare i calcoli per i bacini di attrazione (molto onerosi)
#define RATEO 137./150.         // rapporto tra le x e le y per produrre un'immagine quadrata
// #define VARIABILI 1          //se definita, alcuni parametri vengono impostati "hard-coded"
// #define PI 3.1415926535898   // definisco un pi greco con un dato numero di cifre decimali così da usare lo stesso nelle condizioni iniziali
#define PI M_PI  
#define PNG 1                   // per plottare su png (maggiore risoluzione)
#define PROG 1                  // scrive su terminale la percentuale di avanzamento del processo (rallenta un po')               
#define pSize 1200.0            // pixel sull'asse x dei file .png (la quantita di pixel sull'asse y viene calcolata automaticamente per rendere l'immagine quadrata)

/******************** strutture dati ********************/
typedef struct spazioDelleFasi{
    double x;
    double v;
}spazioFasi;

struct dEdt{
    double dt;
    double dE;
};
/******************** strutture dati ********************/


/******************** variabili globali ********************/
int indice=0, run=0;
double x0=0., v0=0., w2=0., tmax=0., dt=0., E0=0., E=0., gammma=0., fe=0., we=0.;
spazioFasi (*metodo)(spazioFasi ps, double t);
/******************** variabili globali ********************/


/******************** prototipi funzioni ********************/
spazioFasi eulero(spazioFasi ps, double t);
spazioFasi eulero_cromer(spazioFasi ps, double t);
spazioFasi auto_verlet(spazioFasi ps, double t);
spazioFasi rk4(spazioFasi ps, double t);
void trova_radici(double tempX, double x, double tempT, double t, double *radici);
void init_pipe(FILE *pipe, char *xlabel, char *ylabel, char *titolo, char *nome_grafico, char *plotCmd);
void plot_solAnalitiche(FILE *pipe);
void plot_traiettoria(double t, double x, FILE *pipe);
void plot_spazioFasi3D(double x, double v, double t, FILE *pipe);
void plot_spazioFasi(double x, double v, FILE *pipe);
void plot_energia(double t, FILE *pipe);
void plot_dE_dt(struct dEdt *lin_dE_dt, FILE *pipe);
void plot_periodi(double *radici, double *periodi, FILE *pipe);
void plot_poincare(double x, double v, FILE *pipe);
void plot_bacini(double v, FILE *pipe);
void plot_biforc(double fe, double v, FILE *pipe);
void controllo_gnuplot(void);
void controllo_metodo(int algo);
void chiudiPipe(FILE *pipe);
spazioFasi prodPS(spazioFasi a, double k);
spazioFasi sumPS(spazioFasi a, spazioFasi b);
void checkPtr(void *ptr);
/******************** prototipi funzioni ********************/


/******************** main ********************/
int main(int argc, char **argv){
    int algo, analisiCaos;
    double DT, DX, DV, DF;
    FILE *input = fopen("condizioni_iniziali.dat","r");
    checkPtr(input);
    
    // controllo l'esistenza di gnuplot sulla macchina
    controllo_gnuplot();
    
    // scansione e chiusura dello stream di input
    fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf\n",&x0,&v0,&w2,&tmax,&dt,&gammma,&fe,&we,&run,&algo,&analisiCaos,&DT,&DX,&DV,&DF);
    fclose(input);
    
    #ifdef VARIABILI
    we = 2./3.;
    x0 = M_PI/2.;
    v0 = 8.*M_PI/10.;
    dt = 2.*M_PI/(we*100.);
    #endif   
    
    // controllo esistenza metodo
    controllo_metodo(algo);
    
    // scelta del metodo di integrazione
    switch (algo){
        case 1:
            
            printf("Metodo scelto: Eulero\n");
            metodo = eulero;
            
            break;
            
        case 2:
            
            printf("Metodo scelto: Eulero-Cromer\n");
            metodo = eulero_cromer;
            
            break;
            
        case 3:
            
            if (gammma!=0.) {
                printf("Non è possibile usare Auto-Verlet se l'accelerazione dipende dalla velocità\n");
                exit(0);
            }
            
            printf("Metodo scelto: Auto-Verlet\n");
            metodo = auto_verlet;
            
            break;
            
        case 4:
            
            printf("Metodo scelto: Runge-Kutta 4\n");
            metodo = rk4;
            
            break;
    }
    
    // allocazione array
    double *radici = malloc((int)((tmax/dt) + 1.)*sizeof(double));
    double *periodi = malloc((int)(((tmax/dt)/2.) + 1.)*sizeof(double));
    struct dEdt *lin_dE_dt = (struct dEdt *)(malloc(run * sizeof(struct dEdt))); // array di strutture dinamico (meno di 2 MB in RAM per run < 100000)
    checkPtr(radici);
    checkPtr(periodi);
    checkPtr(lin_dE_dt);
    
    // variabili per il calcolo
    spazioFasi ps;
    ps.x = x0;
    ps.v = v0;
    double t=0., tempX=0., tempT=0., v00=v0, x00 =x0, dt0=dt, fe0=fe, progress;
    double T; //periodo della forzante
    int i, j=0, k, l=0;
    bool flag_giro1 = true; // ottimizzazione while
    bool flag_periodi = false; // plotta periodi quando validi
    bool flag_dedt = false; // flag per fare il plot di dE vs dt
    bool flag_bacini = false; // flag per i bacini di attrazione
    bool flag_biforc = false; // flag per i diagrammi di biforcazione
    bool flag_soluzan = false; // flag per plottare le soluzioni analitiche
    
    if(analisiCaos) progress = 2.*run+pow(run,2);
    else progress = run;
    
    if(fe == 0.) flag_periodi = true; // se non c'e' una forzante --> plotta periodi
    
    if(gammma == 0. && fe == 0.){ // le soluzioni analitiche e la divergenza dell'energia all'aumentare del passo funzionano solo per pendolo non smorzato/forzato
        flag_soluzan = true;
        if(run > 1) flag_dedt = true;
    } else if(run > 1 && analisiCaos) flag_bacini = true;
    
    if(we != 0) T = 2.*PI/we;
    else T = 0;
    
    E0 = pow(v0,2) + w2*pow(x0,2);
    
    // apro un'istanza di gnuplot per ogni grafico attraverso le pipe
    FILE *TRAIETT = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *SP3D = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *SP = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *ENERGIA = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *SOLUZAN = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *DEDT = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *PERIODI = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *POINCARE = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *BACINI = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    FILE *BIFORC = popen("/usr/bin/gnuplot -persist >/dev/null 2>/dev/null","w");
    
    checkPtr(TRAIETT);
    checkPtr(SP3D);
    checkPtr(SP);
    checkPtr(ENERGIA);
    checkPtr(SOLUZAN);
    checkPtr(DEDT);
    checkPtr(PERIODI);
    checkPtr(POINCARE);
    checkPtr(BACINI);
    checkPtr(BIFORC);
    
    // inizializzo tutti gli stream delle pipe
    init_pipe(TRAIETT,"tempo","spazio (angolo)","traiettoria","traiettoria.png","");
    
    init_pipe(SP3D,"spazio (angolo)","velocita","spazio delle fasi vs tempo","spazio_fasi_3d.png","splot '-' u 1:2:3 w l");
    
    init_pipe(SP,"spazio (angolo)","velocita","spazio delle fasi","spazio_fasi.png","plot '-' u 1:2 pt 7 ps .2");
    
    init_pipe(ENERGIA,"tempo","energia","conservazione energia","conservazione_energia.png","");
    
    char buf[120] = "plot '-' u 1:2:3 with points pt 7 ps 0.4 palette ";
    #ifdef PNG
    sprintf(buf,"set terminal png size %d, %d \n plot '-' u 1:2:3 with points pointtype 7 ps 0.1 palette ",(int)(run*1.44) + 1,(int)(run*1.44*RATEO) + 1); // 1.44 è il rapporto giusto
    #endif
    init_pipe(BACINI,"x_0","v_{0}","bacini di attrazione (colore in base a v(t=t_{max}))","bacini.png",buf);
    
    init_pipe(BIFORC,"f_e","v(t_{max})","diagramma di biforcazione (velocita al tempo finale in funzione del modulo della forzante esterna)","biforcazione.png","plot '-' u 1:2 pt 7 ps 0.1 lc rgb 'black' ");
    
    if(tmax>T) init_pipe(POINCARE,"spazio (angolo)","velocita(n*periodo forzante)","sezione di Poincare di velocita e spazio","poincare.png","plot '-' u 1:2 pt 7 ps .2 lc rgb 'black' ");
    
    if(flag_soluzan) init_pipe(SOLUZAN,"tempo","spazio (angolo)","traiettoria analitica pendolo semplice (piccole oscillazioni)","traiettoria_analitica.png","plot '-' u 1:2 w l lc rgb 'red' ");
    
    if(flag_periodi) init_pipe(PERIODI,"n","periodo pendolo","periodi","periodi.png","");
    
    if(flag_dedt) init_pipe(DEDT,"log_{10}(dt)","log_{10}(dE)","dE(t_{max}) vs dt","dEdt.png","");
    
    
    for(i=1; i<=run; i++) { 
        
        #ifdef PROG
        j++;
        printf("\rCalcolo in corso:\t\t%2d%%",(int)(j*100./progress));
        #endif
        
        t = 0.;
        ps.x = x0;
        ps.v = v0;
        E = E0;
        
        for(k=1; t < tmax; k++){ // calcola i valori del moto
            if(flag_giro1){ // eseguito solo al primo giro del for
                plot_traiettoria(t,ps.x,TRAIETT);
                plot_spazioFasi3D(ps.x,ps.v,t,SP3D);
                plot_spazioFasi(ps.x,ps.v,SP);
                plot_energia(t,ENERGIA);
                if(tmax > T && fmod(t,T) < dt) plot_poincare(ps.x,ps.v,POINCARE); //sezioni di poincare
                if(flag_periodi){
                    trova_radici(tempX,ps.x,tempT,t,radici);
                    tempX = ps.x;
                    tempT = t;
                }
            }
            
            /********integrazione********/
            ps = metodo(ps,t); // ritorna lo spazio delle fasi aggiornato
            //t += dt;
            /********integrazione********/
            t = dt*k;
            
            if(flag_dedt || flag_giro1) E = pow(ps.v,2) + w2*pow(ps.x,2);
        }
        
        flag_giro1 = false; // ottimizzazione su esecuzione codice non dipendende dal for
        
        #ifndef NOBACINI
        if(flag_bacini){ // bacini di attrazione (per primi poichè vanno eseguiti run*run volte per evitare branch mispredictions)
            plot_bacini(ps.v,BACINI); // plotta il colore in base a v(tmax) e sugli assi x0 e v0
            x0 = x00 + DX*i;
            if(i == run){
                x0 = x00;
                l++;
                v0 = v00 + DV*l;
                if(v0 >= v00 + DV*run){
                    v0 = v00;
                    flag_bacini = false;
                    flag_biforc = true;
                }
                i = 1;
            }
            continue; // ottimizzazione del ciclo for
        }
        #else
        flag_biforc = true;
        #endif
        
        if(flag_dedt){ // scrittura sul vettore di strutture dE_dt
            lin_dE_dt[i].dt = log10(dt);
            lin_dE_dt[i].dE = log10(E - E0);
            dt = dt0 + DT*i;
            if(analisiCaos && i == run){ // se bisogna fare l'analisi del caos fa ripartire il while
                dt = dt0;
                i = 1;
                flag_dedt = false;
                flag_bacini = true;
            }
            continue;
        }
        
        if(flag_biforc){ // diagrammi di biforcazione
            plot_biforc(fe,ps.v,BIFORC);
            fe = fe0 + DF*i;
        }
        
    }
    
    // plotta le soluzioni analitiche
    if(flag_soluzan) plot_solAnalitiche(SOLUZAN);
    
    // plotto dE vs dt
    if(run > 1 && (gammma == 0. && fe == 0.)) plot_dE_dt(lin_dE_dt,DEDT);
    
    // trova e plotta i periodi
    if(flag_periodi) plot_periodi(radici,periodi,PERIODI);
    
    // libero la memoria allocata
    free(radici);
    free(periodi);
    free(lin_dE_dt);
    
    // pulisco e chiudo gli stream
    chiudiPipe(SP3D);
    chiudiPipe(SP);
    chiudiPipe(ENERGIA);
    chiudiPipe(TRAIETT);
    chiudiPipe(SOLUZAN);
    chiudiPipe(DEDT);
    chiudiPipe(PERIODI);
    chiudiPipe(POINCARE);
    chiudiPipe(BACINI);
    chiudiPipe(BIFORC);
    
    printf("\rOperazione completata:\t\t100%%");
    
    return 0;
}  
/******************** main ********************/


/******************** definizioni funzioni ********************/
spazioFasi eulero(spazioFasi ps, double t){
    spazioFasi tempPS = ps;    // tempPS e' ps(t=tn+dt), ps e' ps(t=tn)
    
    tempPS.x += F1(ps,t)*dt;
    tempPS.v += F2(ps,t)*dt;
    
    return tempPS;
}

spazioFasi eulero_cromer(spazioFasi ps, double t){
    spazioFasi tempPS = ps;
    
    tempPS.v += F2(ps,t)*dt;
    tempPS.x += F1(tempPS,t + dt)*dt;
    
    return tempPS;
}

spazioFasi auto_verlet(spazioFasi ps, double t){
    spazioFasi tempPS=ps;
    
    tempPS.x += F1(ps,t)*dt + F2(ps,t)*pow(dt,2)*0.5;
    tempPS.v += (F2(ps,t) + F2(tempPS,t + dt))*dt*0.5;
    
    return tempPS;
}

spazioFasi rk4(spazioFasi ps, double t){
    spazioFasi k[5];
    
    k[1].x = dt * F1(ps,t); 
    k[1].v = dt * F2(ps,t);
    k[2].x = dt * F1(sumPS(ps,prodPS(k[1],0.5)),t+dt/2.);
    k[2].v = dt * F2(sumPS(ps,prodPS(k[1],0.5)),t+dt/2.);
    k[3].x = dt * F1(sumPS(ps,prodPS(k[2],0.5)),t+dt/2.);
    k[3].v = dt * F2(sumPS(ps,prodPS(k[2],0.5)),t+dt/2.);
    k[4].x = dt * F1(sumPS(ps,k[3]),t + dt); // y4
    k[4].v = dt * F2(sumPS(ps,k[3]),t + dt); // y4
    
    return sumPS(ps,prodPS(sumPS(sumPS(k[1],k[4]),prodPS(sumPS(k[2],k[3]),2.)),1./6.));
}

spazioFasi prodPS(spazioFasi a, double k){
    
    a.x *= k;
    a.v *= k;
    
    return a;
}

spazioFasi sumPS(spazioFasi a, spazioFasi b){
    
    a.x += b.x;
    a.v += b.v;
    
    return a;
}

// ritorna il numero di radici trovate
void trova_radici(double tempX, double x, double tempT, double t, double *radici){
    
    if((tempX-x0)*(x-x0) < 0){
        radici[indice] = (t + tempT)/2.;
        indice++;
    }
    
    return;
}

void init_pipe(FILE *pipe, char *xlabel, char *ylabel, char *titolo, char *nome_grafico, char *plotCmd){
    
    fprintf(pipe,"set autoscale \n");
    
    #ifdef PNG
    fprintf(pipe,"set terminal png size %d,%d\n",(int)pSize,(int)(pSize*RATEO)); // imposta la risoluzione dell'immagine in base alla dimensione della griglia di incrementi
    //fprintf(pipe,"set terminal png size 1920,1080\n");
    //fprintf(pipe,"set terminal png size 2560,2048\n");
    #else
    fprintf(pipe,"set terminal x11 \n");
    #endif  
    
    fprintf(pipe,"set xlabel '%s' \n",xlabel);
    fprintf(pipe,"set ylabel '%s' \n",ylabel);
    fprintf(pipe,"set title '%s' \n",titolo);
    
    #ifdef PNG
    fprintf(pipe,"set output '%s' \n",nome_grafico);
    #endif  
    
    if(!strcmp(plotCmd,"")) fprintf(pipe,"plot '-' u 1:2 w l notitle \n");
    else fprintf(pipe,"%s notitle \n", plotCmd);
    
    return;
}

void plot_solAnalitiche(FILE *pipe){
    double t=0, x;
    double A = sqrt(pow(x0,2.) + (pow(v0,2.)/w2));
    double fi = atan(-((v0)/(x0*sqrt(w2))));
    //double periodo = (2*M_PI)/w;
    
    while(t <= tmax){
        x = A*cos(sqrt(w2)*t + fi);
        fprintf(pipe,"%.15lf %.15lf\n",t,x);
        t += dt;
    }
    
    return;
}

void plot_poincare(double x, double v, FILE *pipe){
    
    //restringo x tra -PI e +PI
    if(x > M_PI) x = fmod(x + M_PI,2*M_PI) - M_PI;
    else if(x < M_PI) x = fmod(x - M_PI,2*M_PI) + M_PI;
    
    fprintf(pipe,"%.15lf %.15lf\n",x,v);
    
    return;
}

void plot_bacini(double v, FILE *pipe){
    
    fprintf(pipe,"%.15lf %.15lf %.15lf\n",x0,v0, v);
    
    return;
}

void plot_biforc(double fe, double v, FILE *pipe){
    
    fprintf(pipe,"%.15lf %.15lf\n",fe,v);
    
    return;
}

void plot_traiettoria(double t, double x, FILE *pipe){
    
    fprintf(pipe,"%.15lf %.15lf\n",t,x);
    
    return;
}

void plot_spazioFasi3D(double x, double v, double t, FILE *pipe){
    
    fprintf(pipe,"%.15lf %.15lf %.15lf\n",x,v,t);
    
    return;
}

void plot_spazioFasi(double x, double v, FILE *pipe){
    
    //restringo x tra -PI e +PI
    if(x>M_PI) x=fmod(x+M_PI,2*M_PI)-M_PI;
    else if(x<M_PI) x=fmod(x-M_PI,2*M_PI)+M_PI;
    
    fprintf(pipe,"%.15lf %.15lf \n",x,v);
    
    return;
}

void plot_energia(double t, FILE *pipe){
    
    fprintf(pipe,"%.15lf %.15lf\n",t,E);
    
    return;
}

void plot_dE_dt(struct dEdt *lin_dE_dt, FILE *pipe){
    int i;
    
    for(i=0; i<run; i++){
        fprintf(pipe,"%.15lf %.15lf\n",lin_dE_dt[i].dt,lin_dE_dt[i].dE);
    }
    
    return;
}

void plot_periodi(double *radici, double *periodi, FILE *pipe){
    int i, j=0;
    
    // per trovare i periodi
    for(i=1; i<indice/2.-1; i++){
        periodi[i-1] = (radici[2*i+1] - radici[2*i-1]);
        j++;
    }
    
    // grafico i periodi
    for(i=0; i<j; i++){
        fprintf(pipe,"%d %.15lf\n",i+1,periodi[i]);
    }
    
    return;
}

void controllo_gnuplot(void){
    FILE *output = fopen("/usr/bin/gnuplot","r");
    
    if(output == NULL){
        printf("\nInstalla gnuplot con: sudo apt install gnuplot\n\n");
        exit(1);
    }
    
    return;
}

void controllo_metodo(int algo){
    
    if(algo <= 0 || algo >= 5){
        printf("\nMetodo non trovato!\n\n");
        exit(1);
    }
    
    return;
}

void chiudiPipe(FILE *pipe){
    
    fflush(pipe);
    fprintf(pipe,"quit\n");
    pclose(pipe);
    
    return;
}

void checkPtr(void *ptr){
    
    if (ptr==NULL){
        perror("\nERROR");
        fprintf(stderr,"\n");
        exit(0);
    }
    
}
/******************** definizioni funzioni ********************/
