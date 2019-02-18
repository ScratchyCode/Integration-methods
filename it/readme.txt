#Intro
I seguenti codici simulano il moto di un pendolo forzato smorzato con parametri prestabiliti.
È possibile plottare la traiettoria, lo spazio delle fasi (anche in 3D con un asse temporale),
e grafici per l'analisi del caos deterministico, come i bacini di attrazione,
le sezioni di Poincare e i diagrammi di biforcazione.

#Metodi
Le integrazioni numeriche possibili sono eseguite attraverso i metodi:
1) Eulero
2) Eulero-Cromer
3) Verlet (autoavviante)
4) Rouge-Kutta 4

#How to
Per avviare il programma lanciare lo script "start.sh"

Per cambiare i parametri del moto o metodo modificare il file "condizioni_iniziali.dat"

Per cambiare l'equazione differenziale da integrare modificare la macro nel file .c

#Nota
È stata necessaria l'introduzione di vari livelli di ottimizzazione dei calcoli;
per questo vedere i commenti nel file sorgente .c

I plot vengono salvati in formato .png nella directory di esecuzione
(se è definita la variabile da precompilatore "PNG")
