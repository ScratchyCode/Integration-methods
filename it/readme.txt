I seguenti codici simulano il moto di un pendolo forzato smorzato con parametri prestabiliti.
È possibile plottare la traiettoria, lo spazio delle fasi (anche in 3D con un asse temporale),
e grafici per l'analisi del caos deterministico, come i bacini di attrazione,
le sezioni di Poincare e i diagrammi di biforcazione.

Per avviare il programma lanciare lo script "start.sh"

Per cambiare i parametri del moto modificare il file "condizioni_iniziali.dat"

Per cambiare l'equazione differenziale da integrare modificare la macro nel file .c

È stata necessaria l'introduzione di vari livelli di ottimizzazione dei calcoli;
per questo vedere i commenti nel file sorgente .c

I plot vengono salvati in formato .png nella directory di esecuzione
(se è definita la variabile da precompilatore "PNG")
