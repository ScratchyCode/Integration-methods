I seguenti codici simulano il moto di un pendolo forzato smorzato con parametri prestabiliti.
È possibile plottare la traiettoria, lo spazio delle fasi (anche in 3D con un asse temporale),
e grafici per l'analisi del caos deterministico, come i bacini di attrazione,
le sezioni di Poincare e i diagrammi di biforcazione.

Per avviare il programma lanciare lo script "start.sh"
Per cambiare le condizioni iniziali ed i parametri modificare il file "condizioni_iniziali.dat"
Per cambiare l'equazione differenziale da integrare modificare la macro nel file sorgente "integratore.c"

I plot vengono salvati in formato .png nella directory di esecuzione
(se è definita la variabile da precompilatore "PNG")
