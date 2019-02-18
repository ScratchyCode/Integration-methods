# Intro
I seguenti codici simulano il moto di un pendolo forzato smorzato con parametri prestabiliti.
È possibile plottare la traiettoria, lo spazio delle fasi (anche in 3D con un asse temporale), e grafici per l'analisi del caos deterministico come i bacini di attrazione, le sezioni di Poincare e i diagrammi di biforcazione. Nella cartella "Plot" si possono trovare degli output delle elaborazioni svolte.

Attenzione: molte di queste elaborazioni (come il calcolo dei bacini di attrazione) sono molto onerose computazionalmente.
Leggere attentamente le note.

# Metodi
Le integrazioni numeriche possibili sono eseguite attraverso i metodi:

1) Eulero
2) Eulero-Cromer
3) Verlet (autoavviante)
4) Rouge-Kutta 4

# How to
Per avviare il programma lanciare lo script "start.sh".

Per cambiare i parametri del moto o il metodo di integrazione modificare il file "condizioni_iniziali.dat".

Per cambiare l'equazione differenziale da integrare modificare le macro nel file sorgente  "integratore.c".

# Note
È stato necessario introdurre vari livelli di ottimizzazione dei calcoli;
vedere i commenti nel file sorgente "integratore.c" per personalizzare l'elaborazione attraverso le flag booleane e le variabili da preprocessore.

I plot vengono salvati in formato .png nella directory di esecuzione (se è definita la variabile da precompilatore "PNG").
