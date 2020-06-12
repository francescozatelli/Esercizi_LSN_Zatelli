# Esercizi LSN
Francesco Zatelli, 907321
francesco.zatelli@studenti.unimi.it

Tutte le cartelle sono dotate di makefile. In alcuni casi il flag nel makefile è imposato per utilizzare il compilatore Intel.
Gli eseguibili generati sono sempre chiamati "main" a parte quando ho utilizzato due eseguibili per lo stesso esercizio (ad esempio nell'esercizio 8.2 c'è "main_annealing" per l'ottimizzazione dei parametri variazionali e "main_sampling" per il calcolo dell'energia della funzione d'onda di prova).

Segnalo che nell'esercizio 10.2 e nelle esercitazioni con TensorFlow ho dovuto apportare delle modifiche rispetto alle indicazioni date a lezione, che immagino possano dare problemi nel caso provaste a far andare il codice con versioni e/o configurazioni diversi di OpenMPI e TensorFlow:
* Nell'esercizio 10.2 non mi funzionavano i tipi di dati indicati a lezione, ad esempio MPI_INT, ho dovuto usare MPI_INTEGER.
* Per utilizzare TensorFlow ho dovuto apportare le modifiche discusse a lezione agli import. Inoltre, ho dovuto settare il flag "KMP_DUPLICATE_LIB_OK" a true (a inizio notebook). Non so esattamente cosa faccia, ho trovato la soluzione online e funziona.