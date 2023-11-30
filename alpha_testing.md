# fdaPDE v2.0 - alpha testing

**Status del testing** : potete controllare lo stato di tutte le issues aperte e bug noti qui: [dashboard](https://github.com/orgs/fdaPDE/projects/2/views/1)

> Nel seguito indicherò con `fdaPDE-CRAN` la versione 1.1-16 della libreria attualmente disponibile su CRAN ([link](https://cran.r-project.org/web/packages/fdaPDE/index.html)). `fdaPDE-2.0` si riferisce alla nuova versione cha andrà a sostituire quella attulamente ufficiale.

Questa pagina riporta le istruzioni da seguire durante la procedura di alpha-testing di `fdaPDE-2.0`. Contiene inoltre una prima documentazione **non ufficiale** (ed informale) per approcciare la libreria e la sua nuova interfaccia.

Questa pagina sarà inoltre aggiornata periodicamente per tutta la fase di testing, per tenere traccia dello sviluppo in corso della libreria. Controllate il [changelog](#changelog) per gli utlimi aggiornamenti. 

> `fdaPDE-2.0` è ancora in attivo sviluppo, reinstallate il pacchetto spesso. Inoltre ignorate evenutali warnings in installazione. Cambiamenti di interfaccia potrebbero verificarsi per tutto il periodo di testing (virtualmente, fin quando il pacchetto non arriva su CRAN)

* [Obbiettivo](#Obbiettivo)
* [Testing](#testing)
* [Bug reports](#Bug-reports)
* [Interfaccia](#Interfaccia)
  * [Functional Space](#Functional-space)
  * [PDE - esprimere penalty generiche](#PDEs)
  * [SRPDE - spatial regression](#srpde)
  * [STRPDE -  spatio-temporal regression](#strpde)
  * [Decorare il fitting di un modello](#decorare-il-fitting-di-un-modello)
  * [GCV - Generalized Cross Validation](#gcv---generalized-cross-validation)
  
## ChangeLog

* 30/11/23: inizio della fase di alpha-testing. Prima versione usabile dell'interfaccia `R`.

## Obbiettivo
Lo scopo di questa prima fase di testing è quella di verificare la stabilità della libreria, di suggerire eventuali migliorie, e di arricchire la suite di testing (sia lato R, ma soprattutto C++).
* *correttezza*: questo è il momento per fare check intesivi di correttezza numerica. Rispetto alla versione `fdaPDE-CRAN`, le differenze in norma $L^\infty$ possono essere nell'ordine di 10-7 (quindi non aspettatevi differenze zero macchina). Differenze più alte, possono esserci, e possono essere dovute a
  * bug presenti nella versione CRAN che sono stati corretti nella versione 2.0
  * valori di default differenti per gli algoritmi interni
  
  Se incontrate differenze superiori a 10^-7 notificatemelo (vedi [bug reports](#bug-reports)). Vi chiedo di testare estensivamente soprattutto geometrie diverse dal 2D (test su 2D sono ben accetti, ma l'aspetto 3D e manifold non è stato testato formalmente ed in maniera così estensiva come il 2D).
* *stabilità*: usate il pacchetto, e fatelo crashare. Insomma, usate male il pacchetto, date input sbagliati, fornite dati volutamente sbagliati. Stressate il pacchetto dandogli dataset di grandi dimensioni, mesh molto grandi, etc. Date valori per gli argomenti a funzioni che possono essere in conflitto, fate cose senza senso, etc. Dovete cercare quei casi limite che la libreria deve riuscire a gestire, ma senza crashare (sperabilmente ritornando un messaggio di errore). Il fatto che il codice non crashi a seguito di azioni sbagliate, è comunque da segnalare.
   
   Sarebbe anche interessante testare l'installabilità del pacchetto `R` su sistemi operativi diversi (`windows`, `mac-os`, `linux`, dove con `linux` intendo sia sistemi debian-based, come `ubuntu`, ma anche `fedora`, `openSUSE`, `archlinux`, etc.), e magari con compiler diversi (principalmente `gcc`/`g++`, `clang`, `MVSC`). In caso di problemi di installazione (il pacchetto non è mai stato installato su sistemi `mac-os` ad esempio), notificate.
* *usabilità*: notificate se mancano dei parametri lato `R` che sono presenti in `fdaPDE-CRAN` ma non risultano accessibili in `fdaPDE-2.0`. Se la libreria ritorna messaggi di errore troppo criptici, notificate (magari suggerite quale potrebbe essere un messaggio più esplicativo). Segnalate inoltre se i parametri dei metodi, o le funzioni stesse, hanno nomi poco intuitivi.
* *migloramenti*: se usando il pacchetto vi vengono in mente delle idee che potrebbero migliorare la sua interfaccia, e che magari avete visto usare in altri pacchetti, notificate. Vale anche per features che vorreste avere nel pacchetto (magari di interesse per specifici settori che state trattando ora, o avete trattato in passato). In questo caso le idee saranno prima discusse e poi inserite (allo stesso modo proposte per eventuali cambiamenti di interfaccia). Segnalate inoltre se il codice rende inaccessibili alcune proprietà di alcuni oggetti, o se ricavare quella proprietà è troppo complesso.

## Testing

Una volta che avete testato una funzionalità, vi chiedo di scrivere uno script di test ben isolato (MWE: Minimal Working Example) che riporta la generazione dei dati, così come il test sul modello vero e proprio. Il test non deve contenere nessun plot o output su console. 

> Sarebbe l'ideale accettare una certa funzionalità come funzionante se è dichiarata tale da almeno 2 persone.

<!--Il test deve avvalersi delle funzionalità offerte dal pacchetto R [testthat](https://testthat.r-lib.org/index.html).-->

<!--Una parte di questi test (quelli più significativi), saranno inseriti nel pacchetto `R` ufficiale.--> 
Lo scopo primario è quello di generare un equivalente test C++, che fossilizza la correttezza numerica dei metodi, in ogni sfumatura possibile.

## Bug reports

Se incontrate un bug/crash, notificatemelo. Vi chiedo di 
* spiegare con due parole cosa succede (crash / problema numerico / etc.) 
* è **fondamentale** che voi mi diate uno script R minimale che genera il problema. In generale, non posso fixare un bug se non riesco a riprodurlo.

## Interfaccia
Questa sezione contiene una presentazione non-ufficiale dell'interfaccia attualmente disponibile, in particolare, non vuole essere una documentazione dettagliata. Inoltre sottolineo che della documentazione ufficiale non è ancora disponibile.
### Functional space

Rapprsenta il concetto di base funzionale su un certo dominio. 
Al momento, supporta solo basi FEM (per la fine di questo alpha-testing, supporto anche per basi B-splines).
```R
library(fdaPDE2)
data("unit_square", package = "fdaPDE2")
unit_square <- Mesh(unit_square)

## the functional space of finite element functions of order 1, over the unit square
Vh <- FunctionSpace(unit_square, fe_order = 1)
basis <- Vh$get_basis() ## recover the (lagrangian) basis system of Vh

## evaluate the basis system on a given set of locations (compute Psi matrix)
locations <- unit_square$nodes
Psi <- basis$eval(type = "pointwise", locations)

## integrate a function (expressed as basis expansion on Vh) over the domain
f <- function(p) { p[,1]^2 + p[,2]^2 } ## x^2 + y^2
Vh$integrate(f)
```

### PDEs

`fdaPDE-2.0` permette la scrittura esplicita di PDEs in forma forte. Nella definizione di un modello statistico, questa procedura non è necessaria nel caso più comune di regolarizzazione tramite Laplaciano semplice.
```R
library(fdaPDE2)
data("unit_square", package = "fdaPDE2")
unit_square <- Mesh(unit_square)

## the functional space of finite element functions of order 1, over the unit square
Vh <- FunctionSpace(unit_square, fe_order = 1)
f  <- Function(Vh) ## a generic element of the space Vh

## compose the differential operator (in strong form)
Lf <- -laplace(f) + dot(c(1,1), grad(f)) ## a costant coefficients advection-diffusion problem
## define the forcing term
u <- function(points) { return(rep(1, times = nrow(points))) }

## create your penalty
penalty <- pde(Lf, u)
```

Operatori supportati
| operatore  | codice           | note                                                                                                       |
|------------|------------------|------------------------------------------------------------------------------------------------------------|
| laplaciano | `laplacian(f)`   | scrivere `laplacian(f)` è diverso da scrivere `-laplacian(f)`                                              |
| divergenza | `div(K*grad(f))` | il tensore di diffusione `K` può essere una matrice o una funzione che ritorna una matrice (space-varying) |
| trasporto  | `dot(b,grad(f))` | `b` deve essere un vettore o un campo vettoriale                                                           |
| reazione   | `c*f`            | `c` è una costante scalare o una funzione scalare                                                          |
|            | `dt(f)`          | notifica che il problema è tempo-dipendente (parabolico)                                                   |

```R
## general linear second order parabolic operator
Lf <- dt(f) - div(K * grad(f)) + dot(b, grad(f)) + c * f
```

Ulteriori esempi possono essere trovati nella documentazione di [femR](https://fdapde.github.io/femR/articles/Introduction.html) o negli [script di test](https://github.com/fdaPDE/femR/tree/stable/tests).

### SRPDE

Questa sezione mostra l'interfaccia per la definizione di problemi di regressione per dati nel solo spazio.

```R
library(fdaPDE2)
data("unit_square", package = "fdaPDE2")
unit_square <- Mesh(unit_square)

data_frame <- ## obtain your data in some way...
## currently, the only requirement is data_frame to be a data.frame object, e.g., 

##            y         x1         x2
## 1 -0.04044303  0.1402058 0.00000000
## 2  0.15079619  1.1989599 0.03447593
## 3  0.02391597 -2.3299685 0.06891086
## 4  0.38927632  0.5709451 0.10326387
## 5  0.39417457  2.7482761 0.13749409
## 6  0.33297548  1.7080400 0.17156085

## a model is first described in an "abstract" way, think for instance to the continuous
## functional J(f, \beta) we minimize when solving a smoothing problem

## a nonparametric spatial regression model
model <- SRPDE(y ~ f, domain = unit_square, data = data_frame, lambda = 1e-6)

## this will inject in the current environment a Function object named f, representing the 
## unknown spatial field and defined on a FunctionSpace(unit_square, fe_order = 1)
## the name of the spatial field supplied in formula can be any

model$fit() ## fits the model

## we could have equivalently written...
f <- Function(FunctionSpace(unit_square, fe_order = 1))
Lf <- -laplace(f) ## simple laplacian penalty
u <- function(points) { return(rep(0, times = nrow(points))) }
model <- SRPDE(y ~ f, penalty = PDE(Lf, u), data = data_frame, lambda = 1e-6)
model$fit()

## note that the last writing is more general, allowing for the definition of general penalties

## we can describe a semi-parametric spatial regression model as
model <- SRPDE(y ~ x1 + f, domain = unit_square, data = data_frame, lambda = 1e-6)
model$fit()

## plot esteimated spatial field (requires plotly)
contour(f)
```

`SRPDE` espone un parametro `family`, il cui default è settato su `gaussian` e corrisponde ad un modello di Linear Spatial Regression, come descritto ad esempio in _Sangalli, L.M. (2021),
Spatial regression with partial differential equation regularization, International Statistical Review_.
| family                                      | note                           |
|-----------------------------------------------|--------------------------------|
| uno tra `poisson`, `exponential`, `gamma`, `bernulli` | Il modello implementa una regressione spaziale generalizzata, risolto tramite applicazione di FPIRLS, come descritto in _Wilhelm, M., Sangalli, L.M. (2016), Generalized Spatial Regression with Differential Regularization, Journal of Statistical Computation and Simulation_ |
| `quantile`                                    | Il modello implementa una reqressione spaziale quantilica, risolto tramite applicazione di FPIRLS, come descritto in _De Sanctis, M., Di Battista, I., Spatial Quantile regression with Partial Differential Equation Regularization, PACS report_. **Al momento non supportato** |


### SRTPDE

Ancora non disponbile.

### Decorare il fitting di un modello

Una riga come
```R
model <- SRPDE(y ~ f, penalty = PDE(Lf, u), data = data_frame, lambda = 1e-6)
```
descrive il modello astraendolo da eventuali dettagli computazionali/algoritmici. In particolare, stiamo immaginando di fissare il funzionale
```math
\sum_{i=1}^n (y_i - f(\boldsymbol{p}_i))^2 + \lambda_{\mathcal{D}} \int_{\mathcal{D}} (Lf-u)^2 d\mathcal{D}
```
Possiamo decorare il modo in cui risolviamo questo modello tramite il metodo `fit()`, ad esempio:
```R
library(fdaPDE2)
data("unit_square", package = "fdaPDE2")
unit_square <- Mesh(unit_square)

data_frame <- ## obtain your data in some way...
model <- SRPDE(y ~ f, domain = unit_square, data = data_frame) ## do not fix any lambda...

## ... because you might want to select it via GCV minimization!
lambda_grid = 10^seq(-4, -3, by = 0.1)
model$fit(
   lambda = gcv(optimizer = "grid", lambda = lambda_grid)
)
```
Si veda [GCV - Generalized Cross Validation](#gcv---generalized-cross-validation) per ulteriori dettagli sull'uso di `gcv()`. 

Il metodo `fit()`, indipendentemente dal modello considerato, possiede sempre degli opportuni default:
```R
library(fdaPDE2)
data("unit_square", package = "fdaPDE2")
unit_square <- Mesh(unit_square)

data_frame <- ## obtain your data in some way...
model <- SRPDE(y ~ f, domain = unit_square, data = data_frame) ## do not fix any lambda
## fit will default to some strategy to select the smoothing parameter
model$fit()
```

L'idea si estende a qualunque altra tipologia di modello, ad esempio, metodi basati su FPIRLS potranno configurare i parametri dell'algoritmo come segue:
```R
library(fdaPDE2)
data("unit_square", package = "fdaPDE2")
unit_square <- Mesh(unit_square)

data_frame <- ## obtain your data in some way...
model <- SRPDE(y ~ f, domain = unit_square, data = data_frame, family = "poisson")

lambda_grid = 10^seq(-4, -3, by = 0.1)
model$fit(
   lambda = gcv(optimizer = "grid", lambda = lambda_grid),
   fprils_params = list( ## customize FPIRLS execution
      max_iterations = 100,
      tolerance = 1e-6
   )
)

## you are not forced to set all the parameters, if you do not set some (or all) of them,
## there will always be a default
```

### GCV - Generalized Cross Validation
Questa sezione mostra l'API offerta per la selezione del parametro di smoothing in un modello di regressione mediante GCV.
```R
## stochastic approximation of degrees of freedom, fast but approximate
model$fit(
   lambda = gcv(
      edf_computation = "stochastic",
      seed = 143547654,    ## set seed in internal stochastic engine (for reproducibility)
      n_mc_samples = 500,  ## number of mc samples in stochastic algorithm
      optimizer = "...",
      ... ## additional arguments forwarded to the optimizer
)

## exact approximation of degrees of freedom, slow but exact
model$fit(
   lambda = gcv(
      edf_computation = "exact",
      optimizer = "...",
      ... ## additional arguments forwarded to the optimizer
)
```

Faccio notare che, in caso `edf_computation = "stochastic"`:
* `seed = NULL` (default) porta ad una inizializzazione casuale dell'algoritmo
* `n_mc_samples = NULL` (default) setta il numero di realizzazioni monte carlo a 100

Il parametro `optimizer` può essere uno dei seguenti ottimizzatori. La colonna note riporta dettagli aggiuntivi che decorano il comportamento dello specifico algoritmo di ottimizzazione:

| optimizer | note |
|-----------|------|
| `grid`    | Ottimizzazione su una griglia fissata di valori. Attualmente l'unica possibilità per modelli di regressione spazio-tempo. Richiede <ul><li> `lambda`: vettore di parametri da esplorare. </li></ul>     |
| `newton` | Metodo di Newton. Possibili parametri <ul> <li> `lambda`: punto iniziale del metodo iterativo </li> <li> `exact_derivative` : un flag booleano TRUE/FALSE per abilitare o meno il calcolo approssimato della derivata del GCV. (`exact_derivative = TRUE` attualmente non supportato) </li> <li> `max_iter`: numero massimo di iterazioni prima dell'arresto forzato </li> <li> `tolerance` : tolleranza sull'errore per avere l'arresto del metodo </li> <li> `step` : valore dello step utilizzato nel passo di Newton (passi adattivi `backtracking`, `wolfe` attualmente non supportati) </li> </ul> |
| `gradient_descent` | **mai testato nel contesto del GCV**. Gradient Descent. Possibili parametri <ul> <li> `lambda`: punto iniziale del metodo iterativo </li> <li> `exact_derivative` : un flag booleano TRUE/FALSE per abilitare o meno il calcolo approssimato della derivata del GCV. (`exact_derivative = TRUE` attualmente non supportato) </li> <li> `max_iter`: numero massimo di iterazioni prima dell'arresto forzato </li> <li> `tolerance` : tolleranza sull'errore per avere l'arresto del metodo </li> <li> `step` : valore dello step utilizzato nel passo di Newton (passi adattivi `backtracking`, `wolfe` attualmente non supportati) </li> </ul> |
| `bfgs` | **mai testato nel contesto del GCV**. Metodo di Broyden per ottimizzazione. Possibili parametri <ul> <li> `lambda`: punto iniziale del metodo iterativo </li> <li> `exact_derivative` : un flag booleano TRUE/FALSE per abilitare o meno il calcolo approssimato della derivata del GCV. (`exact_derivative = TRUE` attualmente non supportato) </li> <li> `max_iter`: numero massimo di iterazioni prima dell'arresto forzato </li> <li> `tolerance` : tolleranza sull'errore per avere l'arresto del metodo </li> <li> `step` : valore dello step utilizzato nel passo di Newton (passi adattivi `backtracking`, `wolfe` attualmente non supportati) </li> </ul> |
