### Discussão sobre o porte da função pgain do Streamcluster para OpenMP

- #pragma omp parallel for private() shared()
Região de execução paralela de um for com variáveis que não podem ser compartilhadas pelos processos e as que devem ser de acesso compartilhado.

- #pragma omp flush
Essa diretiva aponta uma região para o compilador garantir que todas as threads em uma região paralela tem a mesma visão de um objeto específico na memória.

- #pragma omp single
Apenas uma thread executará determinada trecho de código, não necessariamente a thread 0

- #pragma omp sections + #pragma omp section
Determina seções de código que serão executadas por threads diferentes de forma paralela

- #pragma omp critical
Determina uma região crítica do código, na qual o acesso à uma região da memória só pode ser feito por 1 thread por vez

- #pragma omp atomic
Essa diretiva garante que as condições de corrida serão evitadas a partir de um controle direto da concorrência entre threads que podem ler ou escrever em uma determinada região da memória. Ela é interessante por garantir a eficiência de algoritmos concorrentes com menos locks.

#### Atomic vs. Critical
- O efeito de ambas seriam o mesmo em uma região como:
#pragma omp atomic
    ++number_of_centers_to_close;

An OpenMP critical section is completely general - it can surround any arbitrary block of code. You pay for that generality, however, by incurring significant overhead every time a thread enters and exits the critical section (on top of the inherent cost of serialization).

(In addition, in OpenMP all unnamed critical sections are considered identical (if you prefer, there's only one lock for all unnamed critical sections), so that if one thread is in one [unnamed] critical section as above, no thread can enter any [unnamed] critical section. As you might guess, you can get around this by using named critical sections).

An atomic operation has much lower overhead. Where available, it takes advantage on the hardware providing (say) an atomic increment operation; in that case there's no lock/unlock needed on entering/exiting the line of code, it just does the atomic increment which the hardware tells you can't be interfered with.

The upsides are that the overhead is much lower, and one thread being in an atomic operation doesn't block any (different) atomic operations about to happen. The downside is the restricted set of operations that atomic supports.