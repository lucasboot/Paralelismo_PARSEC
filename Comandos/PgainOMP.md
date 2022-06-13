### Discussão sobre o porte da função pgain do Streamcluster para OpenMP

- #pragma omp parallel for private() shared()
Região de execução paralela de um for com variáveis que não podem ser compartilhadas pelos processos e as que devem ser de acesso compartilhado.

- #pragma omp flush
Essa diretiva aponta uma região para o compilador garantir que todas as threads em uma região paralela tem a mesma visão de um objeto específico na memória.

