#include "Graph.h"
#include <iostream>
#include <chrono>
#include <random> 
#include <algorithm> 

int main() {
    // Escolha o tipo de representacao
    GraphRepresentation repType = ADJACENCY_MATRIX; // Ou ADJACENCY_MATRIX

    Graph graph(repType);

    // 1. Leitura do grafo
    std::string inputFilename = "grafo_1.txt"; // Nome do arquivo de entrada
    graph.readFromFile(inputFilename);
    std::cout << "\nGrafo lido do arquivo: " << inputFilename << std::endl;

    // ... após a leitura do grafo, PAUSE aqui
    std::cout << "\nPressione Enter para continuar e medir memória..." << std::endl;
    std::cin.get(); // Pausa a execução

    if (graph.getNumVertices() == 0) {
        std::cout << "Nenhum grafo lido ou erro no arquivo de entrada." << std::endl;
        return 1;
    }

    std::cout << "Grafo lido com " << graph.getNumVertices() << " vertices e " << graph.getNumEdges() << " arestas." << std::endl;
    
    // 2. Escrita das informacoes do grafo
    std::string outputInfoFilename = "graph_info.txt";
    graph.writeGraphInfo(outputInfoFilename);
    std::cout << "Informacoes do grafo escritas em: " << outputInfoFilename << std::endl;

    // 3. Busca em Largura (BFS)
    int bfsStartVertex = 1; // Vertice inicial para BFS (1-based)
    std::string bfsOutputFilename = "bfs_tree_output.txt";
    std::cout << "\nExecutando BFS a partir do vertice " << bfsStartVertex << "..." << std::endl;
    graph.BFS(bfsStartVertex, bfsOutputFilename);
    std::cout << "Arvore BFS escrita em: " << bfsOutputFilename << std::endl;

    // 4. Busca em Profundidade (DFS)
    int dfsStartVertex = 1; // Vertice inicial para DFS (1-based)
    std::string dfsOutputFilename = "dfs_tree_output.txt";
    std::cout << "\nExecutando DFS a partir do vertice " << dfsStartVertex << "..." << std::endl;
    graph.DFS(dfsStartVertex, dfsOutputFilename);
    std::cout << "Arvore DFS escrita em: " << dfsOutputFilename << std::endl;

    // 5. Distancia entre vertices
    int u_dist = 1;
    int v_dist = 5;
    int distance = graph.getDistance(u_dist, v_dist);
    if (distance != -1) {
        std::cout << "\nDistancia entre " << u_dist << " e " << v_dist << ": " << distance << std::endl;
    } else {
        std::cout << "\n\nNao ha caminho entre " << u_dist << " e " << v_dist << "." << std::endl;
    }

    // 6. Diametro do grafo
    std::cout << "\nCalculando diametro do grafo" << std::endl;
    auto start_approx_time = std::chrono::high_resolution_clock::now();
    int diameter = graph.getDiameter();
    auto end_approx_time = std::chrono::high_resolution_clock::now();
    auto approx_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_approx_time - start_approx_time);
    std::cout << "Diametro do grafo: " << diameter << std::endl;
    std::cout << "Tempo para calcular diametro aproximado: " << approx_duration.count() << " ms." << std::endl;

    //Medicao de tempo para BFS (para a parte de estudos de caso)
    std::cout << "\nMedindo tempo para 100 BFSs aleatorias..." << std::endl;
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            int v = (i % graph.getNumVertices()) + 1;
            graph.BFS_search(v);      // <<< usa versão rápida
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Tempo medio para 100 BFSs: " 
                << (double)duration.count() / 100.0 << " ms.\n";
    }

    //Medicao de tempo para DFS (para a parte de estudos de caso)
    std::cout << "\nMedindo tempo para 100 DFSs aleatorias..." << std::endl;
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            int v = (i % graph.getNumVertices()) + 1;
            graph.DFS_search(v);      // <<< usa versão rápida
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Tempo medio para 100 DFSs: " 
                << (double)duration.count() / 100.0 << " ms.\n";
    }
    
    //Pai dos v´ertices 10, 20, 30 na ´arvore geradora induzida pela BFS quando iniciamos a busca nos v´ertices 1, 2, 3.
    std::vector<int> starts = {1, 2, 3};
    std::vector<int> targets = {10, 20, 30};

    for (int s : starts) {
        std::cout << "\n=== BFS iniciando em " << s << " ===\n";
        auto parentBFS = graph.BFS(s, "bfs_temp.txt");

        for (int t : targets) {
            auto it = parentBFS.find(t);
            if (it != parentBFS.end()) {
                std::cout << "Pai de " << t << " na BFS = " << it->second << "\n";
            } else {
                std::cout << "Vertice " << t << " nao alcancavel a partir de " << s << " (BFS)\n";
            }
        }
    }

    //Pai dos v´ertices 10, 20, 30 na ´arvore geradora induzida pela DFS quando iniciamos a busca nos v´ertices 1, 2, 3.
    for (int s : starts) {
        std::cout << "\n=== DFS iniciando em " << s << " ===\n";
        auto parentDFS = graph.DFS(s, "dfs_temp.txt");

        for (int t : targets) {
            auto it = parentDFS.find(t);
            if (it != parentDFS.end()) {
                std::cout << "Pai de " << t << " na DFS = " << it->second << "\n";
            } else {
                std::cout << "Vertice " << t << " nao alcancavel a partir de " << s << " (DFS)\n";
            }
        }
    }

    //Distˆancia entre os pares de v´ertices (10,20), (10,30), (20,30).
    std::vector<std::pair<int,int>> pares = {
        {10,20}, {10,30}, {20,30}
    };
    for (auto [a,b] : pares) {
        int d = graph.getDistance(a, b);
        if (d != -1)
            std::cout << "\nDistancia entre " << a << " e " << b << " = " << d << std::endl;
        else
            std::cout << "Nao ha caminho entre " << a << " e " << b << std::endl;
    }

    //Componentes conexas
    auto comps = graph.getConnectedComponents();
    std::cout << "\nNumero de componentes conexas: " << comps.size() << "\n";
    size_t menor = comps[0].size();
    size_t maior = comps[0].size();
    for (const auto& c : comps) {
        menor = std::min(menor, c.size());
        maior = std::max(maior, c.size());
    }
    std::cout << "Tamanho da menor componente: " << menor << "\n";
    std::cout << "Tamanho da maior componente: " << maior << "\n";

    return 0;
}