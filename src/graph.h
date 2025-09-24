#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <list>
#include <string>
#include <map>
#include <algorithm> // Para std::sort

// Enumeracao para o tipo de representacao do grafo
enum GraphRepresentation {
    ADJACENCY_MATRIX,
    ADJACENCY_LIST
};

class Graph {
private:
    int numVertices;
    int numEdges;
    GraphRepresentation representationType;

    // Representacao por Matriz de Adjacencia
    std::vector<std::vector<bool>> adjMatrix;

    // Representacao por Lista de Adjacencia
    std::vector<std::list<int>> adjList;

    // Metodos auxiliares para BFS e DFS (para componentes conexas e outros)
    void BFSUtil(int startVertex, std::vector<bool>& visited, std::vector<int>& componentVertices);
    void DFSUtil(int startVertex, std::vector<bool>& visited, std::vector<int>& componentVertices);

public:
    // Construtor
    Graph(GraphRepresentation rep);

    // Metodos de Leitura e Escrita
    void readFromFile(const std::string& filename);
    void writeGraphInfo(const std::string& filename);

    // Metodos para manipulacao do grafo
    void addEdge(int u, int v); // Adiciona uma aresta entre u e v
    int getNumVertices() const;
    int getNumEdges() const;
    
    // Metodos para informacoes do grafo
    std::vector<int> getDegrees();
    int getMinDegree();
    int getMaxDegree();
    double getAverageDegree();
    double getMedianDegree();

    // Algoritmos de Busca
    std::map<int, int> BFS(int startVertex, const std::string& outputFilename); // Retorna pais
    std::map<int, int> DFS(int startVertex, const std::string& outputFilename); // Retorna pais

    // Metodos para Distancias e Diametro
    int getDistance(int u, int v);
    int getDiameter();
    int getApproximateDiameter(int numSamples);
    int getMaxDistance_from_BFS(int startVertex);
    
    // Metodos para Componentes Conexas
    std::vector<std::vector<int>> getConnectedComponents();

    // BFS sem escrita em arquivo, retorna apenas o vetor de pais
    std::vector<int> BFS_search(int startVertex);
    std::vector<int> DFS_search(int startVertex);

};

#endif // GRAPH_H