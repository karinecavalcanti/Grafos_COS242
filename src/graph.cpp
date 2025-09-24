#include "Graph.h"
#include <fstream>
#include <iostream>
#include <queue>
#include <stack>
#include <iomanip> // Para formatar a saida de numeros decimais
#include <random> // Para geracao de numeros aleatorios
#include <limits> // Para numeric_limits

// Construtor
Graph::Graph(GraphRepresentation rep) : numVertices(0), numEdges(0), representationType(rep) {}

// Adiciona uma aresta entre u e v
void Graph::addEdge(int u, int v) {
    if (u < 1 || u > numVertices || v < 1 || v > numVertices) {
        std::cerr << "Erro: Tentativa de adicionar aresta com vertice invalido (" << u << ", " << v << ")." << std::endl;
        return;
    }

    // Ajusta para 0-based index para as estruturas internas
    --u; 
    --v;

    if (representationType == ADJACENCY_MATRIX) {
        if (!adjMatrix[u][v]) { // So adiciona se nao existir
            adjMatrix[u][v] = true;
            adjMatrix[v][u] = true; // Grafo nao-direcionado
            numEdges++;
        }
    } else { // ADJACENCY_LIST
        bool edgeExists = false;
        for (int neighbor : adjList[u]) {
            if (neighbor == v) {
                edgeExists = true;
                break;
            }
        }
        if (!edgeExists) {
            adjList[u].push_back(v);
            adjList[v].push_back(u); // Grafo nao-direcionado
            numEdges++;
        }
    }
}

// Leitura do grafo a partir de um arquivo
void Graph::readFromFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        return;
    }
    inputFile >> numVertices;

    // Inicializa as estruturas de dados com o numero correto de vertices
    if (representationType == ADJACENCY_MATRIX) {
        adjMatrix.assign(numVertices, std::vector<bool>(numVertices, false));
    } else { // ADJACENCY_LIST
        adjList.assign(numVertices, std::list<int>());
    }

    int u, v;
    while (inputFile >> u >> v) {
        addEdge(u, v);
    }

    inputFile.close();
}

// Retorna o numero de vertices
int Graph::getNumVertices() const {
    return numVertices;
}

// Retorna o numero de arestas
int Graph::getNumEdges() const {
    return numEdges;
}

// Implementacao dos graus
std::vector<int> Graph::getDegrees() {
    std::vector<int> degrees(numVertices, 0);
    if (representationType == ADJACENCY_MATRIX) {
        for (int i = 0; i < numVertices; ++i) {
            for (int j = 0; j < numVertices; ++j) {
                if (adjMatrix[i][j]) {
                    degrees[i]++;
                }
            }
        }
    } else { // ADJACENCY_LIST
        for (int i = 0; i < numVertices; ++i) {
            degrees[i] = adjList[i].size();
        }
    }
    return degrees;
}

int Graph::getMinDegree() {
    std::vector<int> degrees = getDegrees();
    if (degrees.empty()) return 0;
    return *std::min_element(degrees.begin(), degrees.end());
}

int Graph::getMaxDegree() {
    std::vector<int> degrees = getDegrees();
    if (degrees.empty()) return 0;
    return *std::max_element(degrees.begin(), degrees.end());
}

double Graph::getAverageDegree() {
    if (numVertices == 0) return 0.0;
    return (double)(2 * numEdges) / numVertices;
}

double Graph::getMedianDegree() {
    std::vector<int> degrees = getDegrees();
    if (degrees.empty()) return 0.0;
    std::sort(degrees.begin(), degrees.end());
    int mid = degrees.size() / 2;
    if (degrees.size() % 2 == 0) {
        return (double)(degrees[mid - 1] + degrees[mid]) / 2.0;
    } else {
        return (double)degrees[mid];
    }
}

// Escreve informacoes do grafo em um arquivo
void Graph::writeGraphInfo(const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo de saida: " << filename << std::endl;
        return;
    }

    outputFile << "Numero de Vertices: " << numVertices << std::endl;
    outputFile << "Numero de Arestas: " << numEdges << std::endl;
    outputFile << "Grau Minimo: " << getMinDegree() << std::endl;
    outputFile << "Grau Maximo: " << getMaxDegree() << std::endl;
    outputFile << std::fixed << std::setprecision(2); // Formata para 2 casas decimais
    outputFile << "Grau Medio: " << getAverageDegree() << std::endl;
    outputFile << "Mediana do Grau: " << getMedianDegree() << std::endl;
    outputFile << std::endl;

    // Informacoes sobre Componentes Conexas
    std::vector<std::vector<int>> components = getConnectedComponents();
    outputFile << "Numero de Componentes Conexas: " << components.size() << std::endl;
    outputFile << "Componentes Conexas (em ordem decrescente de tamanho):" << std::endl;

    // Ordena as componentes por tamanho decrescente
    std::sort(components.begin(), components.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        return a.size() > b.size();
    });

    for (size_t i = 0; i < components.size(); ++i) {
        outputFile << "  Componente " << (i + 1) << " (Tamanho: " << components[i].size() << "): ";
        for (size_t j = 0; j < components[i].size(); ++j) {
            outputFile << (components[i][j] + 1) << (j == components[i].size() - 1 ? "" : ", ");
        }
        outputFile << std::endl;
    }

    outputFile.close();
}

// Implementacao da Busca em Largura (BFS)
std::map<int, int> Graph::BFS(int startVertex, const std::string& outputFilename) {
    std::map<int, int> parent; // Armazena o pai de cada vertice na arvore BFS
    std::map<int, int> level;  // Armazena o nivel de cada vertice na arvore BFS
    std::vector<bool> visited(numVertices, false);
    std::queue<int> q;

    // Ajusta para 0-based index
    int startIdx = startVertex - 1;

    visited[startIdx] = true;
    level[startVertex] = 0; // O nivel da raiz e 0
    q.push(startIdx);
    
    // Para o arquivo de saida
    std::ofstream outputFile(outputFilename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo de saida para BFS: " << outputFilename << std::endl;
        return parent;
    }
    outputFile << "Arvore BFS a partir do vertice " << startVertex << ":" << std::endl;
    outputFile << "Vertice (Pai, Nivel)" << std::endl;
    outputFile << startVertex << " (N/A, 0)" << std::endl; // Vertice inicial nao tem pai

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        auto processNeighbor = [&](int v) {
            if (!visited[v]) {
                visited[v] = true;
                parent[v + 1] = u + 1; // Armazena pai (1-based)
                level[v + 1] = level[u + 1] + 1; // Calcula nivel (1-based)
                outputFile << (v + 1) << " (" << (u + 1) << ", " << level[v + 1] << ")" << std::endl;
                q.push(v);
            }
        };

        if (representationType == ADJACENCY_MATRIX) {
            for (int v = 0; v < numVertices; ++v) {
                if (adjMatrix[u][v]) {
                    processNeighbor(v);
                }
            }
        } else { // ADJACENCY_LIST
            for (int v : adjList[u]) {
                processNeighbor(v);
            }
        }
    }
    outputFile.close();
    return parent;
}

// Implementacao da Busca em Profundidade (DFS)
std::map<int, int> Graph::DFS(int startVertex, const std::string& outputFilename) {
    std::map<int, int> parent; // Armazena o pai de cada vertice na arvore DFS
    std::map<int, int> level; // Armazena o nivel de cada vertice na arvore DFS
    std::vector<bool> visited(numVertices, false);
    std::stack<int> s;

    // Ajusta para 0-based index
    int startIdx = startVertex - 1;

    s.push(startIdx);
    visited[startIdx] = true;
    level[startVertex] = 0;

    // Para o arquivo de saida
    std::ofstream outputFile(outputFilename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo de saida para DFS: " << outputFilename << std::endl;
        return parent;
    }
    outputFile << "Arvore DFS a partir do vertice " << startVertex << ":" << std::endl;
    outputFile << "Vertice (Pai, Nivel)" << std::endl;
    outputFile << startVertex << " (N/A, 0)" << std::endl; // Vertice inicial nao tem pai

    while (!s.empty()) {
        int u = s.top();
        s.pop();

        // Para garantir que os vizinhos sejam processados em ordem para a saida
        // e melhor iterar e adicionar na pilha na ordem inversa ou usar um conjunto
        // A lista de adjacencia pode nao estar ordenada, entao vamos coletar e ordenar
        std::vector<int> neighbors;
        if (representationType == ADJACENCY_MATRIX) {
            for (int v = 0; v < numVertices; ++v) {
                if (adjMatrix[u][v]) {
                    neighbors.push_back(v);
                }
            }
        } else {
            for (int v : adjList[u]) {
                neighbors.push_back(v);
            }
        }
        std::sort(neighbors.rbegin(), neighbors.rend()); // Para empilhar e desempilhar em uma ordem previsivel (menor primeiro)

        for (int v : neighbors) {
            if (!visited[v]) {
                visited[v] = true;
                parent[v + 1] = u + 1;
                level[v + 1] = level[u + 1] + 1;
                outputFile << (v + 1) << " (" << (u + 1) << ", " << level[v + 1] << ")" << std::endl;
                s.push(v);
            }
        }
    }
    outputFile.close();
    return parent;
}


// Implementacao da Distancia entre dois vertices (usando BFS)
int Graph::getDistance(int u, int v) {
    // Ajusta para 0-based index
    int startIdx = u - 1;
    int targetIdx = v - 1;

    if (startIdx < 0 || startIdx >= numVertices || targetIdx < 0 || targetIdx >= numVertices) {
        std::cerr << "Erro: Vertices invalidos para calculo de distancia (" << u << ", " << v << ")." << std::endl;
        return -1; // Retorna -1 para indicar que a distancia nao foi encontrada ou vertices invalidos
    }

    std::vector<int> dist(numVertices, -1); // dist[i] armazena a distancia de startIdx ate i
    std::queue<int> q;

    dist[startIdx] = 0;
    q.push(startIdx);

    while (!q.empty()) {
        int curr = q.front();
        q.pop();

        if (curr == targetIdx) {
            return dist[curr];
        }

        auto processNeighbor = [&](int neighbor) {
            if (dist[neighbor] == -1) { // Se nao foi visitado
                dist[neighbor] = dist[curr] + 1;
                q.push(neighbor);
            }
        };

        if (representationType == ADJACENCY_MATRIX) {
            for (int i = 0; i < numVertices; ++i) {
                if (adjMatrix[curr][i]) {
                    processNeighbor(i);
                }
            }
        } else { // ADJACENCY_LIST
            for (int neighbor : adjList[curr]) {
                processNeighbor(neighbor);
            }
        }
    }
    return -1; // Nao ha caminho entre u e v
}

// Metodos auxiliares para Componentes Conexas (BFS ou DFS)
void Graph::BFSUtil(int startIdx, std::vector<bool>& visited, std::vector<int>& componentVertices) {
    std::queue<int> q;
    q.push(startIdx);
    visited[startIdx] = true;
    componentVertices.push_back(startIdx);

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        auto processNeighbor = [&](int v) {
            if (!visited[v]) {
                visited[v] = true;
                q.push(v);
                componentVertices.push_back(v);
            }
        };

        if (representationType == ADJACENCY_MATRIX) {
            for (int v = 0; v < numVertices; ++v) {
                if (adjMatrix[u][v]) {
                    processNeighbor(v);
                }
            }
        } else { // ADJACENCY_LIST
            for (int v : adjList[u]) {
                processNeighbor(v);
            }
        }
    }
}

void Graph::DFSUtil(int u, std::vector<bool>& visited, std::vector<int>& componentVertices) {
    visited[u] = true;
    componentVertices.push_back(u);

    auto processNeighbor = [&](int v) {
        if (!visited[v]) {
            DFSUtil(v, visited, componentVertices);
        }
    };

    if (representationType == ADJACENCY_MATRIX) {
        for (int v = 0; v < numVertices; ++v) {
            if (adjMatrix[u][v]) {
                processNeighbor(v);
            }
        }
    } else { // ADJACENCY_LIST
        for (int v : adjList[u]) {
            processNeighbor(v);
        }
    }
}

std::vector<int> Graph::BFS_search(int startVertex) {
    std::vector<int> parent(numVertices, -1);
    std::vector<bool> visited(numVertices, false);
    std::queue<int> q;

    int s = startVertex - 1;
    if (s < 0 || s >= numVertices) return parent; // Proteção adicional
    
    visited[s] = true;
    q.push(s);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        
        // Processar vizinhos conforme a representação
        if (representationType == ADJACENCY_MATRIX) {
            for (int v = 0; v < numVertices; ++v) {
                if (adjMatrix[u][v] && !visited[v]) {
                    visited[v] = true;
                    parent[v] = u;
                    q.push(v);
                }
            }
        } else { // ADJACENCY_LIST
            for (int v : adjList[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    parent[v] = u;
                    q.push(v);
                }
            }
        }
    }
    return parent;
}

std::vector<int> Graph::DFS_search(int startVertex) {
    std::vector<int> parent(numVertices, -1);
    std::vector<bool> visited(numVertices, false);
    std::stack<int> st;

    int s = startVertex - 1;
    if (s < 0 || s >= numVertices) return parent; // Proteção adicional
    
    visited[s] = true;
    st.push(s);

    while (!st.empty()) {
        int u = st.top(); st.pop();
        
        // Processar vizinhos conforme a representação
        if (representationType == ADJACENCY_MATRIX) {
            for (int v = 0; v < numVertices; ++v) {
                if (adjMatrix[u][v] && !visited[v]) {
                    visited[v] = true;
                    parent[v] = u;
                    st.push(v);
                }
            }
        } else { // ADJACENCY_LIST
            for (int v : adjList[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    parent[v] = u;
                    st.push(v);
                }
            }
        }
    }
    return parent;
}


// Implementacao de Componentes Conexas
std::vector<std::vector<int>> Graph::getConnectedComponents() {
    std::vector<std::vector<int>> components;
    std::vector<bool> visited(numVertices, false);

    for (int i = 0; i < numVertices; ++i) {
        if (!visited[i]) {
            std::vector<int> currentComponent;
            // Escolha entre BFSUtil ou DFSUtil
            BFSUtil(i, visited, currentComponent); 
            // DFSUtil(i, visited, currentComponent); 
            components.push_back(currentComponent);
        }
    }
    return components;
}

// Nova funcao: Executa uma BFS completa e retorna a maior distancia encontrada a partir de startVertex
int Graph::getMaxDistance_from_BFS(int startVertex) {
    // startVertex e 0-based
    if (startVertex < 0 || startVertex >= numVertices) {
        // Lidar com erro ou retornar 0 para vertice invalido
        return 0;
    }

    std::vector<int> dist(numVertices, -1);
    std::queue<int> q;

    dist[startVertex] = 0;
    q.push(startVertex);
    int currentMaxDistance = 0;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        currentMaxDistance = std::max(currentMaxDistance, dist[u]);

        auto processNeighbor = [&](int v_neighbor) { // Renomeado para evitar conflito com 'v' em loops anteriores
            if (dist[v_neighbor] == -1) {
                dist[v_neighbor] = dist[u] + 1;
                q.push(v_neighbor);
            }
        };

        if (representationType == ADJACENCY_MATRIX) {
            for (int v_idx = 0; v_idx < numVertices; ++v_idx) {
                if (adjMatrix[u][v_idx]) {
                    processNeighbor(v_idx);
                }
            }
        } else { // ADJACENCY_LIST
            for (int v_neighbor : adjList[u]) {
                processNeighbor(v_neighbor);
            }
        }
    }
    return currentMaxDistance;
}

int Graph::getDiameter() {
    // Função auxiliar BFS que funciona para AMBAS as representações
    auto bfs = [&](int start) {
        std::vector<int> dist(numVertices, -1);
        std::queue<int> q;
        dist[start] = 0;
        q.push(start);

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            // Processar vizinhos conforme a representação escolhida
            if (representationType == ADJACENCY_MATRIX) {
                for (int v = 0; v < numVertices; ++v) {
                    if (adjMatrix[u][v] && dist[v] == -1) {
                        dist[v] = dist[u] + 1;
                        q.push(v);
                    }
                }
            } else { // ADJACENCY_LIST
                for (int v : adjList[u]) {
                    if (dist[v] == -1) {
                        dist[v] = dist[u] + 1;
                        q.push(v);
                    }
                }
            }
        }
        return dist;
    };

    // Se grafo vazio
    if (numVertices == 0) return 0;

    auto components = getConnectedComponents();

    // Caso o grafo seja conexo → calcula diâmetro exato
    if (components.size() == 1) {
        int overallMaxDistance = 0;
        for (int i = 0; i < numVertices; ++i) {
            std::vector<int> dist = bfs(i);
            int localMax = *std::max_element(dist.begin(), dist.end());
            overallMaxDistance = std::max(overallMaxDistance, localMax);
        }
        return overallMaxDistance;
    }

    // Caso grafo desconexo → faz aproximação
    int numSamples = std::min(10, (int)components.size());
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> compDist(0, (int)components.size() - 1);

    int approxDiameter = 0;
    for (int i = 0; i < numSamples; ++i) {
        const auto& comp = components[compDist(gen)];
        if (comp.empty()) continue;

        std::uniform_int_distribution<> vertexDist(0, (int)comp.size() - 1);
        int start = comp[vertexDist(gen)];

        std::vector<int> dist = bfs(start);
        int localMax = 0;
        for (int v : comp) {
            if (dist[v] > localMax) localMax = dist[v];
        }
        approxDiameter = std::max(approxDiameter, localMax);
    }

    return approxDiameter;
}
