#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>

using namespace std;

vector<double> solveEquations(vector<vector<double>>& matrix, vector<double>& treeCurrents) {
    size_t n = matrix.size();

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; i++) {
        double maxElem = abs(matrix[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix[k][i]) > maxElem) {
                maxElem = abs(matrix[k][i]);
                maxRow = k;
            }
        }
        for (int k = i; k < n + 1; k++) {
            double tmp = matrix[maxRow][k];
            matrix[maxRow][k] = matrix[i][k];
            matrix[i][k] = tmp;
        }

        for (int j = i + 1; j < n; j++) {
            double factor = matrix[j][i] / matrix[i][i];
            for (int k = i; k < n + 1; k++) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
    }

    // Обратный ход метода Гаусса
    vector<double> result(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += matrix[i][j] * result[j];
        }
        result[i] = (treeCurrents[i] - sum) / matrix[i][i];
    }
    return result;
}


struct Edge {
    int src, dest;
    double resistance;
    double emf;
};

struct Graph {
    int V, E;
    vector<Edge> edges;
};

// Создание графа с V вершинами и E ребрами
Graph createGraph(int V, int E) {
    Graph graph;
    graph.V = V;
    graph.E = E;
    return graph;
}

void addEdge(Graph& graph, int src, int dest, double resistance, double emf = 0) {
    Edge edge = {src, dest, resistance, emf};
    graph.edges.push_back(edge);
}

int findRoot(vector<int> &parent, int i) {
    while (parent[i] != i) {
        i = parent[i];
    }
    return i;
}

void unionTrees(vector<int> &parent, int x, int y) {
    int rootX = findRoot(parent, x);
    int rootY = findRoot(parent, y);
    parent[rootY] = rootX;
}

bool compareEdges(Edge a, Edge b) {
    return a.resistance < b.resistance;
}

// Функция для поиска максимального остовного дерева
vector<Edge> mstKruskal(Graph graph) {
    vector<Edge> mst;
    vector<int> parent(graph.V);
    // Инициализируем вектор родительских вершин так, чтобы каждая вершина была своим родителем
    for (int i = 0; i < graph.V; i++) {
        parent[i] = i;
    }
    // Сортируем ребра по возрастанию веса
    sort(graph.edges.begin(), graph.edges.end(), compareEdges);
    // Проходимся по отсортированным ребрам и добавляем их в остовное дерево, если они не создают цикл
    for (int i = 0; i < graph.E; i++) {
        Edge currentEdge = graph.edges[i];
        int rootSrc = findRoot(parent, currentEdge.src);
        int rootDest = findRoot(parent, currentEdge.dest);
        if (rootSrc != rootDest) {
            mst.push_back(currentEdge);
            unionTrees(parent, rootSrc, rootDest);
        }
    }
    return mst;
}


void solveCircuit(Graph& graph) {
    int V = graph.V;
    int E = graph.E;

    vector<Edge> mstEdges = mstKruskal(graph);

    // Список смежности для остовного дерева
    vector<vector<pair<int, double>>> adjList(V);
    for (Edge &edge: mstEdges) {
        adjList[edge.src].emplace_back(edge.dest, edge.resistance);
        adjList[edge.dest].emplace_back(edge.src, edge.resistance);
    }

    // Токи в остовном дереве
    vector<double> treeCurrents(V, 0);
    for (int i = 0; i < V; i++) {
        for (auto &neighbor: adjList[i]) {
            int dest = neighbor.first;
            double resistance = neighbor.second;
            treeCurrents[i] += resistance * (i < dest ? 1 : -1); // Используем знак, чтобы учитывать направление тока
        }
    }

    // Создаем СЛАУ
    vector<vector<double>> matrix(V, vector<double>(V, 0));
    for (int i = 0; i < V; i++) {
        for (pair<int, double> &neighbor: adjList[i]) {
            int dest = neighbor.first;
            double resistance = neighbor.second;
            if (i < dest) {
                matrix[i][i] += 1 / resistance;
                matrix[i][dest] -= 1 / resistance;
            } else {
                matrix[i][i] += 1 / resistance;
                matrix[i][dest] += 1 / resistance;
                matrix[i][i] -= 1 / resistance;
            }
        }
    }

    // Добавляем информацию об ЭДС
    for (int i = 0; i < E; i++) {
        Edge &edge = graph.edges[i];
        int src = edge.src;
        int dest = edge.dest;
        double emf = edge.emf;
        if (emf != 0) {
            matrix[src][dest] += 1;
            matrix[dest][src] -= 1;
            treeCurrents[src] -= emf;
            treeCurrents[dest] += emf;
        }
    }

    vector<double> currents = solveEquations(matrix, treeCurrents);

    for (int i = 0; i < E; i++) {
        Edge &edge = graph.edges[i];
        int src = edge.src;
        int dest = edge.dest;
        double current = (src < dest ? currents[src] - currents[dest] : currents[dest] - currents[src]);
        cout << "I" << i << " = " << current << " A" << endl;
    }
}

Graph readGraph(const char* filename) {
    try {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Не удалось открыть файл.");
        }
        int V, E;
        if (!(file >> V) || !(file >> E)) {
            throw std::runtime_error("Некорректный формат файла.");
        }
        if (V < 0 || E < 0) {
            throw std::runtime_error("V и E должны быть неотрицательными.");
        }
        Graph graph = createGraph(V, E);
        int u, v;
        double r, emf;
        while (file >> u >> v >> r >> emf) {
            addEdge(graph, u, v, r, emf);
        }
        file.close();
        return graph;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}



int main() {
#define N 1
    auto times = new long long [N];
    for (size_t i = 0; i < N; ++i) {
        auto graph = readGraph(R"(D:/mkt/graph10.txt)");

        auto start = std::chrono::high_resolution_clock::now();

        solveCircuit(graph);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        times[i] = duration.count();
    }
    long long sum = 0;
    for (int i = 0; i < N; ++i) {
        sum += times[i];
    }
    std::cout << "Время работы функции: " << sum / N << " нс" << std::endl;

    return 0;
}
